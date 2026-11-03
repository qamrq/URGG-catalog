#!/usr/bin/env bash
#                 ############################
#               ##                            ##
#             ##  ####                    ####  ##
#           ##   ######                  ######   ##
#         ##    ########                ########    ##
#       ##      ##########            ##########      ##
#     ##        ####  ####################  ####        ##
#   ##          #####  ##################  #####          ##
# ##            ################################            ##
#   ##            ############################            ##
#     ##              ####################              ##
#       ##             ##  ##########  ##             ##
#         ##           ##   ########   ##           ##
#           ##          #    ######    #          ##
#             ##              ####              ##
#               ##                            ##
#                 ############################
#
# ==========================================================
#  PROJECT : URGG catalog
#  SCRIPT  : 07_phabox2.sh  -  vOTU taxonomy, host and lifestyle (PhaBOX2)
#  AUTHOR  : QRQ
#  INPUT   : votu_representatives.fa
#            GTDB-Tk summaries
#            $URGG_DB/phabox2/phabox_db_v2_2
#  OUTPUT  : final_prediction_summary.tsv, vOTU_taxonomy_host.tsv
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"

source "$(dirname "$0")/config.sh"

# --- confirmed real env name ---
ENV_PHABOX2="${ENV_PHABOX2:-phabox2}"
PHABOX2_DB="${URGG_DB}/phabox2/phabox_db_v2_2"

activate_env "${ENV_PHABOX2}"

log "Step 7: PhaBOX2 taxonomy/host/lifestyle on vOTU representatives"

VOTU_FA="${VOTU_OUT}/votu_representatives.fa"
if [[ ! -s "${VOTU_FA}" ]]; then
    log "ERROR: ${VOTU_FA} not found. Run 06_votu_clustering.sh first."
    exit 1
fi
n_votus=$(grep -c '^>' "${VOTU_FA}")
log "Input: ${n_votus} vOTU representatives"

if [[ ! -d "${PHABOX2_DB}" ]]; then
    log "ERROR: PhaBOX2 database not found at ${PHABOX2_DB}"
    log "       Check the real path -- this was pulled from an earlier session's setup."
    exit 1
fi

GTDBTK_SUMMARY=$(find "${GTDBTK_OUTPUT_DIR}" -maxdepth 1 -iname "gtdbtk.all.summary.tsv" 2>/dev/null | head -1)

if [[ -z "${GTDBTK_SUMMARY}" ]]; then
    log "No pre-combined gtdbtk.all.summary.tsv found -- falling back to bac120+ar53 concat"
    GTDBTK_SUMMARY="${VOTU_OUT}/_gtdbtk_combined_summary.tsv"
    bac_summary=$(find "${GTDBTK_OUTPUT_DIR}" -iname "*bac120*summary.tsv" 2>/dev/null | head -1)
    ar_summary=$(find "${GTDBTK_OUTPUT_DIR}" -iname "*ar53*summary.tsv" 2>/dev/null | head -1)
    first=1
    : > "${GTDBTK_SUMMARY}"
    for f in "${bac_summary}" "${ar_summary}"; do
        [[ -f "${f}" ]] || continue
        if [[ "${first}" -eq 1 ]]; then cat "${f}" >> "${GTDBTK_SUMMARY}"; first=0
        else tail -n +2 "${f}" >> "${GTDBTK_SUMMARY}"; fi
    done
    [[ -s "${GTDBTK_SUMMARY}" ]] || GTDBTK_SUMMARY=""
fi

if [[ -n "${GTDBTK_SUMMARY}" ]]; then
    log "GTDB-Tk summary: ${GTDBTK_SUMMARY} ($(($(wc -l < "${GTDBTK_SUMMARY}") - 1)) genomes)"
    BGTDB_ARG=""
else
    log "WARNING: no GTDB-Tk summary found under ${GTDBTK_OUTPUT_DIR}"
    log "         Proceeding WITHOUT --bgtdb."
    BGTDB_ARG=""
fi

log "Running PhaBOX2 (this may take a while -- watch ${LOG_DIR}/phabox2.log)"
phabox2 --task end_to_end \
    --dbdir "${PHABOX2_DB}" \
    --outpth "${PHAGCN_OUT}" \
    --contigs "${VOTU_FA}" \
    --threads "${THREADS}" \
    --bfolder "${PHASE5_DIR}/00_mags_short_headers" \
    --magonly Y \
    ${BGTDB_ARG} \
    > "${LOG_DIR}/phabox2.log" 2>&1

RESULTS="${PHAGCN_OUT}/final_prediction/final_prediction_summary.tsv"
if [[ ! -s "${RESULTS}" ]]; then
    log "ERROR: PhaBOX2 produced no final_prediction_summary.tsv -- see ${LOG_DIR}/phabox2.log"
    exit 1
fi

n_results=$(($(wc -l < "${RESULTS}") - 1))
n_called_virus=$(awk -F'\t' 'NR>1 && $3=="virus"' "${RESULTS}" | wc -l)
log "PhaBOX2 finished: ${n_results} vOTUs processed, ${n_called_virus} confirmed as virus by PhaMer"
if [[ "${n_called_virus}" -lt $((n_results * 9 / 10)) ]]; then
    log "  NOTE: <90% PhaMer-confirmed on an already-CheckV-filtered set is unusual -- worth a look"
fi

log "Joining taxonomy/host/lifestyle onto vOTU_index.tsv..."
python3 << PYEOF
import csv

votu_idx = {}
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    rows = list(r)
    fieldnames = r.fieldnames

phabox = {}
with open("${RESULTS}") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        phabox[row["Accession"]] = row

extra_cols = ["Pred", "Lineage", "PhaGCNScore", "TYPE", "PhaTYPScore",
              "Host", "CHERRYScore", "Method", "Host_NCBI_lineage", "Host_GTDB_lineage"]

out_path = "${PHAGCN_OUT}/vOTU_taxonomy_host.tsv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(fieldnames + extra_cols)
    n_joined = 0
    for row in rows:
        rep = row["representative"]
        p = phabox.get(rep, {})
        if p:
            n_joined += 1
        w.writerow([row[c] for c in fieldnames] + [p.get(c, "NA") for c in extra_cols])

print(f"Joined {n_joined}/{len(rows)} vOTUs to PhaBOX2 results -> {out_path}")
PYEOF

log "Step 7 finished -> ${PHAGCN_OUT}/vOTU_taxonomy_host.tsv"
log "Next: 09 (PropagAtE), 10a/10b (CRISPR host linkage), 11 (prophage ARGs)"
