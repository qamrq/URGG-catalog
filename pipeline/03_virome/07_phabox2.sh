#!/usr/bin/env bash
###############################################################################
# 07_phabox2.sh — Phase 5: Virome pipeline, Step 7 (taxonomy/host/lifestyle)
#
# Runs PhaBOX2 on the vOTU representative catalog (8,980 seqs, from
# 06_votu_clustering.sh) to assign taxonomy (PhaGCN, ICTV 2024), host
# (CHERRY), and lifestyle -- lytic/temperate (PhaTYP).
#
# WHY THE FULL end_to_end TASK, NOT --task contamination:
#   `--task contamination` (used separately in Phase 3, on whole MAGs) is for
#   FINDING proviruses within genomic sequence -- that job is already done
#   here by VirSorter2/VIBRANT/CheckV (steps 1-3). vOTU representatives are
#   already-confirmed, already-boundary-trimmed viral sequences, so what's
#   needed now is classification of known viruses, which is PhaBOX2's default
#   end_to_end pipeline (PhaMer virus-call + PhaGCN taxonomy + PhaTYP
#   lifestyle + CHERRY host, in one run -- matches the column set seen in the
#   confirmed test run: Pred, Lineage, PhaGCNScore, TYPE, Host, CHERRYScore).
#   PhaMer's virus/non-virus call is expected to be near-100% "virus" here
#   since input is already CheckV-confirmed -- useful as a sanity check, not
#   a filtering step.
#
# TWO ITEMS THAT NEED VERIFICATION BEFORE RUNNING (not confirmed in this
# session -- pulled from an earlier chat that set up PhaBOX2 separately):
#   1. ENV_PHABOX2 below is a placeholder guess. Check the real env name:
#        conda env list | grep -i phabox
#      and correct ENV_PHABOX2 in config.sh if it differs.
#   2. GTDBTK_SUMMARY below assumes standard GTDB-Tk output filenames. Confirm:
#        ls ${GTDBTK_OUTPUT_DIR}/*.summary.tsv
#      and adjust if your actual filenames differ (e.g. classify_wf vs
#      classify subcommand can name outputs slightly differently across
#      GTDB-Tk versions).
#
# Input : VOTU_OUT/votu_representatives.fa
# Output: PHAGCN_OUT/final_prediction/final_prediction_summary.tsv
#         PHAGCN_OUT/vOTU_taxonomy_host.tsv  (joined back onto vOTU_index.tsv)
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

# --- confirmed real env name ---
ENV_PHABOX2="${ENV_PHABOX2:-phabox2}"
PHABOX2_DB="/disk252/qamar/databases/phabox2/phabox_db_v2_2"

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

###############################################################################
# Locate GTDB-Tk summary output for --bgtdb (host-taxonomy-informed CHERRY).
# Bacterial and archaeal summaries are separate files -- concatenate if both
# exist so archaeal hosts (relevant for rumen methanogens) aren't dropped.
###############################################################################
###############################################################################
# Locate GTDB-Tk summary for --bgtdb. A pre-combined bac+archaeal summary
# (gtdbtk.all.summary.tsv) already exists in this catalog's output -- use it
# directly rather than reconstructing it from the separate bac120/ar53 files.
###############################################################################
GTDBTK_SUMMARY=$(find "${GTDBTK_OUTPUT_DIR}" -maxdepth 1 -iname "gtdbtk.all.summary.tsv" 2>/dev/null | head -1)

if [[ -z "${GTDBTK_SUMMARY}" ]]; then
    # fall back to concatenating bac120 + ar53 if no pre-combined file exists
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
    BGTDB_ARG=""  # dropped -- 99%% strain catalog vs 95%% species GTDB-Tk mismatch, see chat
else
    log "WARNING: no GTDB-Tk summary found under ${GTDBTK_OUTPUT_DIR}"
    log "         Proceeding WITHOUT --bgtdb."
    BGTDB_ARG=""
fi

###############################################################################
# Run PhaBOX2 end_to_end, with --bfolder/--magonly pointed at the 99% ANI
# strain catalog so host prediction can resolve to specific genomes in your
# own catalog, not just generic external taxa.
###############################################################################
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

###############################################################################
# Join PhaBOX2 results back onto vOTU_index.tsv by representative header
###############################################################################
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
