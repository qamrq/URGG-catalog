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
#  SCRIPT  : 10p_votus.sh  -  vOTU clustering, PhaBOX2 taxonomy/host/lifestyle, membership and host-sharing tiers
#  AUTHOR  : QRQ
#  INPUT   : confirmed_viral.fa
#            GTDB-Tk summaries
#            $URGG_DB/phabox2/phabox_db_v2_2
#  OUTPUT  : votu_representatives.fa, vOTU_index.tsv, vOTU_taxonomy_host.tsv
#            phase_5_virome/08_results/votu_membership_full.tsv, cross_host_votu_matrix.tsv, votu_host_sharing_tiers.tsv
#  USAGE   : bash 10p_votus.sh [cluster phabox2 membership tables tiers]   (no argument: all steps in this order)
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

step_cluster() (
source "$(dirname "$0")/config_virome.sh"

log "Step 6: vOTU clustering (95% ANI / 85% AF) on confirmed viral set"

CONFIRMED_FA="${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa"
if [[ ! -s "${CONFIRMED_FA}" ]]; then
    log "ERROR: ${CONFIRMED_FA} not found or empty. Run 8p_virus_detection.sh checkv first."
    exit 1
fi
n_input=$(grep -c '^>' "${CONFIRMED_FA}")
log "Input: ${n_input} confirmed viral sequences"

ANICALC="${ANICALC:-$(command -v anicalc.py || true)}"
ANICLUST="${ANICLUST:-$(command -v aniclust.py || true)}"

if [[ ! -f "${ANICALC}" || ! -f "${ANICLUST}" ]]; then
    log "ERROR: could not locate or fetch anicalc.py/aniclust.py."
    log "       Set ANICALC and ANICLUST to their paths."
    log "       Source: https://github.com/snayfach/CheckV/tree/master/scripts"
    exit 1
fi
log "Using anicalc: ${ANICALC}"
log "Using aniclust: ${ANICLUST}"

BLAST_DB="${VOTU_OUT}/confirmed_viral_blastdb"
BLAST_OUT="${VOTU_OUT}/blast_all_vs_all.tsv"

if [[ ! -s "${BLAST_OUT}" ]]; then
    log "Building BLAST db..."
    makeblastdb -in "${CONFIRMED_FA}" -dbtype nucl -out "${BLAST_DB}" \
        > "${LOG_DIR}/votu_makeblastdb.log" 2>&1

    log "Running all-vs-all megablast (this is O(n^2)-ish; ${n_input} seqs -- watch runtime)..."
    blastn -query "${CONFIRMED_FA}" -db "${BLAST_DB}" \
        -outfmt '6 std qlen slen' \
        -max_target_seqs 10000 \
        -perc_identity 90 \
        -num_threads "${THREADS}" \
        -out "${BLAST_OUT}" \
        > "${LOG_DIR}/votu_blastn.log" 2>&1
else
    log "Existing BLAST output found -- skipping (delete ${BLAST_OUT} to force a rerun)"
fi

n_pairs=$(wc -l < "${BLAST_OUT}")
log "  ${n_pairs} pairwise hits"

log "Computing pairwise ANI/AF (anicalc.py)..."
ANI_OUT="${VOTU_OUT}/ani.tsv"
python "${ANICALC}" -i "${BLAST_OUT}" -o "${ANI_OUT}" \
    > "${LOG_DIR}/votu_anicalc.log" 2>&1

if [[ ! -s "${ANI_OUT}" ]]; then
    log "ERROR: anicalc.py produced no output -- see ${LOG_DIR}/votu_anicalc.log"
    exit 1
fi

log "Clustering into vOTUs (95% ANI, 85% AF)..."
CLUSTERS_OUT="${VOTU_OUT}/votu_clusters.tsv"
python "${ANICLUST}" \
    --fna "${CONFIRMED_FA}" \
    --ani "${ANI_OUT}" \
    --out "${CLUSTERS_OUT}" \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0 \
    > "${LOG_DIR}/votu_aniclust.log" 2>&1

if [[ ! -s "${CLUSTERS_OUT}" ]]; then
    log "ERROR: aniclust.py produced no output -- see ${LOG_DIR}/votu_aniclust.log"
    exit 1
fi

n_votus=$(wc -l < "${CLUSTERS_OUT}")
log "  ${n_votus} vOTUs from ${n_input} input sequences"

log "Extracting vOTU representative sequences..."
cut -f1 "${CLUSTERS_OUT}" | sort -u > "${VOTU_OUT}/.rep_ids.txt"
sed 's/^/^/' "${VOTU_OUT}/.rep_ids.txt" > "${VOTU_OUT}/.rep_ids_regex.txt"

if ! command -v seqkit >/dev/null 2>&1; then
    log "ERROR: seqkit required on PATH"
    exit 1
fi
seqkit grep -n -r -f "${VOTU_OUT}/.rep_ids_regex.txt" "${CONFIRMED_FA}" \
    > "${VOTU_OUT}/votu_representatives.fa"

n_reps=$(grep -c '^>' "${VOTU_OUT}/votu_representatives.fa" 2>/dev/null || echo 0)
log "  ${n_reps} representative sequences extracted"
if [[ "${n_reps}" -lt "${n_votus}" ]]; then
    log "  WARNING: extracted ${n_reps} but expected ${n_votus} -- check for CheckV suffix"
    log "  mismatches (same class of issue documented in 8p_virus_detection.sh checkv)"
fi

log "Building vOTU index with member MAG provenance..."
awk -F'\t' '
    BEGIN { OFS="\t"; print "votu_id","representative","n_members","member_ids" }
    {
        rep = $1
        n = split($2, members, ",")
        votu_id = sprintf("vOTU_%06d", NR)
        print votu_id, rep, n, $2
    }
' "${CLUSTERS_OUT}" > "${VOTU_OUT}/vOTU_index.tsv"

rm -f "${VOTU_OUT}/.rep_ids.txt" "${VOTU_OUT}/.rep_ids_regex.txt"

log "Step 6 finished."
log "  ${n_votus} vOTUs -> ${VOTU_OUT}/votu_representatives.fa"
log "  cluster membership -> ${VOTU_OUT}/vOTU_index.tsv"
log "  source MAG / host species: join member_ids (M###### prefix) against ${MAG_INDEX_TSV}"
)

step_phabox2() (
source "$(dirname "$0")/config_virome.sh"

# --- confirmed real env name ---
PHABOX2_DB="${URGG_DB}/phabox2/phabox_db_v2_2"

log "Step 7: PhaBOX2 taxonomy/host/lifestyle on vOTU representatives"

VOTU_FA="${VOTU_OUT}/votu_representatives.fa"
if [[ ! -s "${VOTU_FA}" ]]; then
    log "ERROR: ${VOTU_FA} not found. Run 10p_votus.sh cluster first."
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

log "PhaBOX2 finished -> ${PHAGCN_OUT}/vOTU_taxonomy_host.tsv"
)

step_membership() (
source "$(dirname "$0")/config_virome.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
VOTU_INDEX="${VOTU_OUT}/vOTU_index.tsv"

if [[ ! -s "${VOTU_INDEX}" ]]; then
    log "ERROR: ${VOTU_INDEX} not found."
    exit 1
fi

log "Exploding vOTU_index.tsv into one row per sequence..."

python3 << PYEOF
import csv

out_path = "${RESULTS_DIR}/votu_membership_full.tsv"
n_rows = 0
n_clusters = 0

with open("${VOTU_INDEX}") as fin, open(out_path, "w", newline="") as fout:
    r = csv.DictReader(fin, delimiter='\t')
    w = csv.writer(fout, delimiter='\t')
    w.writerow(["sequence_id", "votu_id", "representative", "is_representative", "n_members"])

    for row in r:
        n_clusters += 1
        votu_id = row["votu_id"]
        rep = row["representative"]
        n_members = row["n_members"]
        members = [m.strip() for m in row["member_ids"].split(",") if m.strip()]

        for member in members:
            is_rep = (member == rep)
            w.writerow([member, votu_id, rep, is_rep, n_members])
            n_rows += 1

print(f"Exploded {n_clusters} clusters -> {n_rows} sequence-level rows")
print(f"-> {out_path}")
PYEOF

n_check=$(($(wc -l < "${RESULTS_DIR}/votu_membership_full.tsv") - 1))
log "Check: ${n_check} sequence-level rows written (expected about 15,301)"
log "Step 07b finished."
)

step_tables() (
source "$(dirname "$0")/config_virome.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
mkdir -p "${RESULTS_DIR}"
log "Writing results tables -> ${RESULTS_DIR}"

log "[1/6] Discovery summary table..."
n_mags=$(find "${MAGS_99ANI_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l)
n_vs2=$(cat "${VS2_OUT}"/*/final-viral-score.tsv 2>/dev/null | grep -v "^seqname" | wc -l)
n_vib=$(cat "${VIBRANT_OUT}"/*/VIBRANT_batch_*/VIBRANT_phages_batch_*/*.phages_combined.fna 2>/dev/null | grep -c '^>')
n_union=$(grep -c '^>' "${CONFIRMED_VIRAL_OUT}/all_viral_union.fa" 2>/dev/null || echo NA)
n_confirmed=$(grep -c '^>' "${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa" 2>/dev/null || echo NA)
n_votus=$(($(wc -l < "${VOTU_OUT}/vOTU_index.tsv") - 1))

{
  printf 'stage\tcount\n'
  printf 'MAGs_screened\t%s\n' "${n_mags}"
  printf 'VirSorter2_calls\t%s\n' "${n_vs2}"
  printf 'VIBRANT_calls\t%s\n' "${n_vib}"
  printf 'Union_pre_filter\t%s\n' "${n_union}"
  printf 'CheckV_confirmed\t%s\n' "${n_confirmed}"
  printf 'Nonredundant_vOTUs\t%s\n' "${n_votus}"
} > "${RESULTS_DIR}/table1_discovery_summary.tsv"
cat "${RESULTS_DIR}/table1_discovery_summary.tsv"

log "[2/6] CheckV quality tier distribution..."
awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="checkv_quality") q=i; next }
    { print $q }
' "${CHECKV_OUT}/quality_summary.tsv" | sort | uniq -c | sort -rn \
  | awk '{print $2"\t"$1}' > "${RESULTS_DIR}/checkv_quality_tiers.tsv"
sed -i '1i quality_tier\tcount' "${RESULTS_DIR}/checkv_quality_tiers.tsv"
cat "${RESULTS_DIR}/checkv_quality_tiers.tsv"

log "[3/6] vOTU size distribution..."
awk -F'\t' 'NR>1 {print $3}' "${VOTU_OUT}/vOTU_index.tsv" \
  | sort -n | uniq -c | awk '{print $2"\t"$1}' \
  > "${RESULTS_DIR}/votu_size_distribution.tsv"
sed -i '1i cluster_size\tn_votus' "${RESULTS_DIR}/votu_size_distribution.tsv"
n_singletons=$(awk -F'\t' 'NR>1 && $3==1' "${VOTU_OUT}/vOTU_index.tsv" | wc -l)
log "  ${n_singletons} / $(($(wc -l < "${VOTU_OUT}/vOTU_index.tsv")-1)) vOTUs are singletons"

log "[4/6] Cross-host vOTU distribution..."
python3 << PYEOF
import csv, sys
from collections import defaultdict

RAW_LABEL_TO_SPECIES = {
    "Cattle": "Bos_taurus", "Cow": "Bos_taurus", "Beef_cattle": "Bos_taurus",
    "Beef_Heifers": "Bos_taurus", "Dairy_cattle": "Bos_taurus",
    "Nellore_cattle": "Bos_taurus", "Boran_cattle": "Bos_taurus",
    "Yak": "Bos_grunniens", "Gaur": "Bos_frontalis", "Buffalo": "Bubalus_bubalis",
    "Sheep": "Ovis_aries", "Goat": "Capra_hircus",
    "Antelope": "Pantholops_hodgsonii", "Gazelle": "Procapra_picticaudata",
    "Camel": "Camelus_dromedarius", "Giraffe": "Giraffa_camelopardalis",
    "Moose": "Alces_alces", "North_American_Moose": "Alces_americanus",
    "Red_deer": "Cervus_elaphus", "Reindeer": "Rangifer_tarandus",
    "Water_deer": "Hydropotes_inermis",
}
def normalize_host(raw_label):
    return RAW_LABEL_TO_SPECIES.get(raw_label, f"UNMAPPED_{raw_label}")

mag_to_host = {}
unmapped_seen = set()
with open("${MAG_INDEX_TSV}") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        mag_id = row["mag_id"]
        parts = mag_id.split("__")
        raw_label = parts[1] if len(parts) > 1 else "UNKNOWN"
        host = normalize_host(raw_label)
        if host.startswith("UNMAPPED_"):
            unmapped_seen.add(raw_label)
        mag_to_host[row["mag_index"]] = host

if unmapped_seen:
    print(f"WARNING: {len(unmapped_seen)} raw label(s) have NO species mapping "
          f"and were left as UNMAPPED_<label> (treated as their own group -- "
          f"add these to RAW_LABEL_TO_SPECIES for a correct result): {sorted(unmapped_seen)}")

all_hosts = sorted(set(mag_to_host.values()))
print(f"Detected {len(all_hosts)} CONSOLIDATED host species: {all_hosts}")

votu_hosts = defaultdict(set)
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        votu_id = row["votu_id"]
        members = row["member_ids"].split(",")
        for m in members:
            mag_idx = m.strip()[:7]
            host = mag_to_host.get(mag_idx, "UNKNOWN")
            votu_hosts[votu_id].add(host)

matrix_path = "${RESULTS_DIR}/cross_host_votu_matrix.tsv"
with open(matrix_path, "w", newline="") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(["votu_id"] + all_hosts + ["n_hosts", "category"])
    n_core = n_shared = n_specific = 0
    for votu_id, hosts in votu_hosts.items():
        presence = [1 if h in hosts else 0 for h in all_hosts]
        n_h = len(hosts)
        if n_h == len(all_hosts):
            cat = "core"; n_core += 1
        elif n_h > 1:
            cat = "shared"; n_shared += 1
        else:
            cat = "host_specific"; n_specific += 1
        w.writerow([votu_id] + presence + [n_h, cat])

summary_path = "${RESULTS_DIR}/cross_host_summary.tsv"
with open(summary_path, "w") as f:
    f.write("category\tn_votus\n")
    f.write(f"core\t{n_core}\n")
    f.write(f"shared\t{n_shared}\n")
    f.write(f"host_specific\t{n_specific}\n")

print(f"core={n_core} shared={n_shared} host_specific={n_specific}")
print(f"-> {matrix_path}")
print(f"-> {summary_path}")
PYEOF

log "[5/6] Length/completeness table..."
awk -F'\t' '
    NR==1 {
        for(i=1;i<=NF;i++){
            if($i=="contig_id") id=i
            if($i=="contig_length") len=i
            if($i=="completeness") comp=i
            if($i=="checkv_quality") q=i
        }
        print "contig_id\tlength\tcompleteness\tquality"; next
    }
    { print $id"\t"$len"\t"$comp"\t"$q }
' "${CHECKV_OUT}/quality_summary.tsv" > "${RESULTS_DIR}/length_completeness_all.tsv"

awk -F'\t' 'NR==FNR{keep[$1]=1;next} FNR==1{print;next} ($1 in keep)' \
    "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt" "${RESULTS_DIR}/length_completeness_all.tsv" \
    > "${RESULTS_DIR}/length_completeness.tsv"
log "  $(($(wc -l < "${RESULTS_DIR}/length_completeness.tsv")-1)) confirmed entries with length/completeness"

log "[6/6] AMG summary across all batches..."
AMG_COMBINED="${RESULTS_DIR}/amg_summary_combined.tsv"
first=1
for f in "${DRAMV_OUT}"/batch_*/distill/amg_summary.tsv; do
    [[ -f "${f}" ]] || continue
    if [[ "${first}" -eq 1 ]]; then
        cp "${f}" "${AMG_COMBINED}"; first=0
    else
        tail -n +2 "${f}" >> "${AMG_COMBINED}"
    fi
done
[[ -f "${DRAMV_OUT}/vibrant_only/distill/amg_summary.tsv" ]] && \
    tail -n +2 "${DRAMV_OUT}/vibrant_only/distill/amg_summary.tsv" >> "${AMG_COMBINED}" 2>/dev/null

n_amgs=$(($(wc -l < "${AMG_COMBINED}" 2>/dev/null || echo 1) - 1))
log "  ${n_amgs} total AMG entries -> ${AMG_COMBINED}"
if [[ "${n_amgs}" -le 0 ]]; then
    log "  EXPECTED (see script comments): AMG scoring lost flanking context during"
    log "  boundary trimming. Treat as a stated methods limitation, not a bug."
fi

log "All 6 result tables written to ${RESULTS_DIR}"
)

step_tiers() (
source "$(dirname "$0")/config_virome.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
MATRIX="${RESULTS_DIR}/cross_host_votu_matrix.tsv"

if [[ ! -s "${MATRIX}" ]]; then
    log "ERROR: ${MATRIX} not found -- run 10p_votus.sh tables first."
    exit 1
fi

log "Building standalone vOTU host-sharing tier file ..."

python3 << PYEOF
import csv

votu_id_to_rep = {}
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        votu_id_to_rep[row["votu_id"]] = row["representative"]

out_path = "${RESULTS_DIR}/votu_host_sharing_tiers.tsv"
n_specific = n_shared = n_core = 0

with open("${MATRIX}") as fin, open(out_path, "w", newline="") as fout:
    r = csv.DictReader(fin, delimiter='\t')
    host_cols = [c for c in r.fieldnames if c not in ("votu_id", "n_hosts", "category")]
    w = csv.writer(fout, delimiter='\t')
    w.writerow(["votu_id", "representative", "n_hosts", "host_sharing_tier", "host_species_list"])

    for row in r:
        n_hosts = int(row["n_hosts"])
        if n_hosts == 1:
            tier = "host_specific"; n_specific += 1
        elif n_hosts <= 3:
            tier = "shared"; n_shared += 1
        else:
            tier = "core_broad"; n_core += 1

        present_species = [h for h in host_cols if row[h] == "1"]
        rep = votu_id_to_rep.get(row["votu_id"], "NA")

        w.writerow([row["votu_id"], rep, n_hosts, tier, ";".join(present_species)])

total = n_specific + n_shared + n_core
print(f"=== vOTU host-sharing tiers ===")
print(f"Host-specific (1 species) : {n_specific} ({100*n_specific/total:.1f}%)")
print(f"Shared (2-3 species)      : {n_shared} ({100*n_shared/total:.1f}%)")
print(f"Core/broad (>3 species)   : {n_core} ({100*n_core/total:.1f}%)")
print(f"Total                     : {total}")
print(f"-> {out_path}")
PYEOF

log "Step 08f finished."
)

for s in ${@:-cluster phabox2 membership tables tiers}; do
  echo "==> 10p_votus.sh: $s"
  "step_$s"
done
