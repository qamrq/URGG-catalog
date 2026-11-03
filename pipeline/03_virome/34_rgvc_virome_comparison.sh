#!/bin/bash
# 34_rgvc_virome_comparison.sh
# Compare this project's 8,980 species-level vOTUs against RGVC
# (Ruminant Gastrointestinal Virome Catalog; 43,981 vOTUs), a more
# comprehensive ruminant-specific virome resource than RVD.
#
# METHODOLOGY: IDENTICAL to the RVD comparison (33_rvd_virome_comparison.sh)
# -- RGVC's own construction used the same 95% ANI / 85% aligned-fraction
# MIUViG species-demarcation threshold as RVD, confirmed directly from the
# RGVC paper's Methods. No new threshold research needed; this script
# reuses anicalc/aniclust exactly as before, pointed at RGVC's reference
# FASTA instead of RVD's.
#
# Tiered classification (species/genus/distant/no-hit) also carries over
# unchanged from classify_rvd_matches.py -- same script, no modification,
# just different input/output paths.

set -euo pipefail

#conda activate checkv  # provides anicalc/aniclust, confirmed on PATH earlier in this project

# ---- paths -------------------------------------------------------------
PROJ_DIR="/disk252/qamar/metagenome_analysis/phase_5_virome"
WORK_DIR="${PROJ_DIR}/10_rgvc_comparison"
QUERY_FA="${PROJ_DIR}/07_votu_clustering/votu_representatives.fa"
RGVC_FA="/disk252/qamar/databases/rgvc/uniq.viral.fa"
RGVC_DB="${WORK_DIR}/rgvc_blastdb/rgvc_db"

THREADS=64

mkdir -p "${WORK_DIR}/rgvc_blastdb"
cd "${WORK_DIR}"

# ---- 0. VERIFY the reference matches the paper's reported count --------
# (43,981 vOTUs per RGVC paper -- confirm before trusting anything
# downstream, same discipline as the RVD.fa count-check earlier in this
# project)
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "[0/4] Verifying RGVC.fa sequence count against the published paper (expect 43,981)..."
RGVC_COUNT=$(grep -c "^>" "${RGVC_FA}")
log "  RGVC.fa contains ${RGVC_COUNT} sequences."
if [[ "${RGVC_COUNT}" -ne 43981 ]]; then
    log "  WARNING: count does not match the paper's reported 43,981 -- confirm this is the correct/complete file before proceeding (could be a different RGVC version, a partial download, or duplicate/malformed headers)."
fi

# ---- 1. Build BLAST db on RGVC (build once, reuse across runs) ---------
if [ ! -f "${RGVC_DB}.nin" ] && [ ! -f "${RGVC_DB}.00.nin" ]; then
    log "[1/4] Building BLAST database from RGVC uniq.viral.fa ..."
    makeblastdb -in "${RGVC_FA}" -dbtype nucl -out "${RGVC_DB}"
else
    log "[1/4] RGVC BLAST database already exists, skipping build."
fi

# ---- 2. All-vs-RGVC BLASTN (query is small: 8,980 seqs) -----------------
BLAST_OUT="${WORK_DIR}/votu_vs_rgvc_blast.tsv"

if [ ! -s "${BLAST_OUT}" ]; then
    log "[2/4] Running BLASTN (query vOTUs vs RGVC) ..."
    blastn -query "${QUERY_FA}" \
           -db "${RGVC_DB}" \
           -outfmt '6 std qlen slen' \
           -max_target_seqs 25 \
           -perc_identity 90 \
           -num_threads "${THREADS}" \
           -out "${BLAST_OUT}"
else
    log "[2/4] BLAST output already exists, skipping."
fi

# ---- 3. anicalc: pairwise ANI + alignment fraction ----------------------
ANI_OUT="${WORK_DIR}/votu_vs_rgvc_ani.tsv"

log "[3/4] Running anicalc ..."
anicalc -i "${BLAST_OUT}" -o "${ANI_OUT}"

# ---- 4. Tiered classification -------------------------------------------
# REUSES classify_rvd_matches.py UNCHANGED -- same species/genus/distant/
# no-hit tiers (95/85 species, 80/50 genus), same script, different
# input/output paths only.
log "[4/4] Classifying vOTUs into species/genus/distant/no-hit tiers vs RGVC..."
python3 "$(dirname "$0")/classify_rvd_matches.py" \
    --query-fa "${QUERY_FA}" \
    --ani "${ANI_OUT}" \
    --species-ani 95 --species-cov 85 \
    --genus-ani 80 --genus-cov 50 \
    --out-table "${WORK_DIR}/votu_rgvc_classification.tsv" \
    --out-plot-csv "${WORK_DIR}/votu_rgvc_for_plotting.csv" \
    --out-summary "${WORK_DIR}/votu_rgvc_summary.txt"

log "Done. See:"
log "  ${WORK_DIR}/votu_rgvc_classification.tsv  (per-vOTU tiered calls)"
log "  ${WORK_DIR}/votu_rgvc_for_plotting.csv     (slim CSV, ready for figure)"
log "  ${WORK_DIR}/votu_rgvc_summary.txt          (tier % summary)"
