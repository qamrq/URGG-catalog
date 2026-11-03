#!/bin/bash
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
#  SCRIPT  : 34_rgvc_virome_comparison.sh  -  vOTU matches to RGVC (95% ANI, 85% AF)
#  AUTHOR  : QRQ
#  INPUT   : votu_representatives.fa
#            $URGG_DB/rgvc/uniq.viral.fa
#  OUTPUT  : phase_5_virome/10_rgvc_comparison/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

PROJ_DIR="${URGG_PROJECT}/phase_5_virome"
WORK_DIR="${PROJ_DIR}/10_rgvc_comparison"
QUERY_FA="${PROJ_DIR}/07_votu_clustering/votu_representatives.fa"
RGVC_FA="${URGG_DB}/rgvc/uniq.viral.fa"
RGVC_DB="${WORK_DIR}/rgvc_blastdb/rgvc_db"

THREADS=64

mkdir -p "${WORK_DIR}/rgvc_blastdb"
cd "${WORK_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "[0/4] Verifying RGVC.fa sequence count against the published paper (expect 43,981)..."
RGVC_COUNT=$(grep -c "^>" "${RGVC_FA}")
log "  RGVC.fa contains ${RGVC_COUNT} sequences."
if [[ "${RGVC_COUNT}" -ne 43981 ]]; then
    log "  WARNING: count does not match the paper's reported 43,981 -- confirm this is the correct/complete file before proceeding (could be a different RGVC version, a partial download, or duplicate/malformed headers)."
fi

if [ ! -f "${RGVC_DB}.nin" ] && [ ! -f "${RGVC_DB}.00.nin" ]; then
    log "[1/4] Building BLAST database from RGVC uniq.viral.fa ..."
    makeblastdb -in "${RGVC_FA}" -dbtype nucl -out "${RGVC_DB}"
else
    log "[1/4] RGVC BLAST database already exists, skipping build."
fi

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

ANI_OUT="${WORK_DIR}/votu_vs_rgvc_ani.tsv"

log "[3/4] Running anicalc ..."
anicalc -i "${BLAST_OUT}" -o "${ANI_OUT}"

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
