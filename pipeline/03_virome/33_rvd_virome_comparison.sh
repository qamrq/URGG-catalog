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
#  SCRIPT  : 33_rvd_virome_comparison.sh  -  vOTU matches to the Rumen Virome Database (95% ANI, 85% AF)
#  AUTHOR  : QRQ
#  INPUT   : votu_representatives.fa
#            $URGG_DB/rvd/RVD.fa
#  OUTPUT  : phase_5_virome/09_rvd_comparison/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

PROJ_DIR="${URGG_PROJECT}/phase_5_virome"
WORK_DIR="${PROJ_DIR}/09_rvd_comparison"
QUERY_FA="${PROJ_DIR}/07_votu_clustering/votu_representatives.fa"
RVD_FA="${URGG_DB}/rvd/RVD.fa"
RVD_DB="${WORK_DIR}/rvd_blastdb/rvd_db"

THREADS=64

mkdir -p "${WORK_DIR}/rvd_blastdb"
cd "${WORK_DIR}"

if [ ! -f "${RVD_DB}.nin" ] && [ ! -f "${RVD_DB}.00.nin" ]; then
    echo "[1/4] Building BLAST database from RVD.fa ..."
    makeblastdb -in "${RVD_FA}" -dbtype nucl -out "${RVD_DB}"
else
    echo "[1/4] RVD BLAST database already exists, skipping build."
fi

BLAST_OUT="${WORK_DIR}/votu_vs_rvd_blast.tsv"

if [ ! -s "${BLAST_OUT}" ]; then
    echo "[2/4] Running BLASTN (query vOTUs vs RVD) ..."
    blastn -query "${QUERY_FA}" \
           -db "${RVD_DB}" \
           -outfmt '6 std qlen slen' \
           -max_target_seqs 25 \
           -perc_identity 90 \
           -num_threads "${THREADS}" \
           -out "${BLAST_OUT}"
else
    echo "[2/4] BLAST output already exists, skipping."
fi

ANI_OUT="${WORK_DIR}/votu_vs_rvd_ani.tsv"

echo "[3/4] Running anicalc ..."
anicalc -i "${BLAST_OUT}" -o "${ANI_OUT}"

echo "[4/4] Classifying matched vs novel vOTUs ..."
python3 "$(dirname "$0")/classify_rvd_matches.py" \
    --query-fa "${QUERY_FA}" \
    --ani "${ANI_OUT}" \
    --min-ani 95 \
    --min-cov 85 \
    --out-table "${WORK_DIR}/votu_rvd_classification.tsv" \
    --out-summary "${WORK_DIR}/votu_rvd_summary.txt"

echo "Done. See:"
echo "  ${WORK_DIR}/votu_rvd_classification.tsv  (per-vOTU matched/novel calls)"
echo "  ${WORK_DIR}/votu_rvd_summary.txt          (matched/novel % summary)"

