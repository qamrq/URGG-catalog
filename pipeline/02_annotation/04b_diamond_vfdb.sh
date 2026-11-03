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
#  SCRIPT  : 04b_diamond_vfdb.sh  -  Virulence factors (DIAMOND vs VFDB)
#  AUTHOR  : QRQ
#  INPUT   : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#  OUTPUT  : annotation/diamond_vfdb/diamond_vfdb_filtered.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/04b_diamond_vfdb_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DREP}"

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/diamond_vfdb
mkdir -p "${OUT}"
RAW=${OUT}/diamond_vfdb_raw.tsv
FILT=${OUT}/diamond_vfdb_filtered.tsv

log "DIAMOND blastp vs VFDB..."
diamond blastp \
    -q "${QUERY}" \
    -d "${VFDB_DMND}" \
    -o "${RAW}" \
    --threads "${THREADS}" \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --id 80 --query-cover 70 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle

awk -F'\t' '$3>=80 && $5<=1e-5' "${RAW}" > "${FILT}"

log "VFDB hits (raw): $(wc -l < "${RAW}")"
log "VFDB hits (>=80% id): $(wc -l < "${FILT}")"
log "Filtered virulence hits -> ${FILT}"
