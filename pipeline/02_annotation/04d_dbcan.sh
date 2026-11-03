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
#  SCRIPT  : 04d_dbcan.sh  -  CAZyme annotation (run_dbcan)
#  AUTHOR  : QRQ
#  INPUT   : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#  OUTPUT  : annotation/dbcan_output/overview.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/04d_dbcan_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DBCAN}"

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/dbcan_output
mkdir -p "${OUT}"

log "Running run_dbcan CAZyme_annotation (dbCAN v5.x subcommand syntax)..."
run_dbcan CAZyme_annotation \
    --input_raw_data "${QUERY}" \
    --mode protein \
    --output_dir "${OUT}" \
    --db_dir "${DBCAN_DB}" \
    --threads "${THREADS}"

log "dbCAN done. Consensus CAZyme table: ${OUT}/overview.txt"
