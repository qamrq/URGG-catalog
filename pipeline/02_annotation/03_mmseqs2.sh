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
#  SCRIPT  : 03_mmseqs2.sh  -  Non-redundant gene catalog (MMseqs2, 95% identity, 90% coverage)
#  AUTHOR  : QRQ
#  INPUT   : annotation/genes/*.faa
#  OUTPUT  : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#            annotation/gene_catalog/nr_gene_catalog_cluster.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/03_mmseqs2_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DREP}"

ALL_PROT=${CATALOG}/all_proteins.faa
TMP=${CATALOG}/tmp
mkdir -p "${TMP}"

log "Concatenating all predicted proteins..."
: > "${ALL_PROT}"
for faa in "${GENES}"/*.faa; do
    base=$(basename "$faa" .faa)
    sed "s/^>/>${base}~/" "$faa" >> "${ALL_PROT}"
done
log "Total proteins: $(grep -c '^>' "${ALL_PROT}")"

log "Running MMseqs2 easy-linclust (95% id, 90% cov)..."
mmseqs easy-linclust \
    "${ALL_PROT}" \
    "${CATALOG}/nr_gene_catalog" \
    "${TMP}" \
    --min-seq-id 0.95 -c 0.9 --cov-mode 1 \
    --threads "${THREADS}"

log "Non-redundant catalog: ${CATALOG}/nr_gene_catalog_rep_seq.fasta"
log "Representatives: $(grep -c '^>' "${CATALOG}/nr_gene_catalog_rep_seq.fasta")"
log "Cluster membership: ${CATALOG}/nr_gene_catalog_cluster.tsv"
rm -rf "${TMP}"
