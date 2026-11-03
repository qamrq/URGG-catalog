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
#  SCRIPT  : 04c_eggnog.sh  -  Functional annotation (eggNOG-mapper)
#  AUTHOR  : QRQ
#  INPUT   : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#  OUTPUT  : annotation/eggnog_output/nr_catalog.emapper.annotations
# ==========================================================

set +u
activate_env "${ENV_EGGNOG}"
set -u
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/04c_eggnog_$(date +%Y%m%d_%H%M%S).log") 2>&1
#activate_env "${ENV_EGGNOG}"

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/eggnog_output
mkdir -p "${OUT}"

log "Running eggNOG-mapper (DIAMOND mode)..."
emapper.py \
    -i "${QUERY}" \
    --itype proteins \
    -m diamond \
    --data_dir "${EGGNOG_DB}" \
    -o nr_catalog \
    --output_dir "${OUT}" \
    --cpu "${THREADS}" \
    --override

log "eggNOG-mapper done."
log "Annotations: ${OUT}/nr_catalog.emapper.annotations"
log "  -> KO column = KEGG orthologs; KEGG_Pathway / KEGG_Module = pathways; CAZy column present too."
