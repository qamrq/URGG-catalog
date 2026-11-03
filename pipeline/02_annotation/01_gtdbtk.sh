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
#  SCRIPT  : 01_gtdbtk.sh  -  Taxonomic classification of species genomes (GTDB-Tk, R232)
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_95_species_final/dereplicated_genomes/*.fa
#  OUTPUT  : annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/01_gtdbtk_$(date +%Y%m%d_%H%M%S).log") 2>&1

export GTDBTK_DATA_PATH="${GTDBTK_DB}"

log "Using GTDB-Tk: $(gtdbtk --version | head -n1)"
log "Using GTDB reference DB: ${GTDBTK_DATA_PATH}"

OUT=${ANNOT}/gtdbtk_output_v232
mkdir -p "${OUT}"

log "Starting GTDB-Tk classify_wf on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes"

gtdbtk classify_wf \
    --genome_dir "${GENOME_DIR}" \
    --out_dir "${OUT}" \
    --extension "${GENOME_EXT}" \
    --cpus "${THREADS}" \
    --pplacer_cpus 8

log "GTDB-Tk complete."

MERGED=${OUT}/gtdbtk.all.summary.tsv
head -n1 "${OUT}/gtdbtk.bac120.summary.tsv" > "${MERGED}"
for f in "${OUT}"/gtdbtk.bac120.summary.tsv "${OUT}"/gtdbtk.ar53.summary.tsv; do
    [[ -f "$f" ]] && tail -n +2 "$f" >> "${MERGED}"
done
log "Combined taxonomy written to ${MERGED}"
