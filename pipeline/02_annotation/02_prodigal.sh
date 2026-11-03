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
#  SCRIPT  : 02_prodigal.sh  -  Gene prediction per genome (Prodigal)
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_95_species_final/dereplicated_genomes/*.fa
#  OUTPUT  : annotation/genes/<genome>.faa, .gff
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/02_prodigal_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DREP}"

predict_one () {
    local f="$1" GENES="$2"
    local base; base=$(basename "$f" ."${GENOME_EXT:-fa}")
    local faa="${GENES}/${base}.faa"
    local gff="${GENES}/${base}.gff"
    [[ -s "$faa" ]] && return 0
    if ! prodigal -i "$f" -a "$faa" -o "$gff" -f gff -p single -q 2>/dev/null; then
        prodigal -i "$f" -a "$faa" -o "$gff" -f gff -p meta -q
    fi
}
export -f predict_one
export GENOME_EXT

log "Predicting genes on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes with ${PARALLEL_JOBS} jobs"

find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | \
    parallel -j "${PARALLEL_JOBS}" --joblog "${LOGS}/02_prodigal.joblog" --bar \
    predict_one {} "${GENES}"

log "Prodigal complete. $(ls "${GENES}"/*.faa | wc -l) protein files in ${GENES}"
