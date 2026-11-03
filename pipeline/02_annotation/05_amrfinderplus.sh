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
#  SCRIPT  : 05_amrfinderplus.sh  -  AMR genes per genome (AMRFinderPlus)
#  AUTHOR  : QRQ
#  INPUT   : species genomes (*.fa)
#  OUTPUT  : annotation/amr_results/*_amr.tsv
#            annotation/amrfinderplus_merged.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/05_amrfinder_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_AMR}"

OUT=${ANNOT}/amr_results
mkdir -p "${OUT}"

amr_one () {
    local f="$1" OUT="$2" DBP="$3"
    local base; base=$(basename "$f" ."${GENOME_EXT:-fa}")
    local res="${OUT}/${base}_amr.tsv"
    [[ -s "$res" ]] && return 0
    amrfinder -n "$f" --plus -d "${DBP}" --threads "${PARALLEL_THREADS:-2}" -o "$res" 2>/dev/null || \
        echo "WARN: amrfinder failed on ${base}" >&2
}
export -f amr_one
export GENOME_EXT PARALLEL_THREADS

log "Running AMRFinderPlus on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes, ${PARALLEL_JOBS} jobs"
log "Using database: ${AMRFINDER_DB}"
find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | \
    parallel -j "${PARALLEL_JOBS}" --joblog "${LOGS}/05_amrfinder.joblog" --bar \
    amr_one {} "${OUT}" "${AMRFINDER_DB}"

# ---- Merge into one table with a 'genome' column ----
MERGED=${ANNOT}/amrfinderplus_merged.tsv
first=$(find "${OUT}" -maxdepth 1 -name "*_amr.tsv" | head -n1)
{ printf "genome\t"; head -n1 "$first"; } > "${MERGED}"
for res in "${OUT}"/*_amr.tsv; do
    base=$(basename "$res" _amr.tsv)
    tail -n +2 "$res" | sed "s/^/${base}\t/" >> "${MERGED}"
done
log "AMR merged table -> ${MERGED}"
