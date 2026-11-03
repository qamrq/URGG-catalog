#!/usr/bin/env bash
# ============================================================
# 02_prodigal.sh  —  Gene prediction per genome (Prodigal)
# Env: drep | Input: GENOME_DIR/*.fa | Output: GENES/{base}.faa + .gff
#
# NOTE on MetaGeneMark: the proposal lists MetaGeneMark, but it needs a
# license key and is awkward to batch. Prodigal is the standard open-source
# equivalent and is what most MAG catalogs use. We run `-p single` because
# dRep reps are near-complete single genomes; add a `-p meta` fallback for
# any tiny/fragmented input that single-mode rejects.
# Resumable: skips genomes already done. Parallelised with GNU parallel.
# ============================================================
set -euo pipefail
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/02_prodigal_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DREP}"

predict_one () {
    local f="$1" GENES="$2"
    local base; base=$(basename "$f" ."${GENOME_EXT:-fa}")
    local faa="${GENES}/${base}.faa"
    local gff="${GENES}/${base}.gff"
    [[ -s "$faa" ]] && return 0   # already done
    # try single mode first, fall back to meta for short/fragmented genomes
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
