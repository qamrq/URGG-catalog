#!/usr/bin/env bash
# ============================================================
# 06_platon.sh  —  Plasmid contig detection (Platon)
# Env: platon | Input: GENOME_DIR/*.fa (nucleotide) | Per-genome, parallel
# Strain-level plasmid content is meaningful (plasmids differ between
# near-identical strains), so run per genome, not on the collapsed catalog.
# Resumable.
# ============================================================
set -euo pipefail
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/06_platon_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_PLATON}"

OUT=${ANNOT}/platon_output
mkdir -p "${OUT}"

platon_one () {
    local f="$1" OUT="$2" DBP="$3"
    local base; base=$(basename "$f" ."${GENOME_EXT:-fa}")
    local odir="${OUT}/${base}"
    [[ -s "${odir}/${base}.tsv" ]] && return 0
    mkdir -p "${odir}"
    platon --db "${DBP}" --output "${odir}" --threads "${PARALLEL_THREADS:-2}" --prefix "${base}" "$f" \
        2>/dev/null || echo "WARN: platon failed on ${base}" >&2
}
export -f platon_one
export GENOME_EXT PARALLEL_THREADS

log "Running Platon on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes, ${PARALLEL_JOBS} jobs"
find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | \
    parallel -j "${PARALLEL_JOBS}" --joblog "${LOGS}/06_platon.joblog" --bar \
    platon_one {} "${OUT}" "${PLATON_DB}"

# ---- Merge plasmid summaries ----
MERGED=${ANNOT}/platon_merged.tsv
: > "${MERGED}"
header_written=0
for tsv in "${OUT}"/*/*.tsv; do
    base=$(basename "$(dirname "$tsv")")
    if [[ $header_written -eq 0 ]]; then
        { printf "genome\t"; head -n1 "$tsv"; } > "${MERGED}"
        header_written=1
    fi
    tail -n +2 "$tsv" | sed "s/^/${base}\t/" >> "${MERGED}"
done
log "Platon merged plasmid table -> ${MERGED}"
