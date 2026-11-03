#!/usr/bin/env bash
# ============================================================
# 01_gtdbtk.sh  —  Taxonomy classification of the 95% ANI species catalog
#                  (25,473 representative genomes)
# Uses EXISTING GTDB-Tk v2.5.2 install (disk192) + R226 DB (disk212).
# No install/download needed for this step.
# RAM: heavy — budget ~320 GB. If constrained, use --batchfile splitting.
# This taxonomy run is also what's needed before any novelty claims vs.
# the comparator catalogs (RUG2, RGMGC, etc.) noted in your dRep export doc.
# ============================================================
set -euo pipefail
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/01_gtdbtk_$(date +%Y%m%d_%H%M%S).log") 2>&1

# GTDB-Tk lives in a DIFFERENT miniconda base (disk192) than the rest of the
# envs (which use the default `conda activate` via activate_env() in
# config.sh). Activate it directly rather than through activate_env().
#source "${GTDBTK_CONDA_BASE}/bin/activate" "${GTDBTK_ENV_PATH}"
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