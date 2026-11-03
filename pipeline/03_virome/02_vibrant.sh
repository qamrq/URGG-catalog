#!/usr/bin/env bash
###############################################################################
# 02_vibrant.sh — Phase 5: Virome pipeline, Step 1 (second caller, for union)
#
# Runs VIBRANT on the SAME batched input as VirSorter2 so the two callers can
# be unioned before CheckV filtering (union > intersection; CheckV quality
# tiers do the real filtering downstream).
#
# Input : BATCH_DIR/batch_*.fa   (headers already tagged M######__)
# Output: VIBRANT_OUT/<batch_id>/VIBRANT_<batch_id>/...
#
# -------------------------------------------------------------------------
# NOTE ON THE -virome FLAG  (changed from the original plan -- please read)
# -------------------------------------------------------------------------
# The original pipeline sketch specified VIBRANT with `-virome`, taken from
# Yan & Yu 2024. But in that paper VIBRANT was applied to sequences VirSorter2
# had ALREADY called as viral -- i.e. the input genuinely was a virome, and
# the flag was appropriate.
#
# Here VIBRANT runs independently on whole MAGs, which are bacterial/archaeal
# genomes with prophages embedded in them: a mixed, overwhelmingly cellular
# input, not a virome. VIBRANT's own documentation scopes `-virome` to virome
# datasets specifically, and it changes how aggressively non-viral scaffolds
# are handled. Applying it to genomic input does not match what it is for.
#
# This script therefore runs WITHOUT `-virome`. If you would rather mirror the
# paper exactly, the coherent way to do it is to change the design -- run
# VIBRANT on VirSorter2 output as a confirmation step (intersection) rather
# than as an independent caller (union) -- rather than to keep the flag on
# genomic input. Worth a decision before the full run, since it affects how
# the two callers combine and how you describe it in the methods.
# -------------------------------------------------------------------------
#
# Install quirk: bioconda names the executable VIBRANT_run.py, not VIBRANT.py.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"
activate_env "${ENV_VIBRANT}"

log "Step 1 (cont.): VIBRANT on batched MAGs (${BATCH_DIR})"

VIB_SHARE=$(find "$(conda info --base)/envs/${ENV_VIBRANT}/share" -maxdepth 1 -iname "vibrant-*" 2>/dev/null | head -1)
VIBRANT_DB_DIR="${VIB_SHARE}/db/databases"
VIBRANT_FILES_DIR="${VIB_SHARE}/db/files"

if [[ ! -d "${VIBRANT_DB_DIR}" ]]; then
    log "ERROR: VIBRANT database dir not found at ${VIBRANT_DB_DIR}"
    exit 1
fi
log "VIBRANT databases: ${VIBRANT_DB_DIR}"

if ! compgen -G "${BATCH_DIR}/batch_*.fa" > /dev/null; then
    log "ERROR: no batch FASTAs in ${BATCH_DIR}. Run 00c_build_batches.sh first."
    exit 1
fi

run_one_batch () {
    local batch_fa="$1"
    local batch_id
    batch_id="$(basename "${batch_fa}" .fa)"
    local out_dir="${VIBRANT_OUT}/${batch_id}"

    if [[ -d "${out_dir}/VIBRANT_${batch_id}" ]]; then
        echo "[skip] ${batch_id}"
        return 0
    fi

    mkdir -p "${out_dir}"

    VIBRANT_run.py \
        -i "${batch_fa}" \
        -folder "${out_dir}" \
        -t "${BATCH_THREADS}" \
        -no_plot \
        -d "${VIBRANT_DB_DIR}" \
        -m "${VIBRANT_FILES_DIR}" \
        > "${LOG_DIR}/vibrant_${batch_id}.log" 2>&1

    if [[ -d "${out_dir}/VIBRANT_${batch_id}" ]]; then
        echo "[done] ${batch_id}"
    else
        echo "[FAIL] ${batch_id} -- see ${LOG_DIR}/vibrant_${batch_id}.log"
    fi
}
export -f run_one_batch
export VIBRANT_OUT LOG_DIR BATCH_THREADS VIBRANT_DB_DIR VIBRANT_FILES_DIR

find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | sort | \
    parallel -j "${BATCH_JOBS}" run_one_batch {}

n_total=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | wc -l)
n_done=$(find "${VIBRANT_OUT}" -mindepth 2 -maxdepth 2 -type d -name "VIBRANT_batch_*" | wc -l)
log "Step 1 (VIBRANT) finished: ${n_done}/${n_total} batches completed."
log "Next: 03_merge_checkv.sh"