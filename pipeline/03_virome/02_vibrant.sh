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
#  SCRIPT  : 02_vibrant.sh  -  Viral sequence detection (VIBRANT)
#  AUTHOR  : QRQ
#  INPUT   : phase_5_virome/batches/batch_*.fa
#  OUTPUT  : phase_5_virome/03_vibrant/<batch>/
# ==========================================================

source "$(dirname "$0")/config.sh"
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
