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
#  SCRIPT  : 01_virsorter2.sh  -  Viral sequence detection (VirSorter2)
#  AUTHOR  : QRQ
#  INPUT   : phase_5_virome/batches/batch_*.fa (strain genomes)
#  OUTPUT  : phase_5_virome/01_virsorter2/<batch>/
# ==========================================================

source "$(dirname "$0")/config.sh"
activate_env "${ENV_VIRSORTER2}"

log "Step 1: VirSorter2 on batched MAGs (${BATCH_DIR})"

VS2_CONDA_PREFIX="${DATABASES_DIR}/virsorter2_db/db/conda_envs"

if ! compgen -G "${BATCH_DIR}/batch_*.fa" > /dev/null; then
    log "ERROR: no batch FASTAs found in ${BATCH_DIR}. Run 00c_build_batches.sh first."
    exit 1
fi

if ! find "${VS2_CONDA_PREFIX}" -maxdepth 2 -name "conda-meta" -type d 2>/dev/null | grep -q .; then
    log "Pre-building VirSorter2 internal conda env (serial, one-time)..."
    first_batch=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | sort | head -1)
    virsorter run -w "${VS2_OUT}/_envbuild" -i "${first_batch}" \
        --min-length 1500 -j 4 all \
        --conda-frontend conda --conda-prefix "${VS2_CONDA_PREFIX}" \
        > "${LOG_DIR}/vs2_envbuild.log" 2>&1 || true
    rm -rf "${VS2_OUT}/_envbuild"
    log "Pre-build finished (log: ${LOG_DIR}/vs2_envbuild.log)"
else
    log "VirSorter2 internal conda env already present -- skipping pre-build."
fi

run_one_batch () {
    local batch_fa="$1"
    local batch_id
    batch_id="$(basename "${batch_fa}" .fa)"
    local out_dir="${VS2_OUT}/${batch_id}"

    if [[ -f "${out_dir}/final-viral-score.tsv" ]]; then
        echo "[skip] ${batch_id}"
        return 0
    fi

    # --prep-for-dramv writes for-dramv/ which 04_dramv.sh consumes.
    virsorter run \
        -w "${out_dir}" \
        -i "${batch_fa}" \
        --min-length 1500 \
        --prep-for-dramv \
        -j "${BATCH_THREADS}" \
        all \
        --conda-frontend conda \
        --conda-prefix "${VS2_CONDA_PREFIX}" \
        > "${LOG_DIR}/vs2_${batch_id}.log" 2>&1

    if [[ -f "${out_dir}/final-viral-score.tsv" ]]; then
        echo "[done] ${batch_id}"
    else
        echo "[FAIL] ${batch_id} -- see ${LOG_DIR}/vs2_${batch_id}.log"
    fi
}
export -f run_one_batch
export VS2_OUT LOG_DIR BATCH_THREADS VS2_CONDA_PREFIX

find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | sort | \
    parallel -j "${BATCH_JOBS}" run_one_batch {}

# ---- completion report ----
n_total=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | wc -l)
n_done=$(find "${VS2_OUT}" -mindepth 2 -maxdepth 2 -name "final-viral-score.tsv" | wc -l)
log "Step 1 finished: ${n_done}/${n_total} batches completed."
if [[ "${n_done}" -lt "${n_total}" ]]; then
    log "Incomplete batches (rerun this script to resume):"
    for b in "${BATCH_DIR}"/batch_*.fa; do
        bid="$(basename "${b}" .fa)"
        [[ -f "${VS2_OUT}/${bid}/final-viral-score.tsv" ]] || log "  ${bid}"
    done
fi
log "Next: 02_vibrant.sh, then 03_merge_checkv.sh"
