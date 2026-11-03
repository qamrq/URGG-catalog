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
#  SCRIPT  : 00b_download_databases.sh  -  Download reference databases
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : $URGG_DB/<database>/
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail
source ./config.sh
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! command -v aria2c &>/dev/null; then
    log "Installing aria2c..."
    conda install -y -n base -c conda-forge aria2 2>/dev/null || \
        conda install -y -c conda-forge aria2
fi
command -v aria2c &>/dev/null && log "aria2c ready: $(aria2c --version | head -n1)"

fast_download () {
    local url="$1" outdir="$2" outfile="${3:-$(basename "$1")}"
    mkdir -p "${outdir}"
    if command -v aria2c &>/dev/null; then
        aria2c --continue \
               --max-connection-per-server=16 \
               --split=16 \
               --min-split-size=1M \
               --dir="${outdir}" \
               --out="${outfile}" \
               "${url}" \
        || { log "aria2c failed for ${url}, falling back to wget"; wget -c -O "${outdir}/${outfile}" "${url}"; }
    else
        wget -c -O "${outdir}/${outfile}" "${url}"
    fi
}

mkdir -p "${DB}"/{ncbi_nr,vfdb,eggnog,dbcan,platon,checkv}
cd "${DB}"

log "GTDB-Tk DB: reusing existing R226 install at ${GTDBTK_DB} — skipping download"

log "Downloading + building NCBI-nr DIAMOND db (large — aria2c-accelerated)..."
if [[ ! -f "${NR_DMND}" ]]; then
    fast_download "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" "${DB}/ncbi_nr" "nr.gz"
    cd "${DB}/ncbi_nr"
    conda activate "${ENV_DREP}"
    diamond makedb --in nr.gz -d "${NR_DMND%.dmnd}" --threads "${THREADS}"
    conda deactivate
else
    log "nr.dmnd already present — skipping"
fi

log "Downloading + building VFDB DIAMOND db (aria2c-accelerated)..."
if [[ ! -f "${VFDB_DMND}" ]]; then
    fast_download "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz" "${DB}/vfdb" "VFDB_setB_pro.fas.gz"
    cd "${DB}/vfdb"
    gunzip -kf VFDB_setB_pro.fas.gz
    conda activate "${ENV_DREP}"
    diamond makedb --in VFDB_setB_pro.fas -d "${VFDB_DMND%.dmnd}" --threads "${THREADS}"
    conda deactivate
else
    log "vfdb.dmnd already present — skipping"
fi

log "Downloading eggNOG-mapper data (~50 GB)..."
if [[ ! -f "${EGGNOG_DB}/eggnog.db" ]]; then
    conda activate "${ENV_EGGNOG}"
    mkdir -p "${EGGNOG_DB}"
    download_eggnog_data.py -y --data_dir "${EGGNOG_DB}"
    conda deactivate
else
    log "eggNOG data already present — skipping"
fi

log "Downloading dbCAN database..."
if [[ ! -d "${DBCAN_DB}" ]] || [[ -z "$(ls -A "${DBCAN_DB}" 2>/dev/null)" ]]; then
    conda activate "${ENV_DBCAN}"
    mkdir -p "${DBCAN_DB}"
    run_dbcan database --db_dir "${DBCAN_DB}"
    conda deactivate
else
    log "dbCAN DB already present — skipping"
fi

log "Downloading AMRFinderPlus database..."
conda activate "${ENV_AMR}"
amrfinder -u
conda deactivate

log "Downloading Platon database (aria2c-accelerated)..."
if [[ ! -d "${PLATON_DB}" ]] || [[ -z "$(ls -A "${PLATON_DB}" 2>/dev/null)" ]]; then
    mkdir -p "$(dirname "${PLATON_DB}")"
    fast_download "https://zenodo.org/record/4066768/files/db.tar.gz" "$(dirname "${PLATON_DB}")" "db.tar.gz"
    cd "$(dirname "${PLATON_DB}")"
    tar -xzf db.tar.gz
else
    log "Platon DB already present — skipping"
fi

log "Downloading CheckV database..."
if [[ ! -d "${CHECKV_DB}" ]]; then
    conda activate "${ENV_CHECKV}"
    checkv download_database "${DB}/checkv"
    conda deactivate
else
    log "CheckV DB already present — skipping"
fi

log "=== All database downloads complete ==="
log "If any folder name differs from config.sh, update the *_DB paths there."
