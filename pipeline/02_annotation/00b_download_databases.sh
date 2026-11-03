#!/usr/bin/env bash
# ============================================================
# 00b_download_databases.sh  —  One-time download of remaining reference DBs
# GTDB-Tk is SKIPPED (existing R226 install reused — see config.sh).
#
# These are large and can run now in the background (tmux/screen). Total
# footprint for what's left to download:
#   NCBI-nr (DIAMOND) ~150-200 GB build  |  eggNOG ~50 GB
#   dbCAN ~ a few GB   |  CheckV ~1.5 GB  |  Platon ~1.5 GB  |  VFDB small
# Make sure /disk252/qamar/databases has room (budget ~250-300 GB).
#
# Comment out any block you already have.
# ============================================================
set -euo pipefail
source ./config.sh
source "$(conda info --base)/etc/profile.d/conda.sh"

# ------------------------------------------------------------
# aria2c — multi-connection downloader, big speedup for large single files
# (nr.gz ~150GB, Platon DB, VFDB). eggNOG/dbCAN/CheckV use their own Python
# downloader scripts internally, so aria2c isn't applicable there.
# ------------------------------------------------------------
if ! command -v aria2c &>/dev/null; then
    log "Installing aria2c..."
    conda install -y -n base -c conda-forge aria2 2>/dev/null || \
        conda install -y -c conda-forge aria2   # falls back to active env if base install isn't permitted
fi
command -v aria2c &>/dev/null && log "aria2c ready: $(aria2c --version | head -n1)"

# fast_download <url> <output_dir> [output_filename]
# Uses aria2c with 16 parallel connections/splits; falls back to wget -c if
# aria2c is unavailable or the server doesn't support ranged requests.
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

# ------------------------------------------------------------
# GTDB-Tk reference — SKIPPED: reusing existing R226 DB at
# /disk212/qamar/databases/release226/ (confirmed present, matched to
# installed GTDB-Tk v2.5.2). No download needed.
# ------------------------------------------------------------
log "GTDB-Tk DB: reusing existing R226 install at ${GTDBTK_DB} — skipping download"

# ------------------------------------------------------------
# NCBI-nr -> DIAMOND db  (VERY large & slow; see note in 04a about cost)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# VFDB (core set) -> DIAMOND db
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# eggNOG-mapper data (~50 GB)
# ------------------------------------------------------------
log "Downloading eggNOG-mapper data (~50 GB)..."
if [[ ! -f "${EGGNOG_DB}/eggnog.db" ]]; then
    conda activate "${ENV_EGGNOG}"
    mkdir -p "${EGGNOG_DB}"
    download_eggnog_data.py -y --data_dir "${EGGNOG_DB}"
    conda deactivate
else
    log "eggNOG data already present — skipping"
fi

# ------------------------------------------------------------
# dbCAN database
# NOTE: dbCAN v5.x (confirmed: 5.2.9) uses `run_dbcan database`, NOT
# `dbcan_build` (that was a v4.x-only command name). The CLI is now a
# subcommand structure (run_dbcan CAZyme_annotation, run_dbcan database, etc).
# ------------------------------------------------------------
log "Downloading dbCAN database..."
if [[ ! -d "${DBCAN_DB}" ]] || [[ -z "$(ls -A "${DBCAN_DB}" 2>/dev/null)" ]]; then
    conda activate "${ENV_DBCAN}"
    mkdir -p "${DBCAN_DB}"
    run_dbcan database --db_dir "${DBCAN_DB}"
    conda deactivate
else
    log "dbCAN DB already present — skipping"
fi

# ------------------------------------------------------------
# AMRFinderPlus database
# NOTE: `amrfinder -u` ignores the AMRFINDER_DB env var and always installs
# to <env>/share/amrfinderplus/data/<version>/. After this runs, copy the
# resulting version-dated path into AMRFINDER_DB in config.sh so
# 05_amrfinderplus.sh points at it explicitly (already done for
# 2026-05-15.1 — update if a newer version downloads here).
# ------------------------------------------------------------
log "Downloading AMRFinderPlus database..."
conda activate "${ENV_AMR}"
amrfinder -u    # downloads/updates the latest AMR database to its own default path
conda deactivate

# ------------------------------------------------------------
# Platon database
# ------------------------------------------------------------
log "Downloading Platon database (aria2c-accelerated)..."
if [[ ! -d "${PLATON_DB}" ]] || [[ -z "$(ls -A "${PLATON_DB}" 2>/dev/null)" ]]; then
    mkdir -p "$(dirname "${PLATON_DB}")"
    fast_download "https://zenodo.org/record/4066768/files/db.tar.gz" "$(dirname "${PLATON_DB}")" "db.tar.gz"
    cd "$(dirname "${PLATON_DB}")"
    tar -xzf db.tar.gz     # unpacks to db/
else
    log "Platon DB already present — skipping"
fi

# ------------------------------------------------------------
# CheckV database
# ------------------------------------------------------------
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