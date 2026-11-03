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
#  SCRIPT  : 00_install_tools.sh  -  Create conda environments for annotation tools
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : conda environments
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"
export GTDBTK_CONDA_BASE="${GTDBTK_CONDA_BASE:-/path/to/gtdbtk_miniconda3}"

set -euo pipefail
source ./config.sh

source "$(conda info --base)/etc/profile.d/conda.sh"

env_exists () { conda env list | awk '{print $1}' | grep -qx "$1"; }

log "=== Phase 3 tool installation starting ==="

log "[1/8] Adding Prodigal, MMseqs2, DIAMOND, parallel, pigz to '${ENV_DREP}'"
conda activate "${ENV_DREP}"
conda install -y -c bioconda -c conda-forge \
    prodigal \
    "mmseqs2>=15" \
    "diamond>=2.1" \
    parallel \
    pigz
conda deactivate

log "[2/8] GTDB-Tk — skipping install, reusing existing env at ${GTDBTK_ENV_PATH:-<not set>}"

if env_exists "${ENV_EGGNOG}"; then
    log "[3/8] '${ENV_EGGNOG}' already exists — skipping"
else
    log "[3/8] Creating '${ENV_EGGNOG}' (eggNOG-mapper)"
    conda create -y -n "${ENV_EGGNOG}" -c bioconda -c conda-forge "eggnog-mapper>=2.1.12"
fi

if env_exists "${ENV_DBCAN}"; then
    log "[4/8] '${ENV_DBCAN}' already exists — skipping"
else
    log "[4/8] Creating '${ENV_DBCAN}' (dbCAN / run_dbcan)"
    conda create -y -n "${ENV_DBCAN}" -c bioconda -c conda-forge "dbcan>=4.1"
fi

if env_exists "${ENV_AMR}"; then
    log "[5/8] '${ENV_AMR}' already exists — skipping"
else
    log "[5/8] Creating '${ENV_AMR}' (NCBI AMRFinderPlus)"
    conda create -y -n "${ENV_AMR}" -c bioconda -c conda-forge ncbi-amrfinderplus
fi

if env_exists "${ENV_PLATON}"; then
    log "[6/8] '${ENV_PLATON}' already exists — skipping"
else
    log "[6/8] Creating '${ENV_PLATON}' (Platon)"
    conda create -y -n "${ENV_PLATON}" -c bioconda -c conda-forge platon
fi

if env_exists "${ENV_CHECKV}"; then
    log "[7/8] '${ENV_CHECKV}' already exists — skipping"
else
    log "[7/8] Creating '${ENV_CHECKV}' (CheckV)"
    conda create -y -n "${ENV_CHECKV}" -c bioconda -c conda-forge checkv
fi

if env_exists "${ENV_PHAGEBOOST}"; then
    log "[8/8] '${ENV_PHAGEBOOST}' already exists — skipping"
else
    log "[8/8] Creating '${ENV_PHAGEBOOST}' (PhageBoost)"
    conda create -y -n "${ENV_PHAGEBOOST}" -c conda-forge python=3.7 pip
    conda activate "${ENV_PHAGEBOOST}"

    pip install "pyrodigal<3.0" || pip install "pyrodigal==2.3.0"

    pip install phageboost
    conda deactivate
fi

log "=== Tool installation complete ==="
log "Next: run 00b_download_databases.sh to fetch reference databases (large, one-time)."
log "Verify with:  conda env list"
