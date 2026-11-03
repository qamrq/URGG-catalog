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
#  SCRIPT  : 00_install_and_databases_virome.sh  -  Install viral tools and download their databases
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : conda environments, $URGG_DB/<tool databases>
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"

source "$(dirname "$0")/config.sh"

# --- Install ---
conda create -y -n virsorter2_env -c conda-forge -c bioconda python=3.10 virsorter=2 checkv
activate_env virsorter2_env
virsorter config --init-source --db-dir="${URGG_DB}/virsorter2_db/"

# --- Verify ---
virsorter --version
checkv --version

# --- Install ---
conda create -y -n vibrant_env -c conda-forge -c bioconda vibrant minced blast
conda activate vibrant_env

VIBRANT_SETUP=$(find "$(conda info --base)/envs/vibrant_env" -name "download-db.sh" 2>/dev/null | head -1)
bash "${VIBRANT_SETUP}"

# --- Verify ---
VIBRANT.py --help | head -5
minced --version
blastn -version

# --- Install ---
conda create -y -n dramv_env -c conda-forge -c bioconda dram
conda activate dramv_env

	DRAM-setup.py prepare_databases --output_dir "${URGG_DB}/dramv_config"

# --- Verify ---
DRAM-setup.py version

# --- Install ---
conda create -y -n phagcn_env python=3.8
conda activate phagcn_env
conda install -y -c pytorch -c nvidia pytorch pytorch-cuda=11.8
pip install torch-geometric

git clone https://github.com/KennthShang/PhaGCN2.0.git "${BASE_DIR}/tools/PhaGCN2.0"
cd "${BASE_DIR}/tools/PhaGCN2.0"
pip install -r requirements.txt

du -sh "${BASE_DIR}/tools/PhaGCN2.0/database"

# --- Verify ---
python run_Speed_up.py --help

# --- Install ---
conda create -y -n propagate_env -c conda-forge -c bioconda propagate
conda activate propagate_env

# --- Verify ---
PropagAtE_run.py --version

mkdir -p "${BASE_DIR}/databases/vogdb"
cd "${BASE_DIR}/databases/vogdb"
wget https://fileshare.lisc.univie.ac.at/vog/latest/vog.hmm.tar.gz
tar -xzf vog.hmm.tar.gz
# produces vog.hmm -- matches VOGDB_HMM in config.sh

mkdir -p "${BASE_DIR}/databases/card"
cd "${BASE_DIR}/databases/card"
wget https://card.mcmaster.ca/latest/data
tar -xjf data card.json
# produces card.json -- matches CARD_DB in config.sh

