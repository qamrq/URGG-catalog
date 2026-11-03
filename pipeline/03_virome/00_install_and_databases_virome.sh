#!/usr/bin/env bash
###############################################################################
# 00_install_and_databases.sh — Phase 5: Virome pipeline
#
# NOT meant to be run top-to-bottom as one job. Each section below is
# self-contained (env creation + tool install + DB download for that one
# tool) so you can copy/paste and run sections one at a time, in whatever
# order/pace you want, and verify each before moving to the next.
#
# Env split rationale (discussed earlier): these tools have real dependency
# conflicts (scikit-learn pins, PyTorch, etc.) so they get separate envs
# rather than one shared "virome_analysis" env. VirSorter2+CheckV share an
# env since CheckV is VirSorter2's companion QC tool.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config_virome.sh"

###############################################################################
# SECTION 1: VirSorter2 + CheckV  ->  env: virsorter2_env
###############################################################################

# --- Install ---
conda create -y -n virsorter2_env -c conda-forge -c bioconda python=3.10 virsorter=2 checkv
activate_env virsorter2_env
virsorter config --init-source --db-dir="/disk252/qamar/databases/virsorter2_db/"
# ^ this downloads the VirSorter2 reference DB (~10GB) to the path shown

# --- CheckV database: already present at /disk252/qamar/databases/checkv/checkv-db-v1.5
# (6.4G, genome_db + hmm_db confirmed) -- no download needed. config.sh
# CHECKV_DB path updated to point here directly.

# --- Verify ---
virsorter --version
checkv --version


###############################################################################
# SECTION 2: VIBRANT + MinCED + BLAST  ->  env: vibrant_env
###############################################################################

# --- Install ---
conda create -y -n vibrant_env -c conda-forge -c bioconda vibrant minced blast
conda activate vibrant_env

# --- VIBRANT database ---
# VIBRANT ships a setup script that downloads/annotates its HMM databases.
# Locate it inside the env and run once:
VIBRANT_SETUP=$(find "$(conda info --base)/envs/vibrant_env" -name "download-db.sh" 2>/dev/null | head -1)
bash "${VIBRANT_SETUP}"
# This can take a while (KEGG/Pfam/VOG HMM profiles). Run in tmux/screen.

# --- Verify ---
VIBRANT.py --help | head -5
minced --version
blastn -version


###############################################################################
# SECTION 3: DRAM-v  ->  env: dramv_env
###############################################################################

# --- Install ---
conda create -y -n dramv_env -c conda-forge -c bioconda dram
conda activate dramv_env

# --- DRAM-v database setup (large: KOfam, Pfam, dbCAN, VOGDB, etc. — reuses
# some of the same source DBs as Phase 3's dbCAN/eggNOG steps conceptually,
# but DRAM keeps its own copies) ---
	DRAM-setup.py prepare_databases --output_dir "/disk252/qamar/databases/dramv_config"
# This is the single slowest DB step in this whole pipeline (multiple hours).
# Run in tmux/screen and let it finish before running 05_dramv.sh.

# --- Verify ---
DRAM-setup.py version


###############################################################################
# SECTION 4: PhaGCN2.1  ->  env: phagcn_env  (PyTorch-based — keep isolated)
###############################################################################

# --- Install ---
conda create -y -n phagcn_env python=3.8
conda activate phagcn_env
conda install -y -c pytorch -c nvidia pytorch pytorch-cuda=11.8   # match your GPU/CUDA driver version -- check with `nvidia-smi` first
pip install torch-geometric

git clone https://github.com/KennthShang/PhaGCN2.0.git "${BASE_DIR}/tools/PhaGCN2.0"
cd "${BASE_DIR}/tools/PhaGCN2.0"
pip install -r requirements.txt

# --- Database ---
# PhaGCN2 ships its reference DB inside the repo (database/ dir) -- no
# separate download step, but confirm size/presence after clone:
du -sh "${BASE_DIR}/tools/PhaGCN2.0/database"

# --- Verify ---
python run_Speed_up.py --help


###############################################################################
# SECTION 5: PropagAtE  ->  env: propagate_env
###############################################################################

# --- Install ---
conda create -y -n propagate_env -c conda-forge -c bioconda propagate
conda activate propagate_env

# --- No separate database -- PropagAtE works directly on your BAM files
# (reads mapped to MAG+prophage coordinates) generated in step 9 of the
# pipeline, so nothing to download here. ---

# --- Verify ---
PropagAtE_run.py --version


###############################################################################
# SECTION 6: VOGDB + CARD  (flat database downloads, no conda env needed)
###############################################################################

# --- VOGDB (for prophage-associated gene flagging: integrase/excisionase/
# repressor HMMs -- step 5 of pipeline) ---
mkdir -p "${BASE_DIR}/databases/vogdb"
cd "${BASE_DIR}/databases/vogdb"
wget https://fileshare.lisc.univie.ac.at/vog/latest/vog.hmm.tar.gz
tar -xzf vog.hmm.tar.gz
# produces vog.hmm -- matches VOGDB_HMM in config.sh

# --- CARD (for prophage ARG identification -- step 11 of pipeline) ---
mkdir -p "${BASE_DIR}/databases/card"
cd "${BASE_DIR}/databases/card"
wget https://card.mcmaster.ca/latest/data
tar -xjf data card.json
# produces card.json -- matches CARD_DB in config.sh


###############################################################################
# QUICK REFERENCE: all envs this pipeline expects (from config.sh)
###############################################################################
#   virsorter2_env   -- VirSorter2, CheckV
#   vibrant_env       -- VIBRANT, MinCED, blastn
#   dramv_env         -- DRAM-v
#   phagcn_env        -- PhaGCN2.1 (PyTorch)
#   propagate_env     -- PropagAtE
###############################################################################