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
#  SCRIPT  : config.sh  -  Shared paths, databases and threads for annotation steps
#  AUTHOR  : QRQ
#  INPUT   : none (set URGG_PROJECT, URGG_DB)
#  OUTPUT  : variables used by every annotation step
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"
export CONDA_BASE_MAIN="${CONDA_BASE_MAIN:-/path/to/miniconda3}"
export GTDBTK_CONDA_BASE="${GTDBTK_CONDA_BASE:-/path/to/gtdbtk_miniconda3}"

export PROJECT=${URGG_PROJECT}
export DREP_DIR=${PROJECT}/drep_analysis
export METADATA=${DREP_DIR}/metadata.tsv
export GENOME_EXT=fa

export GENOME_DIR=${DREP_DIR}/drep_95_species_final/dereplicated_genomes

export SV_STRAIN_DIR=${DREP_DIR}/drep_99_strain/dereplicated_genomes

# ---- Output roots ----
export ANNOT=${PROJECT}/annotation
export GENES=${ANNOT}/genes
export CATALOG=${ANNOT}/gene_catalog
export LOGS=${ANNOT}/logs
export DB=${URGG_DB}

export GTDBTK_DB=${URGG_DB}/gtdbtk_database/release232/
export NR_DMND=${DB}/ncbi_nr/nr.dmnd
export VFDB_DMND=${DB}/vfdb/vfdb.dmnd
export EGGNOG_DB=${DB}/eggnog
export DBCAN_DB=${DB}/dbcan
export AMRFINDER_DB=${CONDA_BASE_MAIN}/envs/amrfinder/share/amrfinderplus/data/latest
export PLATON_DB=${DB}/platon/db
export CHECKV_DB=${DB}/checkv/checkv-db-v1.5

export THREADS=48
export PARALLEL_JOBS=24
export PARALLEL_THREADS=2

# ---- Conda envs (created by 00_install_tools.sh) ----
export ENV_DREP=drep

export GTDBTK_ENV_PATH=${GTDBTK_CONDA_BASE}/envs/gtdbtk-2.5.2

export ENV_EGGNOG=eggnog
export ENV_DBCAN=dbcan
export ENV_AMR=amrfinder
export ENV_PLATON=platon
export ENV_CHECKV=checkv
export ENV_PHAGEBOOST=phageboost

# ---- Helper: activate a conda env non-interactively ----
activate_env () {
    local env_name="$1"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

# ---- Helper: timestamped logging ----
log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---- Create output dirs up front ----
mkdir -p "${ANNOT}" "${GENES}" "${CATALOG}" "${LOGS}" "${DB}"
