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
#  SCRIPT  : config_annotation.sh  -  Paths, databases and threads for annotation (5p-7p)
#  AUTHOR  : QRQ
#  INPUT   : none (set URGG_PROJECT, URGG_DB)
#  OUTPUT  : variables used by 5p-7p
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

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
export AMRFINDER_DB="${AMRFINDER_DB:-${URGG_DB}/amrfinderplus/latest}"
export PLATON_DB=${DB}/platon/db
export CHECKV_DB=${DB}/checkv/checkv-db-v1.5

export THREADS=48
export PARALLEL_JOBS=24
export PARALLEL_THREADS=2

# ---- Helper: timestamped logging ----
log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---- Create output dirs up front ----
mkdir -p "${ANNOT}" "${GENES}" "${CATALOG}" "${LOGS}" "${DB}"
