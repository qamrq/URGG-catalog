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
#  SCRIPT  : config_virome.sh  -  Paths, databases and threads for the virome (8p-13p)
#  AUTHOR  : QRQ
#  INPUT   : none (set URGG_PROJECT, URGG_DB)
#  OUTPUT  : variables used by 8p-13p
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

BASE_DIR="${URGG_PROJECT}"
PHASE5_DIR="${BASE_DIR}/phase_5_virome"

export GENOME_EXT=fa

MAGS_99ANI_DIR="${MAGS_99ANI_DIR:-${URGG_PROJECT}/drep_analysis/drep_99_strain/dereplicated_genomes}"

MAGS_95ANI_DIR="${MAGS_95ANI_DIR:-${URGG_PROJECT}/drep_analysis/drep_95_species_final/dereplicated_genomes}"

MAG_METADATA_TSV="${BASE_DIR}/drep_analysis/metadata.tsv"

SAMPLES_MANIFEST="${BASE_DIR}/phase4_sv/samples_manifest.tsv"

GTDBTK_OUTPUT_DIR="${BASE_DIR}/annotation/gtdbtk_output_v232"

VS2_OUT="${PHASE5_DIR}/01_virsorter2"
CHECKV_OUT="${PHASE5_DIR}/02_checkv"
VIBRANT_OUT="${PHASE5_DIR}/03_vibrant"
CONFIRMED_VIRAL_OUT="${PHASE5_DIR}/04_confirmed_viral"
DRAMV_OUT="${PHASE5_DIR}/05_dramv"
PROPHAGE_FLAG_OUT="${PHASE5_DIR}/06_prophage_flagging"
VOTU_OUT="${PHASE5_DIR}/07_votu_clustering"
PHAGCN_OUT="${PHASE5_DIR}/08_phagcn2"
PROPAGATE_OUT="${PHASE5_DIR}/09_propagate"
CRISPR_OUT="${PHASE5_DIR}/10_crispr_host_linkage"
ARG_OUT="${PHASE5_DIR}/11_prophage_args"
LOG_DIR="${PHASE5_DIR}/logs"

DATABASES_DIR="${URGG_DB}"

CHECKV_DB="${DATABASES_DIR}/checkv/checkv-db-v1.5"
VOGDB_HMM="${DATABASES_DIR}/dramv_config/vog_latest_hmms.txt"
CARD_DB="${DATABASES_DIR}/card/card.json"
DRAMV_CONFIG="${DATABASES_DIR}/dramv_config"
PHAGCN_DB=""

BATCH_SIZE="${BATCH_SIZE:-200}"
BATCH_JOBS="${BATCH_JOBS:-6}"
BATCH_THREADS="${BATCH_THREADS:-8}"

BATCH_DIR="${PHASE5_DIR}/00_batches"
MAG_INDEX_TSV="${BATCH_DIR}/mag_index.tsv"

THREADS="${THREADS:-48}"
PARALLEL_JOBS=24
PARALLEL_THREADS=2

log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "${BATCH_DIR}" "${VS2_OUT}" "${CHECKV_OUT}" "${VIBRANT_OUT}" "${CONFIRMED_VIRAL_OUT}" \
         "${DRAMV_OUT}" "${PROPHAGE_FLAG_OUT}" "${VOTU_OUT}" "${PHAGCN_OUT}" \
         "${PROPAGATE_OUT}" "${CRISPR_OUT}" "${ARG_OUT}" "${LOG_DIR}"
