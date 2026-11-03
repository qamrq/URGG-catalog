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
#  SCRIPT  : config.sh  -  Shared paths and settings for read processing (steps 01-04)
#  AUTHOR  : QRQ
#  INPUT   : none (edit the paths below, or set them as environment variables)
#  OUTPUT  : variables and helper functions sourced by steps 01-04
# ==========================================================

set -euo pipefail

# Root folder of the read-processing run; every other folder is relative to it.
ROOT="${URGG_READS_ROOT:-/path/to/reads_project}"

RAW_DIR="${RAW_DIR:-$ROOT/raw_fastq}"               # <sample>_1.fastq.gz, <sample>_2.fastq.gz
QC_DIR="${QC_DIR:-$ROOT/qc_fastp}"
REF_DIR="${REF_DIR:-$ROOT/hocort_index}"            # host and diet reference FASTA (.fasta.gz)
INDEX_DIR="${INDEX_DIR:-$REF_DIR/Combined_indexes}"
HOCORT_DIR="${HOCORT_DIR:-$ROOT/hocort_clean_reads}"
ASM_DIR="${ASM_DIR:-$ROOT/megahit_assembly}"
BASALT_DIR="${BASALT_DIR:-$ROOT/basalt_binning}"

THREADS="${THREADS:-48}"

log_msg()    { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
ensure_dir() { mkdir -p "$@"; }
