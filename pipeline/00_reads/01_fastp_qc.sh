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
#  SCRIPT  : 01_fastp_qc.sh  -  Read quality filtering (fastp v1.0.1, default parameters)
#  AUTHOR  : QRQ
#  INPUT   : RAW_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#  OUTPUT  : QC_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#            QC_DIR/fastp_reports/<sample>.html, .json
# ==========================================================

source "$(dirname "$0")/config.sh"
FASTP_THREADS="${FASTP_THREADS:-16}"

ensure_dir "$QC_DIR/fastp_reports"
shopt -s nullglob
r1s=("$RAW_DIR"/*_1.fastq.gz)
[ ${#r1s[@]} -gt 0 ] || { echo "No *_1.fastq.gz files in $RAW_DIR" >&2; exit 1; }

for r1 in "${r1s[@]}"; do
  sample=$(basename "$r1" _1.fastq.gz)
  r2="${r1%_1.fastq.gz}_2.fastq.gz"
  out1="$QC_DIR/${sample}_1.fastq.gz"; out2="$QC_DIR/${sample}_2.fastq.gz"
  [ -s "$out2" ] && continue                                   # already done
  fastp -i "$r1" -I "$r2" -o "$out1" -O "$out2" \
        -h "$QC_DIR/fastp_reports/${sample}.html" \
        -j "$QC_DIR/fastp_reports/${sample}.json" \
        -w "$FASTP_THREADS"
  log_msg "QC complete: $sample"
done
log_msg "fastp QC complete: $QC_DIR"
