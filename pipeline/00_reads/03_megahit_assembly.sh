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
#  SCRIPT  : 03_megahit_assembly.sh  -  Per-sample metagenome assembly (MEGAHIT v1.2.9)
#  AUTHOR  : QRQ
#  INPUT   : HOCORT_DIR/<sample>_1_clean.fastq.gz, <sample>_2_clean.fastq.gz
#  OUTPUT  : ASM_DIR/<sample>_megahit_output/final.contigs.fa
# ==========================================================

source "$(dirname "$0")/config.sh"
MEGAHIT_MEM="${MEGAHIT_MEM:-0.8}"     # fraction of machine memory

ensure_dir "$ASM_DIR"
shopt -s nullglob
r1s=("$HOCORT_DIR"/*_1_clean.fastq.gz)
[ ${#r1s[@]} -gt 0 ] || { echo "No *_1_clean.fastq.gz files in $HOCORT_DIR" >&2; exit 1; }

for r1 in "${r1s[@]}"; do
  sample=$(basename "$r1" _1_clean.fastq.gz)
  r2="${r1%_1_clean.fastq.gz}_2_clean.fastq.gz"
  out="$ASM_DIR/${sample}_megahit_output"
  [ -f "$out/final.contigs.fa" ] && continue                    # already done
  megahit -1 "$r1" -2 "$r2" \
          --k-list 27,47,67,87,107,127 \
          --min-contig-len 1000 \
          --presets meta-sensitive \
          --memory "$MEGAHIT_MEM" \
          --num-cpu-threads "$THREADS" \
          --continue \
          --out-dir "$out"
  log_msg "Assembly completed: $sample"
done
log_msg "All MEGAHIT assemblies complete: $ASM_DIR"
