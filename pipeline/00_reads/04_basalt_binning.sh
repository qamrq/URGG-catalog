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
#  SCRIPT  : 04_basalt_binning.sh  -  Binning and refinement (BASALT v1.2.2 autobinning)
#  AUTHOR  : QRQ
#  INPUT   : ASM_DIR/*_megahit_output/final.contigs.fa (step 03)
#            HOCORT_DIR/*_1_clean.fastq.gz, *_2_clean.fastq.gz (step 02)
#  OUTPUT  : BASALT_DIR/ (bins, quality reports, best bin set)
# ==========================================================

source "$(dirname "$0")/config.sh"
BASALT_MEM_GB="${BASALT_MEM_GB:-500}"

ensure_dir "$BASALT_DIR"
shopt -s nullglob

# 1) Link each MEGAHIT assembly as <sample>_final_contigs.fa
for d in "$ASM_DIR"/*_megahit_output; do
  sample=$(basename "$d" _megahit_output)
  ln -sfn "$d/final.contigs.fa" "$BASALT_DIR/${sample}_final_contigs.fa"
done

# 2) Link clean reads as <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
for r1 in "$HOCORT_DIR"/*_1_clean.fastq.gz; do
  sample=$(basename "$r1" _1_clean.fastq.gz)
  ln -sfn "$r1" "$BASALT_DIR/${sample}_R1.fastq.gz"
  ln -sfn "${r1%_1_clean.fastq.gz}_2_clean.fastq.gz" "$BASALT_DIR/${sample}_R2.fastq.gz"
done

# 3) Assembly list (comma-separated) and read pairs "R1,R2/R1,R2/..."
cd "$BASALT_DIR"
asm=$(ls -1 *_final_contigs.fa | sort | paste -sd, -)
pairs=$(ls -1 *_R1.fastq.gz | sort | sed 's/\(.*\)_R1\.fastq\.gz/\1_R1.fastq.gz,\1_R2.fastq.gz/' | paste -sd/ -)
echo "$asm"   > assemblies.txt
echo "$pairs" > shortreads.txt

# 4) Binning
BASALT -a "$asm" -s "$pairs" --module autobinning --sensitive quick \
       -t "$THREADS" -m "$BASALT_MEM_GB" --min-cpn 70 --max-ctn 10
log_msg "BASALT autobinning complete: $BASALT_DIR"
