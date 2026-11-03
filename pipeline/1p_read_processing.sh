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
#  SCRIPT  : 1p_read_processing.sh  -  Read QC, host/diet read removal, assembly and binning
#  AUTHOR  : QRQ
#  INPUT   : raw reads: URGG_READS_ROOT/raw_fastq/<sample>_1/_2.fastq.gz
#            host and diet references (hostremoval)
#  OUTPUT  : clean reads, MEGAHIT assemblies, BASALT bins
#  USAGE   : bash 1p_read_processing.sh [qc hostremoval assembly binning]   (no argument: all steps in this order)
# ==========================================================

set -euo pipefail

step_qc() (
source "$(dirname "$0")/config_reads.sh"
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
)

step_hostremoval() (
source "$(dirname "$0")/config_reads.sh"

mode="${HOCORT_MODE:-map}"
case "$mode" in
  index)
    name="${HOCORT_INDEX_NAME:?set HOCORT_INDEX_NAME}"; fasta="${HOCORT_FASTA:?set HOCORT_FASTA}"
    ensure_dir "$INDEX_DIR"
    hocort index bowtie2 --input "$fasta" --output "$INDEX_DIR/$name"
    ;;
  map)
    index="$INDEX_DIR/${HOCORT_INDEX:-panhomo_matv1}"
    ensure_dir "$HOCORT_DIR/log_files"
    shopt -s nullglob
    r1s=("$QC_DIR"/*_1.fastq.gz)
    [ ${#r1s[@]} -gt 0 ] || { echo "No *_1.fastq.gz files in $QC_DIR" >&2; exit 1; }
    for r1 in "${r1s[@]}"; do
      sample=$(basename "$r1" _1.fastq.gz)
      r2="${r1%_1.fastq.gz}_2.fastq.gz"
      out1="$HOCORT_DIR/${sample}_1_clean.fastq.gz"; out2="$HOCORT_DIR/${sample}_2_clean.fastq.gz"
      [ -s "$out2" ] && continue                               # already done
      hocort map bowtie2 -x "$index" -i "$r1" "$r2" --filter true -o "$out1" "$out2" \
             preset='end-to-end' \
             --log-file="$HOCORT_DIR/log_files/${sample}_map.log"
      log_msg "HoCoRT complete: $sample"
    done
    log_msg "HoCoRT host removal complete: $HOCORT_DIR"
    ;;
  *) echo "HOCORT_MODE must be index or map" >&2; exit 1 ;;
esac
)

step_assembly() (
source "$(dirname "$0")/config_reads.sh"
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
)

step_binning() (
source "$(dirname "$0")/config_reads.sh"
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
)

for s in ${@:-qc hostremoval assembly binning}; do
  echo "==> 1p_read_processing.sh: $s"
  "step_$s"
done
