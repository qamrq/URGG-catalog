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
#  SCRIPT  : 02_hocort_host_removal.sh  -  Host and diet read removal (HoCoRT); run: index <name> <fasta> | map [index]
#  AUTHOR  : QRQ
#  INPUT   : QC_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#            REF_DIR/*.fasta.gz (host genome and dietary plant references)
#  OUTPUT  : INDEX_DIR/<name> (Bowtie2 indexes)
#            HOCORT_DIR/<sample>_1_clean.fastq.gz, <sample>_2_clean.fastq.gz
# ==========================================================

source "$(dirname "$0")/config.sh"

mode="${1:-map}"
case "$mode" in
  index)
    name="$2"; fasta="$3"
    ensure_dir "$INDEX_DIR"
    hocort index bowtie2 --input "$fasta" --output "$INDEX_DIR/$name"
    ;;
  map)
    index="$INDEX_DIR/${2:-panhomo_matv1}"
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
  *) echo "Usage: $0 index <name> <fasta> | map [index_name]" >&2; exit 1 ;;
esac
