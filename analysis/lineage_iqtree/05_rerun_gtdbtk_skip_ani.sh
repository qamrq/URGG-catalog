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
#  SCRIPT  : 05_rerun_gtdbtk_skip_ani.sh  -  GTDB-Tk marker placement for the 76 RUC genomes
#  AUTHOR  : QRQ
#  INPUT   : annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile.tsv
#  OUTPUT  : annotation/gtdbtk_output_v232_lineage_rerun/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

BATCHFILE="${URGG_PROJECT}/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile.tsv"
OUT_DIR="${URGG_PROJECT}/annotation/gtdbtk_output_v232_lineage_rerun"

CPUS=64
PPLACER_CPUS=4

# <<< CONFIRM THIS IS SET, or hardcode your reference data path >>>
: "${GTDBTK_DATA_PATH:?GTDBTK_DATA_PATH is not set — export it (same path used for the original catalogue-wide run) before running this script.}"

echo "Using GTDB-Tk reference data at: $GTDBTK_DATA_PATH"
echo "Batchfile: $BATCHFILE ($(wc -l < "$BATCHFILE") genomes)"
echo ""
echo "NOTE: GTDB-Tk's own release notes for v2.7.0+ (R232 data) specify a"
echo "minimum of ~140 GB RAM for the divide-and-conquer tree placement step,"
echo "higher than earlier releases needed. Confirm you actually have that"
echo "much free before proceeding:"
free -h
echo ""
read -p "Press Enter to continue if the above looks sufficient, or Ctrl+C to abort... "

gtdbtk classify_wf \
  --batchfile "$BATCHFILE" \
  --out_dir "$OUT_DIR" \
  --place_species \
  --cpus "$CPUS" \
  --pplacer_cpus "$PPLACER_CPUS" \
  --extension fa

echo ""
echo "Done. Check for a complete 76-genome alignment at:"
echo "  $OUT_DIR/align/gtdbtk.bac120.user_msa.fasta.gz"
echo ""
echo "Verify count before proceeding:"
echo "  zcat $OUT_DIR/align/gtdbtk.bac120.user_msa.fasta.gz | grep -c '^>'"
echo "(should be 76, or close to it minus any true classify-tree failures —"
echo " check $OUT_DIR/classify/gtdbtk.bac120.summary.tsv 'classification_method'"
echo " column, should now read something other than ani_screen for all 76.)"
