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
#  SCRIPT  : 16a_run_trait_models.sh  -  Run 15a for every trait family and adjustment set
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s54, s50, s36 outputs
#  OUTPUT  : annotation/results_phase3/s49_trait_presence_taxonomy_glmm/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail
SCR=${SCR:-$(cd "$(dirname "$0")" && pwd)}
R3=${URGG_PROJECT}/annotation/results_phase3
S54=$R3/s54_cazyme_propagated
S50=$R3/s50_subpathway_propagated
S36=$R3/s36_cooperative_ko_core_specific
TMP=$R3/s49_trait_presence_taxonomy_glmm/inputs; mkdir -p "$TMP"
CORES=${N_CORES:-24}
run49 () {
  echo -e "\n######## $2  [$3] ########"
  GENOME_COL=$4 TRAIT_COL=$5 N_CORES=$CORES ADJ_KEEP=$3 Rscript "$SCR/15a_trait_presence_glmm.R" "$1" "$2"
}
grep -q '"compl"\|compl =' "$SCR/15a_trait_presence_glmm.R" || { echo "ERROR: use the updated 15a script (with the 'compl' level) to $SCR first"; exit 1; }

run49 "$S54/54_cazy_presence_for_glmm.tsv" cazy_prop unadjusted,full genome cazy_family

# 1. Small trait sets: completeness-only AND full, side by side
awk -F'\t' 'NR==1 || $2 ~ /^(K03147|K00788|K00949)$/' "$S50/50_marker_ko_presence_propagated.tsv" > "$TMP/thiamine_propagated.tsv"
run49 "$TMP/thiamine_propagated.tsv" thiamine_prop_adjcompare compl,full genome_id ko
for t in 50 75 100; do
  run49 "$S50/subpathway_presence_propagated_${t}.tsv" subpathways_prop_${t}_adjcompare compl,full genome_id subpathway
done
awk -F',' 'NR==1 || $3=="TRUE"' "$S36/table_cooperative_category_presence_per_genome.csv" > "$TMP/cooperative_categories_present.csv"
run49 "$TMP/cooperative_categories_present.csv" coop_categories_adjcompare unadjusted,size,compl,full genome cooperative_categories
run49 "$S54/54_cazy_substrate_presence_propagated.tsv" cazy_substrates_adjcompare compl,full genome substrate

run49 "$S54/54_cazy_presence_for_glmm.tsv" cazy_prop_compl compl genome cazy_family

echo -e "\nAll done. Now: Rscript $SCR/33a_manuscript_numbers.R"
