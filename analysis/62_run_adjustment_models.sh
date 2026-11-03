#!/usr/bin/env bash
# ============================================================
# 62_run_adjustment_models.sh
#
# Runs 49 for every trait family with the adjustment sets needed to decide
# the primary model (item 1 of the gap list):
#   compl = taxonomy + MAG completeness            (NEW, no genome size)
#   full  = taxonomy + log genome size + completeness (current primary)
# Also restores the one-tool CAZy file overwritten by the two-tool run.
# Each run has its own label, so nothing is overwritten.
#
#   bash 62_run_adjustment_models.sh 2>&1 | tee 62_run.log
# Needs the updated 49_trait_presence_taxonomy_glmm.R (ADJ_KEEP accepts "compl").
# ============================================================
set -euo pipefail
SCR=${SCR:-/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts}
R3=/disk252/qamar/metagenome_analysis/annotation/results_phase3
S54=$R3/s54_cazyme_propagated
S50=$R3/s50_subpathway_propagated
S36=$R3/s36_cooperative_ko_core_specific
TMP=$R3/s49_trait_presence_taxonomy_glmm/inputs; mkdir -p "$TMP"
CORES=${N_CORES:-24}
run49 () {  # $1 table  $2 label  $3 ADJ_KEEP  $4 GENOME_COL  $5 TRAIT_COL
  echo -e "\n######## $2  [$3] ########"
  GENOME_COL=$4 TRAIT_COL=$5 N_CORES=$CORES ADJ_KEEP=$3 Rscript "$SCR/49_trait_presence_taxonomy_glmm.R" "$1" "$2"
}
grep -q '"compl"\|compl =' "$SCR/49_trait_presence_taxonomy_glmm.R" || { echo "ERROR: copy the updated 49 script (with the 'compl' level) to $SCR first"; exit 1; }

# 0. Restore the primary one-tool CAZy file (overwritten by the two-tool run)
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

# 2. CAZy families: completeness-only (the full model is in cazy_prop from step 0)
run49 "$S54/54_cazy_presence_for_glmm.tsv" cazy_prop_compl compl genome cazy_family

echo -e "\nAll done. Now: Rscript $SCR/61_pull_manuscript_numbers.R"
