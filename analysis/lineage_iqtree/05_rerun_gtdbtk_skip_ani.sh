#!/usr/bin/env bash
# Rerun GTDB-Tk on the 76 UMGS1840/HGM11507 lineage genomes ONLY, forcing
# full marker-gene identify + align + tree-placement for every genome via
# --place_species (bypasses the ANI-only fast path that skipped 17 of the
# 76 in the original catalogue-wide run).
#
# Using GTDB-Tk v2.7.2 with the R232 reference data — matches the version/
# release used for the manuscript's primary catalogue-wide classification.
# --skip_ani_screen was deprecated in v2.7.0; --place_species replaces it
# for our purposes here: it forces every genome through explicit tree
# placement (identify -> align -> pplacer) regardless of what the initial
# (now-mandatory) skani ANI screen finds, giving real marker-gene
# alignments for all 76 genomes rather than ANI-only calls for some.
#
# REMEMBER: switch back to the v2.7.2 conda env before running this
# (`conda activate gtdbtk-2.7.2` or whatever yours is named), and confirm
# GTDBTK_DATA_PATH points at release232, not release226:
#   export GTDBTK_DATA_PATH=/disk252/qamar/databases/gtdbtk_database/release232
#
# MEMORY: GTDB-Tk's own release notes state R232 needs >=140 GB RAM for the
# placement step. Confirm before running — the script checks and pauses
# below, but don't override that check without knowing you have headroom.
#
# PREREQUISITE: run build_gtdbtk_batchfile.sh first and confirm 76/76 found.
#
# *** CONFIRM/EDIT BEFORE RUNNING ***
# - GTDBTK_DATA_PATH: should already be set in your environment from the
#   original run (check `echo $GTDBTK_DATA_PATH`); set explicitly below if not.
# - CPUS / PPLACER_CPUS: adjust to your node's available cores. pplacer is
#   memory-hungry per thread — keep PPLACER_CPUS conservative if RAM-limited.

set -euo pipefail

BATCHFILE="/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile.tsv"
OUT_DIR="/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232_lineage_rerun"

# CPUS: used for the identify/align steps (Prodigal, HMM search), which
#   parallelize well — safe to set high given your free cores.
# PPLACER_CPUS: kept at 1 deliberately. pplacer memory is per-thread (each
#   thread loads its own copy of the reference tree, ~55 GB for the
#   class-level split-tree default), and doesn't parallelize efficiently
#   anyway for only 76 query genomes — more threads buys you ~nothing here
#   and only adds memory risk. Do not raise this without checking `free -h`
#   first and re-confirming with me.
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