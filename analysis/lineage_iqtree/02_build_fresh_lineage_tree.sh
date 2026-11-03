#!/usr/bin/env bash
# Fresh, standalone IQ-TREE build for the UMGS1840/HGM11507 lineage (76 genomes).
#
# IMPORTANT: unlike extend_and_rebuild_tree.sh (previous session), this does
# NOT touch, extend, or reuse the 313-genome / 372-tip tree or its genome
# list. It is a fully independent extract -> trim -> tree cycle using ONLY
# the 76 lineage genomes, run for finer resolution within the lineage itself.
#
# PREREQUISITE: run regenerate_lineage_table.R first, so
# lineage_genome_ids.txt exists.
#
# *** CONFIRM MSA_SOURCE BELOW BEFORE RUNNING ***
# I don't have the path to your GTDB-Tk `align/` output (the concatenated
# bac120 MSA used as input to the original 313-tip tree build). The path
# below assumes GTDB-Tk's standard classify_wf directory layout, as a
# sibling of the `classify/` dir you gave me. Update if your run differs.

set -euo pipefail

GENOME_LIST="/disk252/qamar/metagenome_analysis/annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt"

# <<< CONFIRM/EDIT THIS PATH >>>
MSA_SOURCE="/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232/align/gtdbtk.bac120.user_msa.fasta.gz"

WORKDIR="/disk252/qamar/metagenome_analysis/annotation/iqtree_output/tree_build_lineage_only_fresh"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

N_EXPECTED=$(wc -l < "$GENOME_LIST")
echo "[1/4] Extracting $N_EXPECTED lineage genomes from source MSA..."
echo "      Source: $MSA_SOURCE"

if ! command -v seqkit &> /dev/null; then
  echo "seqkit not found on PATH — falling back to Biopython (see extract_ids_biopython.py)."
  python3 extract_ids_biopython.py "$MSA_SOURCE" "$GENOME_LIST" lineage_raw_msa.fasta
else
  seqkit grep -f "$GENOME_LIST" "$MSA_SOURCE" -o lineage_raw_msa.fasta
fi

N_FOUND=$(grep -c '^>' lineage_raw_msa.fasta || true)
echo "      -> $N_FOUND / $N_EXPECTED genome IDs matched in MSA"
if [ "$N_FOUND" -lt "$N_EXPECTED" ]; then
  echo "      WARNING: fewer matches than expected."
  echo "      GTDB-Tk MSA headers sometimes carry RS_/GB_ prefixes or version"
  echo "      suffixes that differ from the summary's user_genome IDs — check"
  echo "      lineage_raw_msa.fasta headers against lineage_genome_ids.txt."
fi

echo "[2/4] Trimming alignment with trimAl (automated1 heuristic, matches prior pipeline)..."
trimal -in lineage_raw_msa.fasta -out lineage_trimmed_msa.fasta -automated1

echo "[3/4] Running IQ-TREE fresh (ModelFinder + 1000 ultrafast bootstrap)..."
iqtree -s lineage_trimmed_msa.fasta \
  -m MFP \
  -bb 1000 \
  -nt AUTO \
  -pre lineage_only_bac120_tree \
  -redo

echo "[4/4] Done. Fresh 76-tip tree: $WORKDIR/lineage_only_bac120_tree.treefile"
