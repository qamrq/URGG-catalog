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
#  SCRIPT  : 06_build_fresh_lineage_tree.sh  -  RUC lineage tree from the rerun alignment (trimAl, IQ-TREE)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt
#            annotation/gtdbtk_output_v232_lineage_rerun/align/gtdbtk.bac120.user_msa.fasta.gz
#  OUTPUT  : annotation/iqtree_output/tree_build_lineage_only_fresh/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

GENOME_LIST="${URGG_PROJECT}/annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt"

MSA_SOURCE="${URGG_PROJECT}/annotation/gtdbtk_output_v232_lineage_rerun/align/gtdbtk.bac120.user_msa.fasta.gz"

N_EXPECTED_HINT="all 76 should match now — if not 76/76, stop and check the classification_method column in the rerun's gtdbtk.bac120.summary.tsv before trusting the tree"

WORKDIR="${URGG_PROJECT}/annotation/iqtree_output/tree_build_lineage_only_fresh"
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
