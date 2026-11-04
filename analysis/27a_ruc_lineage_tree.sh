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
#  SCRIPT  : 27a_ruc_lineage_tree.sh  -  GTDB-Tk marker placement and IQ-TREE tree of the 76 RUC genomes
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt
#            species genomes
#  OUTPUT  : annotation/gtdbtk_output_v232_lineage_rerun/
#            annotation/iqtree_output/tree_build_lineage_only_fresh/
#  USAGE   : bash 27a_ruc_lineage_tree.sh [batchfile gtdbtk tree]   (no argument: all steps in this order)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

step_batchfile() (
set -euo pipefail

GENOME_LIST="${URGG_PROJECT}/annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt"

# <<< SET THIS >>>
SEARCH_ROOT="${URGG_PROJECT}/drep_analysis/drep_95_species_final/dereplicated_genomes"

OUT_BATCHFILE="${URGG_PROJECT}/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile.tsv"
MISSING_LOG="${URGG_PROJECT}/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile_missing.txt"

mkdir -p "$(dirname "$OUT_BATCHFILE")"
> "$OUT_BATCHFILE"
> "$MISSING_LOG"

# Common bin FASTA extensions, checked in this order per genome ID
EXTS=(fa fasta fna fa.gz fasta.gz fna.gz)

n_found=0
n_missing=0

while IFS= read -r gid; do
  [ -z "$gid" ] && continue
  hit=""
  for ext in "${EXTS[@]}"; do
    candidate=$(find "$SEARCH_ROOT" -maxdepth 6 -type f -name "${gid}.${ext}" 2>/dev/null | head -1)
    if [ -n "$candidate" ]; then
      hit="$candidate"
      break
    fi
  done

  if [ -n "$hit" ]; then
    printf "%s\t%s\n" "$hit" "$gid" >> "$OUT_BATCHFILE"
    n_found=$((n_found + 1))
  else
    echo "$gid" >> "$MISSING_LOG"
    n_missing=$((n_missing + 1))
  fi
done < "$GENOME_LIST"

echo "Matched:  $n_found / $((n_found + n_missing))"
echo "Missing:  $n_missing (see $MISSING_LOG)"
echo "Batchfile written: $OUT_BATCHFILE"

if [ "$n_missing" -gt 0 ]; then
  echo ""
  echo "WARNING: not all 76 genomes were found under SEARCH_ROOT."
  echo "Check $MISSING_LOG and either fix SEARCH_ROOT or the extension list"
  echo "before running GTDB-Tk on this batchfile."
fi
)

step_gtdbtk() (
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
)

step_tree() (
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

command -v seqkit &> /dev/null || { echo "seqkit not found on PATH" >&2; exit 1; }
seqkit grep -f "$GENOME_LIST" "$MSA_SOURCE" -o lineage_raw_msa.fasta

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
)

for s in ${@:-batchfile gtdbtk tree}; do
  echo "==> 27a_ruc_lineage_tree.sh: $s"
  "step_$s"
done
