#!/usr/bin/env bash
# Build a GTDB-Tk batchfile (TSV: <fasta_path>\t<genome_id>) for the 76
# lineage genomes, by searching under SEARCH_ROOT for a file matching each
# genome ID + a set of common bin FASTA extensions.
#
# *** SET SEARCH_ROOT BEFORE RUNNING *** — point it at wherever your bin
# FASTAs live (can be a broad parent directory; this searches recursively).

set -euo pipefail

GENOME_LIST="/disk252/qamar/metagenome_analysis/annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt"

# <<< SET THIS >>>
SEARCH_ROOT="/disk252/qamar/metagenome_analysis/drep_analysis/drep_95_species_final/dereplicated_genomes"

OUT_BATCHFILE="/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile.tsv"
MISSING_LOG="/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232_lineage_rerun/lineage_76_batchfile_missing.txt"

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
    # try exact "<id>.<ext>" first (fast path), then fall back to a
    # recursive find in case directory nesting/prefixing differs
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
