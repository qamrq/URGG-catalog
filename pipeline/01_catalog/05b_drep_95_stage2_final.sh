#!/usr/bin/env bash
###############################################################################
# 05b_drep_95_stage2_final.sh
#
# STAGE 2 of hierarchical dereplication for the 95% species catalog.
#
# Takes all per-subgroup representatives produced by Stage 1 (05a), pools
# them, and runs one final 95% dRep to collapse redundancy ACROSS hosts/
# BioProjects (e.g. the same species recovered in cattle and yak becomes
# one representative here).
#
# The pooled Stage-1 winner set is far smaller than the original 117,145
# genomes, so this final run fits comfortably under the memory ceiling.
#
# Output = the final species-level catalog:
#   - read mapping reference (Bowtie2)
#   - abundance (CoverM)
#   - GTDB-Tk taxonomy
#   - gene catalog + DIAMOND annotation
###############################################################################
set -euo pipefail

# ── EDIT THESE ────────────────────────────────────────────────────────────────
WORKDIR="/disk252/qamar/metagenome_analysis/drep_analysis"
THREADS=48
# ─────────────────────────────────────────────────────────────────────────────

STAGE1_DIR="$WORKDIR/stage1_per_host"
GENOMEINFO="$WORKDIR/genomeInfo.csv"
POOLED_LIST="$WORKDIR/stage2_pooled_reps.txt"
OUTPUT="$WORKDIR/drep_95_species_final"
LOGFILE="$WORKDIR/drep_95_stage2_run.log"

# ── Step 1: collect all Stage-1 representatives into one list ────────────────
echo ">> Collecting Stage-1 representatives..."
: > "$POOLED_LIST"

n_subgroups=0
while IFS=$'\t' read -r subgroup list; do
    repdir="$STAGE1_DIR/${subgroup}_drep/dereplicated_genomes"
    if [[ -d "$repdir" ]]; then
        ls "$repdir"/*.fa 2>/dev/null >> "$POOLED_LIST" || true
        n_subgroups=$((n_subgroups+1))
    else
        echo "   WARN: no representatives found for subgroup '$subgroup' at $repdir" >&2
    fi
done < "$STAGE1_DIR/subgroup_index.txt"

n_pooled=$(wc -l < "$POOLED_LIST")
echo "   Pooled $n_pooled representatives from $n_subgroups subgroups"
echo "   Pooled list: $POOLED_LIST"
echo ""

if [[ "$n_pooled" -eq 0 ]]; then
    echo "ERROR: no Stage-1 representatives found. Did Stage 1 complete?" >&2
    exit 1
fi

# ── Step 2: guard against existing output dir ────────────────────────────────
if [[ -d "$OUTPUT" ]]; then
    echo "ERROR: output directory exists: $OUTPUT" >&2
    echo "       To start fresh: rm -rf $OUTPUT" >&2
    exit 1
fi

echo "========================================================"
echo "STAGE 2: final cross-host dereplication (95% ANI)"
echo "Input (pooled Stage-1 reps): $n_pooled genomes"
echo "Threads: $THREADS"
echo "Output: $OUTPUT"
echo "Started: $(date)"
echo "========================================================"

/usr/bin/time -v dRep dereplicate "$OUTPUT" \
    -p $THREADS \
    --genomeInfo "$GENOMEINFO" \
    --ignoreGenomeQuality \
    -pa 0.9 -sa 0.95 -nc 0.30 \
    -cm larger \
    --S_algorithm skani \
    --multiround_primary_clustering \
    --primary_chunksize 2000 \
    -g "$POOLED_LIST" \
    2>&1 | tee "$LOGFILE"

echo ""
echo "========================================================"
echo "STAGE 2 complete: $(date)"
echo "FINAL species catalog: $OUTPUT/dereplicated_genomes/"
echo "Final rep count: $(ls $OUTPUT/dereplicated_genomes/ | wc -l)"
echo ""
echo "This is your final 95% ANI species-level catalog."
echo "Key tables:"
echo "  $OUTPUT/data_tables/Cdb.csv   final cluster assignments"
echo "  $OUTPUT/data_tables/Wdb.csv   final representatives"
echo ""
echo "NOTE on two-stage cluster tracing: to map an ORIGINAL genome to its"
echo "FINAL species cluster, you must chain Stage-1 Cdb -> Stage-2 Cdb, since"
echo "Stage 2 only saw Stage-1 winners. A helper for this can be written"
echo "later if needed for the shared-vs-host-specific analysis."
echo "========================================================"
