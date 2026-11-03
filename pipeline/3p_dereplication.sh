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
#  SCRIPT  : 3p_dereplication.sh  -  Two-stage species clustering (95% ANI) and strain clustering (99% ANI) with dRep
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/genomes/, metadata.tsv, genomeInfo.csv, genome lists
#  OUTPUT  : drep_analysis/drep_95_species_final/ (25,451 species)
#            drep_analysis/drep_99_strain/ (23,999 strains)
#  USAGE   : bash 3p_dereplication.sh [species_stage1 species_stage2 strains]   (no argument: all steps in this order)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

step_species_stage1() (
set -euo pipefail

WORKDIR="${URGG_PROJECT}/drep_analysis"

CONCURRENT_JOBS=4
THREADS_PER_JOB=24

MAX_PER_GROUP=10000

METADATA="$WORKDIR/metadata.tsv"
GENOMEINFO="$WORKDIR/genomeInfo.csv"
STAGE1_DIR="$WORKDIR/stage1_per_host"
LISTS_DIR="$STAGE1_DIR/genome_lists"
GENOMES_DIR="$WORKDIR/genomes"

mkdir -p "$STAGE1_DIR" "$LISTS_DIR"

# check GNU parallel is available
if ! command -v parallel &>/dev/null; then
    echo "ERROR: GNU parallel not found. Install with:" >&2
    exit 1
fi

echo "========================================================"
echo "STAGE 1 (PARALLEL): per-host dereplication (95% ANI)"
echo "Concurrent jobs: $CONCURRENT_JOBS  x  threads/job: $THREADS_PER_JOB"
echo "  = up to $((CONCURRENT_JOBS * THREADS_PER_JOB)) cores in use"
echo "Max genomes/subgroup: $MAX_PER_GROUP"
echo "Started: $(date)"
echo "========================================================"

> "$STAGE1_DIR/subgroup_index.txt"
hosts=$(awk -F'\t' 'NR>1 {print $3}' "$METADATA" | sort -u)

for host in $hosts; do
    host_count=$(awk -F'\t' -v h="$host" 'NR>1 && $3==h' "$METADATA" | wc -l)
    safe_host=$(echo "$host" | sed 's/[^A-Za-z0-9._-]/_/g')

    if [[ "$host_count" -le "$MAX_PER_GROUP" ]]; then
        list="$LISTS_DIR/host__${safe_host}.txt"
        awk -F'\t' -v h="$host" 'NR>1 && $3==h {print $1}' "$METADATA" \
            | sed "s|^|$GENOMES_DIR/|" > "$list"
        echo -e "host__${safe_host}\t${list}" >> "$STAGE1_DIR/subgroup_index.txt"
    else
        bps=$(awk -F'\t' -v h="$host" 'NR>1 && $3==h {print $2}' "$METADATA" | sort -u)
        for bp in $bps; do
            safe_bp=$(echo "$bp" | sed 's/[^A-Za-z0-9._-]/_/g')
            list="$LISTS_DIR/host__${safe_host}__bp__${safe_bp}.txt"
            awk -F'\t' -v h="$host" -v b="$bp" \
                'NR>1 && $3==h && $2==b {print $1}' "$METADATA" \
                | sed "s|^|$GENOMES_DIR/|" > "$list"
            echo -e "host__${safe_host}__bp__${safe_bp}\t${list}" \
                >> "$STAGE1_DIR/subgroup_index.txt"
        done
    fi
done

n_subgroups=$(wc -l < "$STAGE1_DIR/subgroup_index.txt")
echo ""
echo ">> $n_subgroups subgroups created. Sizes:"
while IFS=$'\t' read -r sg list; do
    printf "   %-45s %d genomes\n" "$sg" "$(wc -l < "$list")"
done < "$STAGE1_DIR/subgroup_index.txt" | sort -t' ' -k2 -rn
echo ""

run_one_subgroup() {
    local subgroup="$1"
    local list="$2"
    local out="$STAGE1_DIR/${subgroup}_drep"
    local n; n=$(wc -l < "$list")

    # skip if already complete
    if [[ -d "$out/dereplicated_genomes" ]] && \
       [[ $(ls "$out/dereplicated_genomes/" 2>/dev/null | wc -l) -gt 0 ]]; then
        echo "[SKIP] $subgroup (already done, $(ls "$out/dereplicated_genomes/" | wc -l) reps)"
        return 0
    fi

    [[ -d "$out" ]] && rm -rf "$out"

    if [[ "$n" -le 1 ]]; then
        mkdir -p "$out/dereplicated_genomes"
        [[ "$n" -eq 1 ]] && cp "$(cat "$list")" "$out/dereplicated_genomes/"
        echo "[COPY] $subgroup ($n genome, no dRep needed)"
        return 0
    fi

    echo "[START] $subgroup ($n genomes) @ $(date +%H:%M:%S)"
    if dRep dereplicate "$out" \
        -p "$THREADS_PER_JOB" \
        --genomeInfo "$GENOMEINFO" \
        --ignoreGenomeQuality \
        -pa 0.9 -sa 0.95 -nc 0.30 \
        -cm larger \
        --S_algorithm skani \
        --multiround_primary_clustering \
        --primary_chunksize 2000 \
        -g "$list" \
        > "$STAGE1_DIR/${subgroup}.log" 2>&1
    then
        echo "[DONE]  $subgroup ($(ls "$out/dereplicated_genomes/" | wc -l) reps) @ $(date +%H:%M:%S)"
    else
        echo "[FAIL]  $subgroup — see $STAGE1_DIR/${subgroup}.log"
        return 1
    fi
}
export -f run_one_subgroup
export STAGE1_DIR GENOMEINFO THREADS_PER_JOB

sort -t$'\t' -k2 "$STAGE1_DIR/subgroup_index.txt" | \
while IFS=$'\t' read -r sg list; do
    printf "%d\t%s\t%s\n" "$(wc -l < "$list")" "$sg" "$list"
done | sort -rn | cut -f2,3 | \
parallel --colsep '\t' --jobs "$CONCURRENT_JOBS" --line-buffer \
    run_one_subgroup {1} {2}

echo ""
echo "========================================================"
echo "STAGE 1 (parallel) finished: $(date)"
n_done=$(grep -l "dereplicated_genomes" "$STAGE1_DIR"/*_drep/data_tables/Cdb.csv 2>/dev/null | wc -l || echo "?")
echo "Check for any failures:"
echo "  grep -L . $STAGE1_DIR/*.log   # empty logs = possible issues"
echo "  ls $STAGE1_DIR/*_drep/dereplicated_genomes/ 2>/dev/null | wc -l"
echo ""
echo "If any subgroup FAILED, just rerun this script — completed subgroups"
echo "are skipped automatically, so only failed ones retry."
echo ""
echo "When all subgroups are done, proceed to Stage 2: 3p_dereplication.sh species_stage2"
echo "========================================================"
)

step_species_stage2() (
set -euo pipefail

WORKDIR="${URGG_PROJECT}/drep_analysis"
THREADS=48

STAGE1_DIR="$WORKDIR/stage1_per_host"
GENOMEINFO="$WORKDIR/genomeInfo.csv"
POOLED_LIST="$WORKDIR/stage2_pooled_reps.txt"
OUTPUT="$WORKDIR/drep_95_species_final"
LOGFILE="$WORKDIR/drep_95_stage2_run.log"

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
)

step_strains() (
set -euo pipefail

WORKDIR="${URGG_PROJECT}/drep_analysis"
THREADS=24

METADATA="$WORKDIR/metadata.tsv"
GENOMEINFO="$WORKDIR/genomeInfo.csv"
GENOME_LIST="$WORKDIR/genome_list_all_HQ.txt"
OUTPUT="$WORKDIR/drep_99_strain"
LOGFILE="$WORKDIR/drep_99_strain_run.log"

if [[ ! -f "$GENOME_LIST" ]]; then
    echo ">> Building all-site HQ genome list from metadata.tsv..."
    awk -F'\t' 'NR>1 && $6=="HQ" {print $1}' "$METADATA" \
        | sed "s|^|$WORKDIR/genomes/|" > "$GENOME_LIST"
    n=$(wc -l < "$GENOME_LIST")
    echo "   Written $n HQ genomes (all sites) to $GENOME_LIST"
    echo ""

    # show site distribution of what's going in
    echo "   Site distribution of HQ genomes:"
    awk -F'\t' 'NR>1 && $6=="HQ" {print $5}' "$METADATA" \
        | sort | uniq -c | sort -rn \
        | awk '{printf "     %-25s %d genomes\n", $2, $1}'
    echo ""
fi

echo "========================================================"
echo "Tier B: Strain-level dereplication (99% ANI)"
echo "Input genomes: $(wc -l < "$GENOME_LIST") (all sites, HQ only)"
echo "Threads: $THREADS"
echo "Output: $OUTPUT"
echo "Started: $(date)"
echo "========================================================"

# validate inputs
for f in "$GENOME_LIST" "$GENOMEINFO"; do
    [[ -f "$f" ]] || { echo "ERROR: required file not found: $f" >&2; exit 1; }
done

/usr/bin/time -v dRep dereplicate "$OUTPUT" \
    -p $THREADS \
    --genomeInfo "$GENOMEINFO" \
    --ignoreGenomeQuality \
    -pa 0.95 -sa 0.99 -nc 0.30 \
    -cm larger \
    --S_algorithm skani \
    --multiround_primary_clustering \
    --primary_chunksize 2000 \
    -g "$GENOME_LIST" \
    2>&1 | tee "$LOGFILE"

echo ""
echo "========================================================"
echo "Tier B complete: $(date)"
echo "Strain representatives: $OUTPUT/dereplicated_genomes/"
echo "Rep count: $(ls $OUTPUT/dereplicated_genomes/ | wc -l)"
echo ""
echo "Key output tables (preserve these):"
echo "  $OUTPUT/data_tables/Cdb.csv   cluster assignments (strain)"
echo "  $OUTPUT/data_tables/Wdb.csv   winning representatives"
echo "  $OUTPUT/data_tables/Sdb.csv   genome scores"
echo "  $OUTPUT/data_tables/Ndb.csv   pairwise ANI cache"
echo ""
echo "Next step: build SGVFinder2 reference database from these representatives:"
echo "  mkdir -p $WORKDIR/sv_reference_db"
echo "  cp $OUTPUT/dereplicated_genomes/*.fa $WORKDIR/sv_reference_db/"
echo "========================================================"
)

for s in ${@:-species_stage1 species_stage2 strains}; do
  echo "==> 3p_dereplication.sh: $s"
  "step_$s"
done
