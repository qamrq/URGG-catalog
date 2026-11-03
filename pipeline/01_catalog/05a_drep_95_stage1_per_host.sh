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
#  SCRIPT  : 05a_drep_95_stage1_per_host.sh  -  Species dereplication stage 1: per host and BioProject (dRep, 95% ANI)
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/genomes/, metadata.tsv, genomeInfo.csv
#  OUTPUT  : drep_analysis/stage1_per_host/<subgroup>_drep/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

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
    echo "  conda install -c conda-forge parallel" >&2
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
echo "When all subgroups are done, proceed to Stage 2: 05b_drep_95_stage2_final.sh"
echo "========================================================"
