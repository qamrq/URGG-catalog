#!/usr/bin/env bash
###############################################################################
# batch_checkm2.sh
#
# Runs CheckM2 on multiple MAG folders in a loop.
# Saves results to each folder as:
#   <magdir>/checkm2_output/quality_report.tsv   (raw CheckM2 output)
#   <magdir>/checkm2_summary.tsv                 (clean 5-column summary)
#
# checkm2_summary.tsv columns:
#   Bin_ID | Genome_size | Completeness | Contamination | N50
#
# After this script: manually add sample_site, bioproject, host columns
# to each checkm2_summary.tsv → rename to Best_binset_quality_report_final.tsv
# → add to prep_basalt_for_drep.sh PROJECTS array
###############################################################################
set -euo pipefail

# ── EDIT THESE ────────────────────────────────────────────────────────────────
THREADS=32
EXTENSION="fa"   # extension of your MAG files (fa / fasta / fna)

MAG_DIRS=(
    "/disk252/qamar/basalt_binning/PRJCA008889_DONE/Gazelle/BestBinset_outlier_refined_filtrated/"
)
# ─────────────────────────────────────────────────────────────────────────────

# confirm checkm2 is available
if ! command -v checkm2 &>/dev/null; then
    echo "ERROR: checkm2 not found. Activate the correct conda environment first." >&2
    exit 1
fi

echo "========================================================"
echo "Batch CheckM2 — ${#MAG_DIRS[@]} folders"
echo "Threads: $THREADS | Extension: .$EXTENSION"
echo "Started: $(date)"
echo "========================================================"
echo ""

PASS=0; FAIL=0; SKIP=0

for magdir in "${MAG_DIRS[@]}"; do

    echo "── $(date '+%H:%M:%S') | $(basename $(dirname $magdir))/$(basename $magdir)"

    # check directory exists
    if [[ ! -d "$magdir" ]]; then
        echo "  ERROR: directory not found — skipping" >&2
        FAIL=$((FAIL+1)); continue
    fi

    # count MAG files
    n_mags=$(find "$magdir" -maxdepth 1 -name "*.${EXTENSION}" | wc -l)
    if [[ "$n_mags" -eq 0 ]]; then
        echo "  WARN: no .${EXTENSION} files found — skipping" >&2
        FAIL=$((FAIL+1)); continue
    fi
    echo "  MAGs found: $n_mags"

    checkm2_out="$magdir/checkm2_output"
    summary_tsv="$magdir/checkm2_summary.tsv"

    # skip if already successfully completed
    if [[ -f "$checkm2_out/quality_report.tsv" && -f "$summary_tsv" ]]; then
        n_done=$(wc -l < "$summary_tsv")
        echo "  Already done ($((n_done-1)) MAGs in summary) — skipping"
        echo "  (delete checkm2_output/ to force rerun)"
        SKIP=$((SKIP+1)); continue
    fi

    # run CheckM2
    echo "  Running CheckM2..."
    checkm2 predict \
        --input "$magdir" \
        --output-directory "$checkm2_out" \
        --threads "$THREADS" \
        --extension "$EXTENSION" \
        --force \
        2>&1 | grep -E "INFO|WARNING|ERROR|complete|Finished" || true

    if [[ ! -f "$checkm2_out/quality_report.tsv" ]]; then
        echo "  ERROR: CheckM2 failed — quality_report.tsv not produced" >&2
        FAIL=$((FAIL+1)); continue
    fi

    # format to clean summary TSV
    python3 - << PYEOF
import csv, sys

infile  = "$checkm2_out/quality_report.tsv"
outfile = "$summary_tsv"

with open(infile, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)

# verify required columns
required = {"Name", "Genome_Size", "Completeness", "Contamination", "Contig_N50"}
missing  = required - set(reader.fieldnames)
if missing:
    sys.exit(f"ERROR: CheckM2 report missing columns: {missing}")

with open(outfile, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Bin_ID", "Genome_size", "Completeness", "Contamination", "N50"])
    for row in rows:
        writer.writerow([
            row["Name"],
            row["Genome_Size"],
            row["Completeness"],
            row["Contamination"],
            row["Contig_N50"],
        ])

# quality distribution
hq = mq = fail_q = 0
for row in rows:
    c = float(row["Completeness"])
    x = float(row["Contamination"])
    if   c > 90  and x <= 5: hq     += 1
    elif c >= 70 and x <= 5: mq     += 1
    else:                    fail_q += 1

print(f"  Total: {len(rows)}  |  HQ (>90/<5): {hq}  |  MQ (>=70/<5): {mq}  |  Fail: {fail_q}")
print(f"  Summary saved: {outfile}")
PYEOF

    PASS=$((PASS+1))
    echo ""

done

echo "========================================================"
echo "Finished: $(date)"
echo "Completed: $PASS  |  Skipped (already done): $SKIP  |  Failed: $FAIL"
echo ""
echo "Output files per folder:"
echo "  checkm2_output/quality_report.tsv   (full CheckM2 output)"
echo "  checkm2_summary.tsv                 (Bin_ID|Genome_size|Completeness|Contamination|N50)"
echo ""
echo "Next steps for each folder:"
echo "  1. Add columns: sample_site, bioproject, host to checkm2_summary.tsv"
echo "  2. Rename to:   Best_binset_quality_report_final.tsv"
echo "  3. Add path to: prep_basalt_for_drep.sh PROJECTS array"
echo ""
echo "For multi-site projects (e.g. PRJNA656389), run:"
echo "  python3 add_sample_site_to_quality_report.py \\"
echo "      --quality  <folder>/checkm2_summary.tsv \\"
echo "      --metadata <project_metadata.csv> \\"
echo "      --output   <folder>/Best_binset_quality_report_final.tsv"
