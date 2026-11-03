#!/usr/bin/env bash
# ============================================================
# 35b_coverm_abundance.sh
#
# RESULTS SECTION 3.1/3.2 — builds the genome x sample relative-abundance
# matrix required for mapping-rate validation and alpha/beta diversity
# (matching RGMGC Fig 1e/1g and the Tibetan Plateau paper's Fig 2B/3A).
# This is the SAME underlying data cluster 3932_1's contamination
# question needed — once this matrix exists, that question becomes a
# single row lookup rather than a separate targeted mapping run (33c is
# superseded by this, once complete).
#
# REUSES the existing samples_manifest.tsv built for SV analysis
# (build_manifest.sh) rather than re-deriving read locations — that
# manifest already resolved every naming/layout inconsistency across
# hocort_results/.
#
# METHOD: CoverM genome mode, matching standard MAG-catalogue abundance
# profiling thresholds (--min-read-percent-identity 95
# --min-read-aligned-percent 75), against the dereplicated 95% ANI
# species catalogue. CoverM's relative_abundance output includes an
# "unmapped" pseudo-genome row, which directly gives the per-sample
# mapping-rate statistic (100 - unmapped%) without a separate step.
#
# SCOPE NOTE: batched by BioProject (using the manifest's existing
# bioproject column) rather than one all-samples invocation, so partial
# progress survives an interruption and batches can run in parallel.
#
# ROBUSTNESS FIX: the ENTIRE script body is wrapped in a single main()
# function, called once at the very bottom. Bash parses a function's
# whole body as one unit the moment it's defined -- so wrapping the
# whole script this way makes it immune to being corrupted if the file
# on disk is edited/overwritten WHILE this script is still running
# (bash normally reads a running script incrementally by file position,
# so an in-place edit mid-run scrambles everything read afterward --
# this is what caused a "line 284: syntax error near unexpected token
# '('" crash on a previous run, on an otherwise-inert Python heredoc
# line). Only the final one-line `main "$@"` call remains exposed to a
# lazy re-read, which is trivial.
# ============================================================
set -euo pipefail

main() {

# SCOPE NOTE (confirmed from samples_manifest.tsv): read-mapping-based
# analyses (this script and 35c) can only cover host species with at
# least one raw-FASTQ source processed through this project's own
# QC pipeline. Confirmed via the manifest's species column: three of the
# catalogue's 17 host species have NO QC'd reads available at all —
# Bos indicus, Cervus elaphus, and Rangifer tarandus — and are absent
# from every output this script produces. Bubalus bubalis DOES have
# coverage (53 samples, PRJNA657455/Water_buffalo) despite its catalogue
# representative genomes originating from a different, MAG-only
# BioProject (PRJNA656389) — read-mapping only requires QC'd reads from
# the species, not that those specific reads were ever assembled/binned.
# State the 14-of-17 scope explicitly in Methods/Results.

MANIFEST=/disk252/qamar/metagenome_analysis/phase_4_SV_analysis/working/samples_manifest.tsv
REF_DIR=/disk252/qamar/metagenome_analysis/annotation/results_phase3/s35_diversity/reference

# Redirect ALL temp files (CoverM's FIFO pipes, samtools sort's
# intermediate chunks) away from the small root-partition /tmp (439GB
# total, only ~95GB free when 3 large parallel batches run
# simultaneously) and onto /disk252, which has ~850GB free. A batch
# failure (PRJCA008889, "samtools sort: truncated file") occurred with
# the default /tmp — plausibly disk-space exhaustion under 3-way
# parallel load, though also check `dmesg` for OOM as an alternative
# cause before assuming this alone fixes it. CoverM's "coverm_fifo.*"
# naming is the standard signature of Rust's std::env::temp_dir(), which
# respects $TMPDIR on Unix — no CoverM-specific flag needed.
export TMPDIR=/disk252/qamar/metagenome_analysis/annotation/results_phase3/s35_diversity/tmp
mkdir -p "${TMPDIR}"
COMBINED_REF="${REF_DIR}/combined_reference.fasta"
GENOME_DEF="${REF_DIR}/genome_definition.tsv"

OUT=/disk252/qamar/metagenome_analysis/annotation/results_phase3/s35_diversity/abundance
WORK="${OUT}/work"
BATCH_DIR="${OUT}/batches"
LOGS=/disk252/qamar/metagenome_analysis/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${BATCH_DIR}" "${LOGS}"

exec > >(tee -a "${LOGS}/35b_coverm_abundance_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-64}          # per-batch thread count (CoverM's own -t)
N_PARALLEL=${N_PARALLEL:-3}     # how many BioProject batches to run AT ONCE
MIN_IDENTITY=95
MIN_ALIGNED_PCT=75

# ------------------------------------------------------------------------
# PARALLELISM / MEMORY NOTE: CoverM loads/builds an index of the full
# genome reference (25,451 genomes) on EVERY invocation. Running
# N_PARALLEL batches simultaneously multiplies memory use by N_PARALLEL,
# not just CPU use. Before raising N_PARALLEL, measure one batch's actual
# peak memory (e.g. `/usr/bin/time -v bash 35b_coverm_abundance.sh
# <small_bioproject>` and check "Maximum resident set size"), then set
# N_PARALLEL conservatively below (free_RAM_GB / peak_RSS_GB). The
# defaults above (64 threads/batch x 3 parallel = up to 192 of 256 cores)
# are a starting point, not a verified-safe setting for this specific
# reference size — tune after the memory check.
# ------------------------------------------------------------------------

# Optional test mode: pass one BioProject ID as the first argument to
# restrict this run to that project only, e.g.:
#   bash 35b_coverm_abundance.sh PRJNA972320
# Omit the argument to process every BioProject in the manifest.
TEST_BIOPROJECT="${1:-}"

[[ -s "${MANIFEST}" ]] || { log "ERROR: manifest not found: ${MANIFEST}"; exit 1; }
command -v coverm >/dev/null 2>&1 || { log "ERROR: coverm not on PATH. conda activate the env that provides it."; exit 1; }

[[ -s "${COMBINED_REF}" ]] || { log "ERROR: ${COMBINED_REF} not found. Run 35b0_build_combined_reference.sh first."; exit 1; }
[[ -s "${GENOME_DEF}" ]] || { log "ERROR: ${GENOME_DEF} not found. Run 35b0_build_combined_reference.sh first."; exit 1; }
if ! find "${REF_DIR}" -maxdepth 1 -name "combined_reference.fasta.r*" -print -quit | grep -q .; then
    log "ERROR: no strobealign index (combined_reference.fasta.r*.sti) found in ${REF_DIR}."
    log "       Run 35b0_build_combined_reference.sh's index-build step first."
    exit 1
fi

log "CoverM: $(coverm --version 2>&1 | head -1)"

# ============================================================
# STEP 1 — batch samples by BioProject (already resolved in the manifest)
# ============================================================
log "Reading manifest..."
N_TOTAL=$(($(wc -l < "${MANIFEST}") - 1))
log "Total samples in manifest: ${N_TOTAL}"

tail -n +2 "${MANIFEST}" | cut -f4 | sort -u > "${WORK}/bioprojects.txt"

# Optional: temporarily exclude specific BioProjects from this run (e.g.
# one under separate investigation) without touching the resume logic —
# set SKIP_BIOPROJECTS as a space-separated list before invoking this
# script, e.g.:  SKIP_BIOPROJECTS="PRJCA008889" bash 35b_coverm_abundance.sh
if [[ -n "${SKIP_BIOPROJECTS:-}" ]]; then
    for skip_bp in ${SKIP_BIOPROJECTS}; do
        grep -vx "${skip_bp}" "${WORK}/bioprojects.txt" > "${WORK}/bioprojects.txt.tmp" && \
            mv "${WORK}/bioprojects.txt.tmp" "${WORK}/bioprojects.txt"
        log "Excluding ${skip_bp} from this run (SKIP_BIOPROJECTS)."
    done
fi

if [[ -n "${TEST_BIOPROJECT}" ]]; then
    if ! grep -qx "${TEST_BIOPROJECT}" "${WORK}/bioprojects.txt"; then
        log "ERROR: '${TEST_BIOPROJECT}' not found in the manifest's BioProject list."
        log "       Available BioProjects:"
        cat "${WORK}/bioprojects.txt" | sed 's/^/         /'
        exit 1
    fi
    echo "${TEST_BIOPROJECT}" > "${WORK}/bioprojects.txt"
    log "TEST MODE: restricted to single BioProject: ${TEST_BIOPROJECT}"
fi

N_BP=$(wc -l < "${WORK}/bioprojects.txt")
log "BioProjects to process: ${N_BP}"

# ============================================================
# STEP 2 — per-BioProject CoverM run (resumable)
# ============================================================
run_batch() {
    local bp="$1"
    local out_tsv="${BATCH_DIR}/${bp}.relabund.tsv"
    local per_bp_log="${LOGS}/coverm_${bp}_$(date +%Y%m%d_%H%M%S).log"

    if [[ -s "${out_tsv}" ]]; then
        log "SKIP ${bp}: already done."
        return 0
    fi

    local coupled_args=()
    while IFS=$'\t' read -r sample_id fq1 fq2 bioproject species animal; do
        [[ "${bioproject}" == "${bp}" ]] || continue
        [[ -s "${fq1}" && -s "${fq2}" ]] || { log "  WARN: missing reads for ${sample_id}, skipping"; continue; }
        coupled_args+=("${fq1}" "${fq2}")
    done < <(tail -n +2 "${MANIFEST}")

    if [[ "${#coupled_args[@]}" -eq 0 ]]; then
        log "SKIP ${bp}: no valid sample pairs found."
        return 0
    fi

    log "=== ${bp}: running CoverM on $((${#coupled_args[@]} / 2)) samples (details: ${per_bp_log}) ==="
    # CoverM's own stderr is redirected to a PER-BIOPROJECT log file, not
    # the shared global log. When multiple batches run concurrently
    # (N_PARALLEL > 1), their multi-line CoverM output previously
    # interleaved in the single shared log — e.g. PRJCA008889's error
    # block appearing to belong to PRJNA362629's log section purely
    # because both were writing to the same file at the same moment.
    # Only single-line status updates (via `log`) go to the shared log;
    # full CoverM diagnostic output is now always traceable to the
    # correct BioProject via its own dedicated file.
    #
    # --strobealign-params "-r 150" is REQUIRED here: without it, strobealign
    # defaults to looking for an index named "...r100.sti" regardless of the
    # read length it estimates from the actual reads, even though the index
    # built in 35b0 is "...r150.sti" (matching this project's dominant read
    # length). Confirmed via a failed run where strobealign correctly logged
    # "Estimated read length: 150 nt" but then tried to open r100.sti anyway.
    # --min-covered-fraction 0 added following diagnosis: CoverM's default
    # 10% breadth-of-coverage threshold was silently zeroing out entire
    # samples where reads mapped well but were spread thinly across many
    # genomes in this unusually large 25,451-genome catalogue -- confirmed
    # via isolated testing (a sample went from 0.00% to 24.53% mapped with
    # ONLY this flag changed, identity/aligned-percent held at 95%/75%).
    # Applied catalogue-wide for consistency, not just to previously-
    # flagged samples, since any sample's true diversity could have been
    # similarly suppressed by this default.
    if ! coverm genome \
        --coupled "${coupled_args[@]}" \
        --reference "${COMBINED_REF}" \
        --strobealign-use-index \
        --strobealign-params "-r 150" \
        --genome-definition "${GENOME_DEF}" \
        --min-read-percent-identity "${MIN_IDENTITY}" \
        --min-read-aligned-percent "${MIN_ALIGNED_PCT}" \
        --min-covered-fraction 0 \
        -m relative_abundance \
        -t "${THREADS}" \
        -o "${out_tsv}.raw" 2> "${per_bp_log}"
    then
        log "ERROR [${bp}]: coverm genome failed -- full details in ${per_bp_log}"
        log "  ---- last 20 lines of that log ----"
        tail -20 "${per_bp_log}" | sed 's/^/    /'
        log "  ------------------------------------"
        # EXPLICIT check, not reliance on `set -e`: run_batch is invoked
        # inside a separately-spawned `bash -c` subshell via xargs -P, and
        # `set -e` set at the top of THIS script does not propagate into
        # that child invocation. Without this explicit check, a failed
        # coverm call was previously falling through into the renaming
        # step (which then failed confusingly on an empty/missing file)
        # and the batch still got logged as "complete" — silently
        # corrupting the resume logic for any future re-run.
        rm -f "${out_tsv}.raw"
        return 1
    fi

    if [[ ! -s "${out_tsv}.raw" ]]; then
        log "ERROR [${bp}]: coverm reported success but produced an empty/missing output file."
        rm -f "${out_tsv}.raw"
        return 1
    fi

    # CoverM names output columns from the FASTQ filename (typically the
    # shared prefix between R1/R2 after stripping paired-end suffixes),
    # NOT the manifest's sample_id. Since the manifest builds unique IDs
    # specifically as <bioproject>__<animal>__<accession> to avoid
    # collisions, these will not match automatically — remap columns to
    # sample_id using the manifest as the lookup, matching on fq1's
    # basename with common suffixes stripped.
    #
    # VERIFY BEFORE TRUSTING DOWNSTREAM RESULTS: after the first batch
    # completes, inspect ${out_tsv}.raw's header row against this
    # renaming logic — if columns are left unrenamed (a "COULD NOT MATCH"
    # warning below), the basename-stripping pattern needs adjusting to
    # match CoverM's actual naming convention in this version.
    python3 - "${out_tsv}.raw" "${MANIFEST}" "${bp}" "${out_tsv}.tmp" <<'PY'
import sys, csv, os, re

raw_f, manifest_f, bp, out_f = sys.argv[1:5]

# Build fastq-basename -> sample_id lookup for this BioProject only.
lookup = {}
with open(manifest_f, newline="") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for row in rd:
        if row["bioproject"] != bp:
            continue
        base = os.path.basename(row["fq1"])
        # Strip common paired-end / cleaning suffixes to get CoverM's
        # likely column-naming basis.
        stripped = re.sub(r"(_1_clean|_clean_R1|_1)\.fastq\.gz$", "", base)
        lookup[stripped] = row["sample_id"]
        lookup[base] = row["sample_id"]  # fallback: exact basename match

with open(raw_f, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    rows = list(rd)

new_header = [header[0]]
unmatched = []
for col in header[1:]:
    # CoverM's --reference mode produces columns like:
    #   "combined_reference.fasta/SRR27361827_1_clean.fastq.gz Relative Abundance (%)"
    # i.e. a "<reference_basename>/" prefix and a " <Method Name>" suffix
    # wrapped around the actual sample filename. Strip both before
    # attempting to match against the manifest lookup — matching the raw
    # column string directly (as the original version of this script
    # did) never succeeds under --reference mode, since the wrapping is
    # always present.
    col_clean = col.split("/", 1)[1] if "/" in col else col
    m = re.search(r"^(.*\.fastq\.gz)", col_clean)
    fastq_name = m.group(1) if m else col_clean

    if fastq_name in lookup:
        new_header.append(lookup[fastq_name])
    else:
        stripped_col = re.sub(r"(_1_clean|_clean_R1|_1)\.fastq\.gz$", "", fastq_name)
        if stripped_col in lookup:
            new_header.append(lookup[stripped_col])
        else:
            new_header.append(col)
            unmatched.append(col)

with open(out_f, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(new_header)
    w.writerows(rows)

if unmatched:
    print(f"WARNING [{bp}]: {len(unmatched)} column(s) could not be matched to a "
          f"sample_id and were left as-is: {unmatched[:5]}{'...' if len(unmatched)>5 else ''}")
    print("  Inspect these against the manifest and adjust the suffix-stripping "
         "pattern above if this list is nonzero.")
else:
    print(f"[{bp}]: all {len(header)-1} columns matched to sample_id successfully.")
PY
    rm -f "${out_tsv}.raw"
    mv "${out_tsv}.tmp" "${out_tsv}"
    log "${bp} complete -> ${out_tsv}"
}

log "Processing BioProjects (${N_PARALLEL} in parallel, ${THREADS} threads each)..."

# xargs -P runs run_batch for multiple BioProjects concurrently. The
# function and every variable it references must be exported so the
# subshells xargs spawns can see them.
export -f run_batch log
export MANIFEST COMBINED_REF GENOME_DEF MIN_IDENTITY MIN_ALIGNED_PCT THREADS BATCH_DIR LOGS

FAILED_LOG="${WORK}/failed_bioprojects.txt"
: > "${FAILED_LOG}"

# `|| true` is intentional: xargs returns nonzero if ANY single invocation
# fails, and since `set -e` IS active in this outer script (unlike inside
# the run_batch subshells it spawns), letting that propagate would abort
# the entire script on the first failing BioProject — skipping the merge
# step and abandoning every other already-completed batch. Individual
# failures are logged by run_batch itself (which correctly leaves no
# output file behind for a failed batch, so it will be retried on the
# next run); this script continues and reports which ones failed, rather
# than crashing on the first one.
xargs -a "${WORK}/bioprojects.txt" -I{} -P "${N_PARALLEL}" bash -c 'run_batch "$@"' _ {} || true

for bp in $(cat "${WORK}/bioprojects.txt"); do
    if [[ ! -s "${BATCH_DIR}/${bp}.relabund.tsv" ]]; then
        echo "${bp}" >> "${FAILED_LOG}"
    fi
done
N_FAILED=$(wc -l < "${FAILED_LOG}")
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "WARNING: ${N_FAILED} BioProject(s) failed or did not complete — see ${FAILED_LOG}"
    log "         Re-run this script to retry them (already-completed batches will be skipped)."
    cat "${FAILED_LOG}" | sed 's/^/  FAILED: /'
else
    log "All BioProjects completed successfully."
fi

log "All BioProject batches finished (or skipped if already done)."

# ============================================================
# STEP 3 — merge all batch outputs into one genome x sample matrix
# ============================================================
log "Merging batch outputs..."
python3 - "${BATCH_DIR}" "${OUT}" <<'PY'
import os, sys, csv
from collections import defaultdict

batch_dir, out_dir = sys.argv[1], sys.argv[2]
files = sorted(f for f in os.listdir(batch_dir) if f.endswith(".relabund.tsv"))

all_samples = []
genome_data = defaultdict(dict)

for fn in files:
    path = os.path.join(batch_dir, fn)
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        samples = header[1:]
        all_samples.extend(samples)
        for row in rd:
            genome = row[0]
            for s, v in zip(samples, row[1:]):
                genome_data[genome][s] = v

out_matrix = os.path.join(out_dir, "genome_x_sample_relative_abundance.tsv")
with open(out_matrix, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Genome"] + all_samples)
    for genome, vals in genome_data.items():
        w.writerow([genome] + [vals.get(s, "0") for s in all_samples])

print(f"Merged matrix: {len(genome_data)} genomes x {len(all_samples)} samples")
print(f"Written: {out_matrix}")

if "unmapped" in genome_data:
    mr_out = os.path.join(out_dir, "mapping_rate_per_sample.tsv")
    with open(mr_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "pct_unmapped", "pct_mapped_to_catalogue"])
        for s in all_samples:
            unmapped = float(genome_data["unmapped"].get(s, 0) or 0)
            w.writerow([s, f"{unmapped:.3f}", f"{100-unmapped:.3f}"])
    print(f"Mapping-rate summary written: {mr_out}")
else:
    print("WARNING: no 'unmapped' row found — mapping-rate summary not produced.")
PY

log "-------------------------------------------------------"
log "Done. Outputs in ${OUT}"
log "  genome_x_sample_relative_abundance.tsv   full abundance matrix"
log "  mapping_rate_per_sample.tsv              per-sample % reads mapped to catalogue"
log "-------------------------------------------------------"

}  # end main

main "$@"