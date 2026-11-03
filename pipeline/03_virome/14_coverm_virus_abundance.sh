#!/usr/bin/env bash
# ============================================================
# 36_coverm_virus_abundance.sh
#
# Builds the vOTU x sample relative-abundance matrix, needed for the
# ruminal virus co-occurrence network (mirrors the design in published
# rumen-virome co-occurrence figures: vOTU-vOTU correlation network,
# node size = differential abundance across host species/breed).
#
# ADAPTED FROM 35b_coverm_abundance.sh (bacterial/archaeal catalogue
# abundance). Reuses the SAME samples_manifest.tsv, BioProject batching,
# resumability, and column-renaming logic. TWO deliberate differences
# from 35b, both noted inline below:
#   1. CoverM CONTIG mode (--reference, single multi-FASTA) instead of
#      GENOME mode (--genome-fasta-directory). Each vOTU representative
#      in votu_representatives.fa is already a single sequence, not a
#      multi-contig assembly requiring genome-level aggregation, so
#      contig mode is the correct unit here.
#   2. trimmed_mean method + --min-covered-fraction 0.7, matching
#      established viral-vOTU abundance profiling convention (guards
#      against a vOTU getting inflated abundance from a short, spurious
#      partial-coverage read match) rather than plain relative_abundance
#      as used for the bacterial/archaeal catalogue.
#
# SCOPE NOTE (carried over from 35b, unchanged): read-mapping-based
# abundance can only cover host species with QC'd raw reads available.
# Confirm the same 14-of-17 host-species scope applies here (three
# species — Bos indicus, Cervus elaphus, Rangifer tarandus — have no
# QC'd reads and will be absent from every output this script produces).
# State this explicitly in Methods/Results, consistent with 35b's
# framing.
# ============================================================
set -euo pipefail

MANIFEST=/disk252/qamar/metagenome_analysis/phase_4_SV_analysis/working/samples_manifest.tsv
VOTU_FA=/disk252/qamar/metagenome_analysis/phase_5_virome/07_votu_clustering/votu_representatives.fa

OUT=/disk252/qamar/metagenome_analysis/annotation/results_phase3/s36_virus_abundance
WORK="${OUT}/work"
BATCH_DIR="${OUT}/batches"
LOGS=/disk252/qamar/metagenome_analysis/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${BATCH_DIR}" "${LOGS}"

exec > >(tee -a "${LOGS}/36_coverm_virus_abundance_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-48}          # per-batch thread count (CoverM's own -t)
N_PARALLEL=${N_PARALLEL:-3}     # how many BioProject batches to run AT ONCE
MIN_IDENTITY=95
MIN_ALIGNED_PCT=75
MIN_COVERED_FRACTION=0.7        # NEW vs. 35b: guards against short spurious-coverage hits

# ------------------------------------------------------------------------
# PARALLELISM / MEMORY NOTE (carried over from 35b): CoverM builds an
# index of the full reference on EVERY invocation. Running N_PARALLEL
# batches simultaneously multiplies memory use by N_PARALLEL. The vOTU
# reference (8,980 sequences) is far smaller than the 25,451-genome
# bacterial reference, so peak memory per batch should be substantially
# lower here than in 35b -- but verify once before trusting the default:
#   /usr/bin/time -v bash 36_coverm_virus_abundance.sh <small_bioproject>
# and check "Maximum resident set size" before raising N_PARALLEL.
# ------------------------------------------------------------------------

# Optional test mode: pass one BioProject ID as the first argument to
# restrict this run to that project only, e.g.:
#   bash 36_coverm_virus_abundance.sh PRJNA972320
# Omit the argument to process every BioProject in the manifest.
TEST_BIOPROJECT="${1:-}"

[[ -s "${MANIFEST}" ]] || { log "ERROR: manifest not found: ${MANIFEST}"; exit 1; }
[[ -s "${VOTU_FA}" ]] || { log "ERROR: vOTU FASTA not found: ${VOTU_FA}"; exit 1; }
command -v coverm >/dev/null 2>&1 || { log "ERROR: coverm not on PATH. conda activate the env that provides it."; exit 1; }

log "CoverM: $(coverm --version 2>&1 | head -1)"

# ============================================================
# STEP 1 — batch samples by BioProject (already resolved in the manifest)
# ============================================================
log "Reading manifest..."
N_TOTAL=$(($(wc -l < "${MANIFEST}") - 1))
log "Total samples in manifest: ${N_TOTAL}"

tail -n +2 "${MANIFEST}" | cut -f4 | sort -u > "${WORK}/bioprojects.txt"

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

    log "=== ${bp}: running CoverM (contig mode, vOTU reference) on $((${#coupled_args[@]} / 2)) samples ==="
    # DIFFERENCE vs. 35b: --reference (single multi-FASTA, contig mode)
    # instead of --genome-fasta-directory (genome mode). DIFFERENCE vs.
    # 35b: trimmed_mean + --min-covered-fraction instead of plain
    # relative_abundance.
    coverm contig \
        --coupled "${coupled_args[@]}" \
        --reference "${VOTU_FA}" \
        --min-read-percent-identity "${MIN_IDENTITY}" \
        --min-read-aligned-percent "${MIN_ALIGNED_PCT}" \
        --min-covered-fraction "${MIN_COVERED_FRACTION}" \
        -m trimmed_mean \
        -t "${THREADS}" \
        -o "${out_tsv}.raw"

    # Same column-renaming logic as 35b, unchanged: CoverM names output
    # columns from the FASTQ filename, not the manifest's sample_id.
    #
    # VERIFY BEFORE TRUSTING DOWNSTREAM RESULTS: after the first batch
    # completes, inspect ${out_tsv}.raw's header row against this
    # renaming logic -- if columns are left unrenamed (a "COULD NOT
    # MATCH" warning below), the basename-stripping pattern needs
    # adjusting to match CoverM's actual naming convention in this
    # version. If 35b already validated this pattern against the same
    # manifest, it should carry over unchanged here.
    python3 - "${out_tsv}.raw" "${MANIFEST}" "${bp}" "${out_tsv}.tmp" <<'PY'
import sys, csv, os, re

raw_f, manifest_f, bp, out_f = sys.argv[1:5]

lookup = {}
with open(manifest_f, newline="") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for row in rd:
        if row["bioproject"] != bp:
            continue
        base = os.path.basename(row["fq1"])
        stripped = re.sub(r"(_1_clean|_clean_R1|_1)\.fastq\.gz$", "", base)
        lookup[stripped] = row["sample_id"]
        lookup[base] = row["sample_id"]

with open(raw_f, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    rows = list(rd)

new_header = [header[0]]
unmatched = []
for col in header[1:]:
    if col in lookup:
        new_header.append(lookup[col])
    else:
        stripped_col = re.sub(r"(_1_clean|_clean_R1|_1)$", "", col)
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

export -f run_batch log
export MANIFEST VOTU_FA MIN_IDENTITY MIN_ALIGNED_PCT MIN_COVERED_FRACTION THREADS BATCH_DIR

xargs -a "${WORK}/bioprojects.txt" -I{} -P "${N_PARALLEL}" bash -c 'run_batch "$@"' _ {}

log "All BioProject batches finished (or skipped if already done)."

# ============================================================
# STEP 3 — merge all batch outputs into one vOTU x sample matrix
# ============================================================
log "Merging batch outputs..."
python3 - "${BATCH_DIR}" "${OUT}" <<'PY'
import os, sys, csv
from collections import defaultdict

batch_dir, out_dir = sys.argv[1], sys.argv[2]
files = sorted(f for f in os.listdir(batch_dir) if f.endswith(".relabund.tsv"))

all_samples = []
votu_data = defaultdict(dict)

for fn in files:
    path = os.path.join(batch_dir, fn)
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        samples = header[1:]
        all_samples.extend(samples)
        for row in rd:
            votu = row[0]
            for s, v in zip(samples, row[1:]):
                votu_data[votu][s] = v

out_matrix = os.path.join(out_dir, "votu_x_sample_trimmed_mean.tsv")
with open(out_matrix, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["vOTU"] + all_samples)
    for votu, vals in votu_data.items():
        w.writerow([votu] + [vals.get(s, "0") for s in all_samples])

print(f"Merged matrix: {len(votu_data)} vOTUs x {len(all_samples)} samples")
print(f"Written: {out_matrix}")

# NOTE: unlike 35b's genome-mode CoverM run, contig-mode output does NOT
# produce an 'unmapped' pseudo-reference row, so there is no equivalent
# mapping-rate summary step here. If a virus-specific mapping-rate
# statistic is needed later, it would require a separate CoverM
# invocation configured for that purpose.
PY

log "-------------------------------------------------------"
log "Done. Outputs in ${OUT}"
log "  votu_x_sample_trimmed_mean.tsv   full vOTU abundance matrix"
log "-------------------------------------------------------"
