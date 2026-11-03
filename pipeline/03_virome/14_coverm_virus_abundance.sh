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
#  SCRIPT  : 14_coverm_virus_abundance.sh  -  vOTU abundance per sample (CoverM)
#  AUTHOR  : QRQ
#  INPUT   : phase_4_SV_analysis/working/samples_manifest.tsv
#            phase_5_virome/07_votu_clustering/votu_representatives.fa
#  OUTPUT  : annotation/results_phase3/s36_virus_abundance/votu_x_sample_trimmed_mean.tsv
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

MANIFEST=${URGG_PROJECT}/phase_4_SV_analysis/working/samples_manifest.tsv
VOTU_FA=${URGG_PROJECT}/phase_5_virome/07_votu_clustering/votu_representatives.fa

OUT=${URGG_PROJECT}/annotation/results_phase3/s36_virus_abundance
WORK="${OUT}/work"
BATCH_DIR="${OUT}/batches"
LOGS=${URGG_PROJECT}/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${BATCH_DIR}" "${LOGS}"

exec > >(tee -a "${LOGS}/36_coverm_virus_abundance_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-48}
N_PARALLEL=${N_PARALLEL:-3}
MIN_IDENTITY=95
MIN_ALIGNED_PCT=75
MIN_COVERED_FRACTION=0.7

TEST_BIOPROJECT="${1:-}"

[[ -s "${MANIFEST}" ]] || { log "ERROR: manifest not found: ${MANIFEST}"; exit 1; }
[[ -s "${VOTU_FA}" ]] || { log "ERROR: vOTU FASTA not found: ${VOTU_FA}"; exit 1; }
command -v coverm >/dev/null 2>&1 || { log "ERROR: coverm not on PATH. conda activate the env that provides it."; exit 1; }

log "CoverM: $(coverm --version 2>&1 | head -1)"

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
    coverm contig \
        --coupled "${coupled_args[@]}" \
        --reference "${VOTU_FA}" \
        --min-read-percent-identity "${MIN_IDENTITY}" \
        --min-read-aligned-percent "${MIN_ALIGNED_PCT}" \
        --min-covered-fraction "${MIN_COVERED_FRACTION}" \
        -m trimmed_mean \
        -t "${THREADS}" \
        -o "${out_tsv}.raw"

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
