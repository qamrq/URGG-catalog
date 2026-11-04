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
#  SCRIPT  : 7a_read_mapping_abundance.sh  -  Sample manifest and genome abundance per sample (CoverM)
#  AUTHOR  : QRQ
#  INPUT   : $URGG_HOCORT clean read pairs
#            species genomes
#  OUTPUT  : phase_4_SV_analysis/working/samples_manifest.tsv
#            annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv
#  USAGE   : bash 7a_read_mapping_abundance.sh [manifest coverm]   (no argument: all steps in this order)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_HOCORT="${URGG_HOCORT:-/path/to/hocort_clean_reads}"

set -euo pipefail

step_manifest() (
set -euo pipefail

HOCORT_ROOT="${URGG_HOCORT}"
WORK_DIR="${URGG_PROJECT}/phase_4_SV_analysis/working"
MANIFEST="${WORK_DIR}/samples_manifest.tsv"
TEST_MANIFEST="${WORK_DIR}/samples_manifest_test.tsv"
UNMAPPED_LOG="${WORK_DIR}/manifest_unmapped_species.txt"
UNPAIRED_LOG="${WORK_DIR}/manifest_unpaired_samples.txt"
TEST_N=5

SKIP_PROJECTS="PRJCA024157"

mkdir -p "${WORK_DIR}"

declare -A SPECIES_MAP=(
  # --- flat single-species projects ---
  ["PRJEB10338"]="Bos_taurus"
  ["PRJNA214227"]="Ovis_aries"
  ["PRJNA526070"]="Bos_taurus"
  ["PRJEB14312"]="Ovis_aries"
  ["PRJNA448333"]="Bos_taurus"
  ["PRJNA603266"]="Camelus_dromedarius"
  ["PRJNA60251"]="Bos_taurus"
  ["PRJNA512570"]="Ovis_aries"
  ["PRJNA485488"]="Bos_taurus"
  ["PRJNA362629"]="Bos_taurus"
  ["PRJNA419351"]="Bos_taurus"
  ["PRJNA1069655"]="Ovis_aries"
  ["PRJNA1034944"]="Giraffa_camelopardalis"
  ["PRJNA260105"]="Bos_taurus"
  ["PRJNA298982"]="Alces_americanus"
  ["PRJEB8939"]="Bos_taurus"
  ["PRJEB12797"]="Alces_alces"
  ["PRJNA723218"]="Bos_taurus"
  ["PRJNA301235"]="Alces_alces"

  # --- nested multi-species: PROJECT/ANIMAL -> species ---
  ["PRJCA008889/Goat"]="Capra_hircus"
  ["PRJCA008889/Sheep"]="Ovis_aries"
  ["PRJCA008889/Antelope"]="Pantholops_hodgsonii"
  ["PRJCA008889/Gazelle"]="Procapra_picticaudata"

  ["PRJEB38078/Camelus_dromedarius"]="Camelus_dromedarius"
  ["PRJEB38078/Camelus_bact"]="Camelus_bactrianus"
  ["PRJEB38078/Gaur"]="Bos_frontalis"
  ["PRJEB38078/Giraffe"]="Giraffa_camelopardalis"

  ["PRJNA486341/Camel_drom"]="Camelus_dromedarius"
  ["PRJNA486341/Goat"]="Capra_hircus"

  ["PRJNA983076/Camelus_bactrianus"]="Camelus_bactrianus"
  ["PRJNA983076/Giraffa_camelopardalis"]="Giraffa_camelopardalis"

  ["PRJNA972320/Goat"]="Capra_hircus"
  ["PRJNA972320/Sheep"]="Ovis_aries"

  ["PRJNA657455/Dairy_cattle"]="Bos_taurus"
  ["PRJNA657455/Yak"]="Bos_grunniens"
  ["PRJNA657455/Sheep"]="Ovis_aries"
  ["PRJNA657455/Goat"]="Capra_hircus"
  ["PRJNA657455/Water_deer"]="Hydropotes_inermis"
  ["PRJNA657455/Water_buffalo"]="Bubalus_bubalis"
  ["PRJNA657455/Roe_deer"]="Capreolus_pygargus"
)

canon_proj() { echo "${1%_DONE}"; }

is_aux_dir() {
  case "$(basename "$1")" in
    *log*|*fastp*|*hocort_log*|*fasp_log*|log_files_*|hocort_done) return 0 ;;
    *) return 1 ;;
  esac
}

echo -e "sample_id\tfq1\tfq2\tbioproject\tspecies\tanimal_label" > "${MANIFEST}"
: > "${UNMAPPED_LOG}"
: > "${UNPAIRED_LOG}"

total=0; paired=0; unpaired=0; unmapped=0

process_r1() {
    local r1="$1"
    local dir; dir="$(dirname "$r1")"
    local base; base="$(basename "$r1")"

    is_aux_dir "$dir" && return 0

    local sample_id="" r2=""
    if [[ "$base" == *_1_clean.fastq.gz ]]; then
        sample_id="${base%_1_clean.fastq.gz}"
        r2="${dir}/${sample_id}_2_clean.fastq.gz"
    elif [[ "$base" == *_clean_R1.fastq.gz ]]; then
        sample_id="${base%_clean_R1.fastq.gz}"
        r2="${dir}/${sample_id}_clean_R2.fastq.gz"
    else
        return 0
    fi

    total=$((total+1))

    local rel="${r1#${HOCORT_ROOT}/}"
    local proj_folder="${rel%%/*}"
    local proj; proj="$(canon_proj "$proj_folder")"

    case " ${SKIP_PROJECTS} " in *" ${proj} "*|*" ${proj_folder} "*) return 0 ;; esac

    local after="${rel#${proj_folder}/}"
    local animal=""
    if [[ "$after" == */* ]]; then
        local cand="${after%%/*}"
        if ! is_aux_dir "${HOCORT_ROOT}/${proj_folder}/${cand}"; then
            animal="$cand"
        fi
    fi

    local species=""
    if [[ -n "$animal" && -n "${SPECIES_MAP[${proj}/${animal}]:-}" ]]; then
        species="${SPECIES_MAP[${proj}/${animal}]}"
    elif [[ -n "${SPECIES_MAP[${proj}]:-}" ]]; then
        species="${SPECIES_MAP[${proj}]}"
    else
        species="UNKNOWN"
        echo -e "${proj}\t${animal}\t${r1}" >> "${UNMAPPED_LOG}"
        unmapped=$((unmapped+1))
    fi

    if [[ ! -f "$r2" ]]; then
        echo -e "${sample_id}\t${r1}\t(missing R2)" >> "${UNPAIRED_LOG}"
        unpaired=$((unpaired+1))
        return 0
    fi

    local uid="${proj}__${animal:+${animal}__}${sample_id}"

    echo -e "${uid}\t${r1}\t${r2}\t${proj}\t${species}\t${animal:-NA}" >> "${MANIFEST}"
    paired=$((paired+1))
}

echo "Scanning ${HOCORT_ROOT} ..."
while IFS= read -r -d '' r1; do
    process_r1 "$r1"
done < <(find "${HOCORT_ROOT}" -type f \
            \( -name '*_1_clean.fastq.gz' -o -name '*_clean_R1.fastq.gz' \) -print0)

head -1 "${MANIFEST}" > "${TEST_MANIFEST}"
tail -n +2 "${MANIFEST}" | head -"${TEST_N}" >> "${TEST_MANIFEST}"

echo "----------------------------------------------"
echo "Total R1 files scanned : ${total}"
echo "Paired samples written : ${paired}"
echo "Unpaired (missing R2)  : ${unpaired}  -> ${UNPAIRED_LOG}"
echo "Unmapped species       : ${unmapped}  -> ${UNMAPPED_LOG}"
echo ""
echo "Full manifest : ${MANIFEST}"
echo "Test manifest : ${TEST_MANIFEST} (${TEST_N} samples)"
echo ""
echo "Per-species sample counts:"
tail -n +2 "${MANIFEST}" | cut -f5 | sort | uniq -c | sort -rn
echo ""
echo "Per-project sample counts:"
tail -n +2 "${MANIFEST}" | cut -f4 | sort | uniq -c | sort -rn
echo ""
if [[ "${unmapped}" -gt 0 ]]; then
    echo "WARNING: ${unmapped} reads had UNKNOWN species - review ${UNMAPPED_LOG}"
    echo "         and add the missing project/animal keys to SPECIES_MAP."
fi
echo ""
echo "Preview:"
head -6 "${MANIFEST}" | column -t -s$'\t'
)

step_coverm() (
set -euo pipefail

main() {

MANIFEST=${URGG_PROJECT}/phase_4_SV_analysis/working/samples_manifest.tsv
REF_DIR=${URGG_PROJECT}/annotation/results_phase3/s35_diversity/reference

export TMPDIR=${URGG_PROJECT}/annotation/results_phase3/s35_diversity/tmp
mkdir -p "${TMPDIR}"
COMBINED_REF="${REF_DIR}/combined_reference.fasta"
GENOME_DEF="${REF_DIR}/genome_definition.tsv"

OUT=${URGG_PROJECT}/annotation/results_phase3/s35_diversity/abundance
WORK="${OUT}/work"
BATCH_DIR="${OUT}/batches"
LOGS=${URGG_PROJECT}/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${BATCH_DIR}" "${LOGS}"

exec > >(tee -a "${LOGS}/35b_coverm_abundance_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-64}
N_PARALLEL=${N_PARALLEL:-3}
MIN_IDENTITY=95
MIN_ALIGNED_PCT=75

TEST_BIOPROJECT="${1:-}"

[[ -s "${MANIFEST}" ]] || { log "ERROR: manifest not found: ${MANIFEST}"; exit 1; }
command -v coverm >/dev/null 2>&1 || { log "ERROR: coverm not on PATH."; exit 1; }

[[ -s "${COMBINED_REF}" ]] || { log "ERROR: ${COMBINED_REF} not found. Run 35b0_build_combined_reference.sh first."; exit 1; }
[[ -s "${GENOME_DEF}" ]] || { log "ERROR: ${GENOME_DEF} not found. Run 35b0_build_combined_reference.sh first."; exit 1; }
if ! find "${REF_DIR}" -maxdepth 1 -name "combined_reference.fasta.r*" -print -quit | grep -q .; then
    log "ERROR: no strobealign index (combined_reference.fasta.r*.sti) found in ${REF_DIR}."
    log "       Run 35b0_build_combined_reference.sh's index-build step first."
    exit 1
fi

log "CoverM: $(coverm --version 2>&1 | head -1)"

log "Reading manifest..."
N_TOTAL=$(($(wc -l < "${MANIFEST}") - 1))
log "Total samples in manifest: ${N_TOTAL}"

tail -n +2 "${MANIFEST}" | cut -f4 | sort -u > "${WORK}/bioprojects.txt"

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
        rm -f "${out_tsv}.raw"
        return 1
    fi

    if [[ ! -s "${out_tsv}.raw" ]]; then
        log "ERROR [${bp}]: coverm reported success but produced an empty/missing output file."
        rm -f "${out_tsv}.raw"
        return 1
    fi

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
        lookup[base] = row["sample_id"]  # fallback: exact basename match

with open(raw_f, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    rows = list(rd)

new_header = [header[0]]
unmatched = []
for col in header[1:]:
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

export -f run_batch log
export MANIFEST COMBINED_REF GENOME_DEF MIN_IDENTITY MIN_ALIGNED_PCT THREADS BATCH_DIR LOGS

FAILED_LOG="${WORK}/failed_bioprojects.txt"
: > "${FAILED_LOG}"

xargs -a "${WORK}/bioprojects.txt" -I{} -P "${N_PARALLEL}" bash -c 'run_batch "$@"' _ {} || true

for bp in $(cat "${WORK}/bioprojects.txt"); do
    if [[ ! -s "${BATCH_DIR}/${bp}.relabund.tsv" ]]; then
        echo "${bp}" >> "${FAILED_LOG}"
    fi
done
N_FAILED=$(wc -l < "${FAILED_LOG}")
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "WARNING: ${N_FAILED} BioProject(s) failed or did not complete - see ${FAILED_LOG}"
    log "         Re-run this script to retry them (already-completed batches will be skipped)."
    cat "${FAILED_LOG}" | sed 's/^/  FAILED: /'
else
    log "All BioProjects completed successfully."
fi

log "All BioProject batches finished (or skipped if already done)."

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
    print("WARNING: no 'unmapped' row found - mapping-rate summary not produced.")
PY

log "-------------------------------------------------------"
log "Done. Outputs in ${OUT}"
log "  genome_x_sample_relative_abundance.tsv   full abundance matrix"
log "  mapping_rate_per_sample.tsv              per-sample % reads mapped to catalogue"
log "-------------------------------------------------------"

}

main "$@"
)

for s in ${@:-manifest coverm}; do
  echo "==> 7a_read_mapping_abundance.sh: $s"
  "step_$s"
done
