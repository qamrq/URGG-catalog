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
#  SCRIPT  : 8p_virus_detection.sh  -  Viral sequences from strain genomes: VirSorter2, VIBRANT, union and CheckV
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_99_strain/dereplicated_genomes/*.fa (strain genomes)
#  OUTPUT  : confirmed_viral.fa, caller_support.tsv, CheckV quality_summary.tsv
#  USAGE   : bash 8p_virus_detection.sh [batches virsorter2 vibrant checkv]   (no argument: all steps in this order)
# ==========================================================

set -euo pipefail

step_batches() (
source "$(dirname "$0")/config_virome.sh"

log "Building batches of ${BATCH_SIZE} MAGs from ${MAGS_99ANI_DIR}"

mkdir -p "${BATCH_DIR}"

# Refuse to silently append to a half-built set from a previous run
if compgen -G "${BATCH_DIR}/batch_*.fa" > /dev/null; then
    log "ERROR: ${BATCH_DIR} already contains batch_*.fa files."
    log "       Remove them first if you intend to rebuild:  rm -f ${BATCH_DIR}/batch_*.fa ${MAG_INDEX_TSV}"
    exit 1
fi

printf 'mag_index\tmag_id\tbatch_id\tsource_path\n' > "${MAG_INDEX_TSV}"

i=0
batch_n=1
n_in_batch=0
batch_file=""

while IFS= read -r mag_fasta; do
    i=$((i+1))

    # start a new batch every BATCH_SIZE genomes
    if [[ "${n_in_batch}" -eq 0 ]]; then
        batch_id="$(printf 'batch_%04d' "${batch_n}")"
        batch_file="${BATCH_DIR}/${batch_id}.fa"
        : > "${batch_file}"
    fi

    mag_id="$(basename "${mag_fasta}" ".${GENOME_EXT}")"
    mag_index="$(printf 'M%06d' "${i}")"

    awk -v idx="${mag_index}" '
        /^>/ {
            name = substr($0, 2)
            sub(/[ \t].*$/, "", name)          # drop description after first space
            gsub(/[^A-Za-z0-9._-]/, "_", name) # sanitise
            print ">" idx "__" name
            next
        }
        { print }
    ' "${mag_fasta}" >> "${batch_file}"

    printf '%s\t%s\t%s\t%s\n' "${mag_index}" "${mag_id}" "${batch_id}" "${mag_fasta}" >> "${MAG_INDEX_TSV}"

    n_in_batch=$((n_in_batch+1))
    if [[ "${n_in_batch}" -ge "${BATCH_SIZE}" ]]; then
        n_in_batch=0
        batch_n=$((batch_n+1))
    fi
done < <(find "${MAGS_99ANI_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | sort)

n_batches=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | wc -l)

log "Done. ${i} MAGs -> ${n_batches} batches in ${BATCH_DIR}"
log "Mapping table: ${MAG_INDEX_TSV}"

log "--- sanity check: header uniqueness ---"
total_headers=$(cat "${BATCH_DIR}"/batch_*.fa | grep -c '^>')
uniq_headers=$(cat "${BATCH_DIR}"/batch_*.fa | grep '^>' | sort -u | wc -l)
log "  total contigs: ${total_headers}   unique headers: ${uniq_headers}"
if [[ "${total_headers}" -ne "${uniq_headers}" ]]; then
    log "  WARNING: duplicate headers detected -- downstream merging will be ambiguous."
fi
log "--- sanity check: first few headers ---"
head -1 "${BATCH_DIR}/batch_0001.fa"
grep '^>' "${BATCH_DIR}/batch_0001.fa" | head -3
)

step_virsorter2() (
source "$(dirname "$0")/config_virome.sh"

log "Step 1: VirSorter2 on batched MAGs (${BATCH_DIR})"

VS2_CONDA_PREFIX="${DATABASES_DIR}/virsorter2_db/db/conda_envs"

if ! compgen -G "${BATCH_DIR}/batch_*.fa" > /dev/null; then
    log "ERROR: no batch FASTAs found in ${BATCH_DIR}. Run 8p_virus_detection.sh batches first."
    exit 1
fi

run_one_batch () {
    local batch_fa="$1"
    local batch_id
    batch_id="$(basename "${batch_fa}" .fa)"
    local out_dir="${VS2_OUT}/${batch_id}"

    if [[ -f "${out_dir}/final-viral-score.tsv" ]]; then
        echo "[skip] ${batch_id}"
        return 0
    fi

    # --prep-for-dramv writes for-dramv/ which 9p_virus_annotation.sh dramv consumes.
    virsorter run \
        -w "${out_dir}" \
        -i "${batch_fa}" \
        --min-length 1500 \
        --prep-for-dramv \
        -j "${BATCH_THREADS}" \
        all \
        --conda-frontend conda \
        --conda-prefix "${VS2_CONDA_PREFIX}" \
        > "${LOG_DIR}/vs2_${batch_id}.log" 2>&1

    if [[ -f "${out_dir}/final-viral-score.tsv" ]]; then
        echo "[done] ${batch_id}"
    else
        echo "[FAIL] ${batch_id} -- see ${LOG_DIR}/vs2_${batch_id}.log"
    fi
}
export -f run_one_batch
export VS2_OUT LOG_DIR BATCH_THREADS VS2_CONDA_PREFIX

find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | sort | \
    parallel -j "${BATCH_JOBS}" run_one_batch {}

# ---- completion report ----
n_total=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | wc -l)
n_done=$(find "${VS2_OUT}" -mindepth 2 -maxdepth 2 -name "final-viral-score.tsv" | wc -l)
log "Step 1 finished: ${n_done}/${n_total} batches completed."
if [[ "${n_done}" -lt "${n_total}" ]]; then
    log "Incomplete batches (rerun this script to resume):"
    for b in "${BATCH_DIR}"/batch_*.fa; do
        bid="$(basename "${b}" .fa)"
        [[ -f "${VS2_OUT}/${bid}/final-viral-score.tsv" ]] || log "  ${bid}"
    done
fi
log "Next: 8p_virus_detection.sh vibrant, then 8p_virus_detection.sh checkv"
)

step_vibrant() (
source "$(dirname "$0")/config_virome.sh"

log "Step 1 (cont.): VIBRANT on batched MAGs (${BATCH_DIR})"

VIBRANT_DB_DIR="${VIBRANT_DB_DIR:-${DATABASES_DIR}/vibrant/databases}"
VIBRANT_FILES_DIR="${VIBRANT_FILES_DIR:-${DATABASES_DIR}/vibrant/files}"

if [[ ! -d "${VIBRANT_DB_DIR}" ]]; then
    log "ERROR: VIBRANT database dir not found at ${VIBRANT_DB_DIR}"
    exit 1
fi
log "VIBRANT databases: ${VIBRANT_DB_DIR}"

if ! compgen -G "${BATCH_DIR}/batch_*.fa" > /dev/null; then
    log "ERROR: no batch FASTAs in ${BATCH_DIR}. Run 8p_virus_detection.sh batches first."
    exit 1
fi

run_one_batch () {
    local batch_fa="$1"
    local batch_id
    batch_id="$(basename "${batch_fa}" .fa)"
    local out_dir="${VIBRANT_OUT}/${batch_id}"

    if [[ -d "${out_dir}/VIBRANT_${batch_id}" ]]; then
        echo "[skip] ${batch_id}"
        return 0
    fi

    mkdir -p "${out_dir}"

    VIBRANT_run.py \
        -i "${batch_fa}" \
        -folder "${out_dir}" \
        -t "${BATCH_THREADS}" \
        -no_plot \
        -d "${VIBRANT_DB_DIR}" \
        -m "${VIBRANT_FILES_DIR}" \
        > "${LOG_DIR}/vibrant_${batch_id}.log" 2>&1

    if [[ -d "${out_dir}/VIBRANT_${batch_id}" ]]; then
        echo "[done] ${batch_id}"
    else
        echo "[FAIL] ${batch_id} -- see ${LOG_DIR}/vibrant_${batch_id}.log"
    fi
}
export -f run_one_batch
export VIBRANT_OUT LOG_DIR BATCH_THREADS VIBRANT_DB_DIR VIBRANT_FILES_DIR

find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | sort | \
    parallel -j "${BATCH_JOBS}" run_one_batch {}

n_total=$(find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" | wc -l)
n_done=$(find "${VIBRANT_OUT}" -mindepth 2 -maxdepth 2 -type d -name "VIBRANT_batch_*" | wc -l)
log "Step 1 (VIBRANT) finished: ${n_done}/${n_total} batches completed."
log "Next: 8p_virus_detection.sh checkv"
)

step_checkv() (
source "$(dirname "$0")/config_virome.sh"

log "Steps 2-3: union VirSorter2 + VIBRANT -> CheckV -> quality filter"

UNION_FA="${CONFIRMED_VIRAL_OUT}/all_viral_union.fa"
CONFIRMED_FA="${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa"
SUPPORT_TSV="${CONFIRMED_VIRAL_OUT}/caller_support.tsv"
VS2_SRC="${CONFIRMED_VIRAL_OUT}/.vs2_source_contigs.txt"
VIB_SRC="${CONFIRMED_VIRAL_OUT}/.vib_source_contigs.txt"
mkdir -p "${CONFIRMED_VIRAL_OUT}"

log "Collecting VirSorter2 sequences..."
: > "${UNION_FA}"; : > "${VS2_SRC}"
vs2_batches=0
while IFS= read -r fa; do
    [[ -s "${fa}" ]] || continue
    grep '^>' "${fa}" | sed -e 's/^>//' -e 's/||.*$//' >> "${VS2_SRC}"
    sed -e '/^>/s/||/__vs2/' "${fa}" >> "${UNION_FA}"
    vs2_batches=$((vs2_batches+1))
done < <(find "${VS2_OUT}" -mindepth 2 -maxdepth 2 -name "final-viral-combined.fa" | sort)
vs2_seqs=$(grep -c '^>' "${UNION_FA}" 2>/dev/null || echo 0)
log "  VirSorter2: ${vs2_seqs} sequences from ${vs2_batches} batches"

log "Collecting VIBRANT sequences..."
: > "${VIB_SRC}"
vib_batches=0
while IFS= read -r fa; do
    [[ -s "${fa}" ]] || continue
    grep '^>' "${fa}" | sed -e 's/^>//' -e 's/[ \t].*$//' -e 's/_fragment_[0-9]*$//' >> "${VIB_SRC}"
    awk '/^>/ { sub(/[ \t].*$/, ""); print $0 "__vib"; next } { print }' "${fa}" >> "${UNION_FA}"
    vib_batches=$((vib_batches+1))
done < <(find "${VIBRANT_OUT}" -name "*.phages_combined.fna" | sort)
total_seqs=$(grep -c '^>' "${UNION_FA}" 2>/dev/null || echo 0)
log "  VIBRANT: $((total_seqs - vs2_seqs)) sequences from ${vib_batches} batches"
log "Union set: ${total_seqs} sequences (pre-CheckV)"

if [[ "${total_seqs}" -eq 0 ]]; then
    log "ERROR: union is empty. Did 8p_virus_detection.sh virsorter2 and 8p_virus_detection.sh vibrant complete?"
    exit 1
fi

log "Building caller support table..."
sort -u "${VS2_SRC}" -o "${VS2_SRC}"
sort -u "${VIB_SRC}" -o "${VIB_SRC}"

awk -F'\t' '
    BEGIN { OFS="\t"; print "source_contig","mag_index","vs2","vibrant","support" }
    FNR==NR { vs2[$1]=1; next }
    { vib[$1]=1 }
    END {
        for (c in vs2) all[c]=1
        for (c in vib) all[c]=1
        for (c in all) {
            v = (c in vs2) ? 1 : 0
            b = (c in vib) ? 1 : 0
            s = (v && b) ? "both" : (v ? "vs2_only" : "vibrant_only")
            print c, substr(c,1,7), v, b, s
        }
    }
' "${VS2_SRC}" "${VIB_SRC}" > "${SUPPORT_TSV}"

read -r n_both n_vs2 n_vib < <(awk -F'\t' '
    NR>1 { if($5=="both")a++; else if($5=="vs2_only")b++; else c++ }
    END { print a+0, b+0, c+0 }' "${SUPPORT_TSV}")
log "  caller support: ${n_both} both, ${n_vs2} VS2-only, ${n_vib} VIBRANT-only"
log "Running CheckV end_to_end on the union set..."
checkv end_to_end \
    "${UNION_FA}" \
    "${CHECKV_OUT}" \
    -d "${CHECKV_DB}" \
    -t "${THREADS}" \
    > "${LOG_DIR}/checkv_union.log" 2>&1

QUALITY_SUMMARY="${CHECKV_OUT}/quality_summary.tsv"
if [[ ! -f "${QUALITY_SUMMARY}" ]]; then
    log "ERROR: no quality_summary.tsv -- see ${LOG_DIR}/checkv_union.log"
    exit 1
fi

log "Filtering by CheckV quality tier (Complete/High/Medium)..."
awk -F'\t' '
    NR==1 { for (i=1;i<=NF;i++){ if($i=="contig_id")id=i; if($i=="checkv_quality")q=i }
            if(!id||!q){ print "ERROR: missing contig_id/checkv_quality" > "/dev/stderr"; exit 1 }
            next }
    $q=="Complete" || $q=="High-quality" || $q=="Medium-quality" { print $id }
' "${QUALITY_SUMMARY}" > "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt"

n_confirmed=$(wc -l < "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt")
log "  ${n_confirmed} sequences pass the filter"

CHECKV_POOL="${CONFIRMED_VIRAL_OUT}/.checkv_pool.fna"
: > "${CHECKV_POOL}"
[[ -s "${CHECKV_OUT}/viruses.fna"    ]] && cat "${CHECKV_OUT}/viruses.fna"    >> "${CHECKV_POOL}"
[[ -s "${CHECKV_OUT}/proviruses.fna" ]] && cat "${CHECKV_OUT}/proviruses.fna" >> "${CHECKV_POOL}"
[[ -s "${CHECKV_POOL}" ]] || cp "${UNION_FA}" "${CHECKV_POOL}"
log "  extraction pool: $(grep -c '^>' "${CHECKV_POOL}") sequences (viruses.fna + proviruses.fna)"

if ! command -v seqkit >/dev/null 2>&1; then
    log "ERROR: seqkit is required for prefix-match extraction (exact-match"
    log "       fallbacks have been shown to silently drop ~30% of provirus"
    exit 1
fi

sed 's/^/^/' "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt" > "${CONFIRMED_VIRAL_OUT}/.confirmed_ids_regex.txt"

for attempt in 1 2; do
    seqkit grep -n -r -f "${CONFIRMED_VIRAL_OUT}/.confirmed_ids_regex.txt" "${CHECKV_POOL}" > "${CONFIRMED_FA}"
    n_final=$(grep -c '^>' "${CONFIRMED_FA}" 2>/dev/null || echo 0)
    if [[ "${n_final}" -ge "${n_confirmed}" ]]; then
        break
    fi
    if [[ "${attempt}" -eq 1 ]]; then
        log "  WARNING: attempt 1 recovered only ${n_final}/${n_confirmed} -- possible filesystem race, sleeping 10s and retrying"
        sleep 10
        sync
    fi
done
log "  extracted ${n_final} sequences -> ${CONFIRMED_FA}"
if [[ "${n_final}" -lt "${n_confirmed}" ]]; then
    log "  WARNING: recovered ${n_final} of ${n_confirmed} expected -- investigate before proceeding"
elif [[ "${n_final}" -gt $((n_confirmed * 12 / 10)) ]]; then
    log "  WARNING: recovered ${n_final}, notably more than ${n_confirmed} expected -- check for over-matching"
    log "  spot-check: grep '^>' ${CONFIRMED_FA} | sed -E 's/_[0-9]+\$//' | sort | uniq -c | sort -rn | head"
fi

awk -F'\t' 'NR>1 && $5=="both" {print $1}' "${SUPPORT_TSV}" | sort -u > "${CONFIRMED_VIRAL_OUT}/.both_contigs.txt"
awk 'NR==FNR { both[$1]=1; next }
     { ctg=$1
       sub(/__vs2.*$/,"",ctg)
       sub(/__vib.*$/,"",ctg)
       sub(/_fragment_[0-9]*$/,"",ctg)
       if (ctg in both) print $1 }' \
    "${CONFIRMED_VIRAL_OUT}/.both_contigs.txt" "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt" \
    > "${CONFIRMED_VIRAL_OUT}/confirmed_both_ids.txt"

n_both_final=$(wc -l < "${CONFIRMED_VIRAL_OUT}/confirmed_both_ids.txt")
log "  of which ${n_both_final} have two-caller support -> confirmed_both_ids.txt"

rm -f "${VS2_SRC}" "${VIB_SRC}" "${CHECKV_POOL}" "${CONFIRMED_VIRAL_OUT}/.both_contigs.txt" \
      "${CONFIRMED_VIRAL_OUT}/.confirmed_ids_regex.txt"

log "Steps 2-3 finished."
log "  primary set     : ${CONFIRMED_FA} (${n_final} seqs, union)"
log "  sensitivity set : confirmed_both_ids.txt (${n_both_final} seqs, two-caller)"
log "  source MAG      : first 7 chars of header, join via ${MAG_INDEX_TSV}"
log "Next: 9p_virus_annotation.sh dramv and 10p_votus.sh cluster"
)

for s in ${@:-batches virsorter2 vibrant checkv}; do
  echo "==> 8p_virus_detection.sh: $s"
  "step_$s"
done
