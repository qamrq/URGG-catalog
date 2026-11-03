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
#  SCRIPT  : 03_merge_checkv.sh  -  Union of viral calls and quality filtering (CheckV)
#  AUTHOR  : QRQ
#  INPUT   : VirSorter2 and VIBRANT outputs
#  OUTPUT  : confirmed_viral.fa, caller_support.tsv, CheckV quality_summary.tsv
# ==========================================================

source "$(dirname "$0")/config.sh"

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
    log "ERROR: union is empty. Did 01_virsorter2.sh and 02_vibrant.sh complete?"
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

activate_env "${ENV_CHECKV}"
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
    log "       calls). Install with: conda install -n ${ENV_CHECKV} -c bioconda seqkit"
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
log "Next: 04_dramv.sh and 06_votu_clustering.sh"
