#!/usr/bin/env bash
###############################################################################
# 04_dramv.sh — Phase 5: Virome pipeline, Step 4 (functional annotation)
#
# Annotates the VS2-associated portion of the confirmed viral set with DRAM-v.
#
# SCOPE (read before running):
#   DRAM-v requires VirSorter2's own per-batch output format (specifically
#   viral-affi-contigs-for-dramv.tab, which encodes viral/non-viral gene calls
#   VirSorter2 made along each contig). VIBRANT-only confirmed sequences have
#   no equivalent file and CANNOT go through `DRAM-v.py annotate` as-is.
#   This script therefore annotates only sequences with vs2 support (support
#   column "vs2_only" or "both" in caller_support.tsv), which is most of the
#   confirmed set (28,926 both + 86,177 vs2_only vs. 11,084 vibrant_only,
#   pre-quality-filter). VIBRANT-only sequences are NOT annotated here --
#   see METHODS_NOTES.md for the options (plain `DRAM.py annotate` without
#   -v, or leave unannotated) and decide before treating annotation coverage
#   as complete.
#
# WHY FILTER BEFORE RUNNING (not annotate everything then subset):
#   DRAM-v is HMM+DIAMOND search per sequence -- comparable cost to CheckV.
#   Running on all ~115k raw VS2 calls instead of the ~15k CheckV-confirmed,
#   quality-passed subset would be roughly 7-8x more compute for sequences
#   that get discarded downstream anyway.
#
# WHY PREFIX MATCH (not exact): identical reasoning to the CheckV extraction
# fix in 03_merge_checkv.sh -- VirSorter2's for-dramv output renames contigs
# by converting "|" to "_" and appending its own suffix (full/partial/lt2gene),
# which does not equal our union script's "__vs2" tagging byte-for-byte. Both
# are DERIVED from the same original contig name, so prefix matching on the
# shared source-contig portion is what actually reconciles them.
#
# Runs per batch (for-dramv files are batch-scoped), in parallel like
# 01/02. Produces one DRAM-v annotation set per batch, then concatenates into
# combined tables.
#
# Inputs :
#   VS2_OUT/<batch_id>/for-dramv/final-viral-combined-for-dramv.fa
#   VS2_OUT/<batch_id>/for-dramv/viral-affi-contigs-for-dramv.tab
#   CONFIRMED_VIRAL_OUT/confirmed_ids.txt   (from 03_merge_checkv.sh)
#   CONFIRMED_VIRAL_OUT/caller_support.tsv  (from 03_merge_checkv.sh)
# Outputs:
#   DRAMV_OUT/<batch_id>/                    per-batch DRAM-v annotate+distill
#   DRAMV_OUT/combined_annotations.tsv       concatenated annotations
#   DRAMV_OUT/combined_distillate/           concatenated distillation summaries
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"
activate_env "${ENV_DRAMV}"

log "Step 4: DRAM-v annotation (VS2-associated confirmed viral sequences only)"

CONFIRMED_IDS="${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt"
SUPPORT_TSV="${CONFIRMED_VIRAL_OUT}/caller_support.tsv"
FILTER_DIR="${DRAMV_OUT}/_filtered_input"
mkdir -p "${FILTER_DIR}"

if [[ ! -f "${CONFIRMED_IDS}" || ! -f "${SUPPORT_TSV}" ]]; then
    log "ERROR: missing ${CONFIRMED_IDS} or ${SUPPORT_TSV} -- run 03_merge_checkv.sh first."
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    log "ERROR: seqkit required (exact-match extraction has been shown to silently"
    log "       drop sequences -- see 03_merge_checkv.sh comments). Install:"
    log "       conda install -n ${ENV_DRAMV} -c bioconda seqkit"
    exit 1
fi

###############################################################################
# 1. Build the list of confirmed source contigs that have VS2 support.
#    caller_support.tsv's source_contig is the ORIGINAL contig name (before
#    any of our tagging or VirSorter2's own suffixing) -- exactly what a
#    prefix match against for-dramv headers needs.
###############################################################################
log "Building VS2-associated confirmed contig list..."
awk -F'\t' '
    NR==FNR { if ($5=="vs2_only" || $5=="both") vs2_ok[$1]=1; next }
    { ctg=$1
      # Reduce a confirmed_id back to its bare source contig name. Order matters:
      # strip the __vs2/__vib tag (and anything CheckV appended after it, e.g. a
      # trailing "_1" segment number) FIRST, then strip a _fragment_N marker that
      # sat before the __vib tag. Do NOT blindly strip a trailing "_[0-9]+$" after
      # that -- a real bug run showed this deletes legitimate zero-padded contig
      # numbers (e.g. "..._contig_000018" -> "..._contig", losing the match
      # entirely) whenever a contig happened to end in digits, which is the
      # normal case for assembler-generated contig names.
      sub(/__vs2.*$/,"",ctg)
      sub(/__vib.*$/,"",ctg)
      sub(/_fragment_[0-9]*$/,"",ctg)
      if (ctg in vs2_ok) print ctg
    }
' "${SUPPORT_TSV}" "${CONFIRMED_IDS}" | sort -u > "${FILTER_DIR}/confirmed_vs2_contigs.txt"

n_targets=$(wc -l < "${FILTER_DIR}/confirmed_vs2_contigs.txt")
log "  ${n_targets} confirmed source contigs have VS2 support -> candidates for DRAM-v"

if [[ "${n_targets}" -eq 0 ]]; then
    log "ERROR: no VS2-associated confirmed contigs found. Check 03_merge_checkv.sh output."
    exit 1
fi

sed 's/^/^/' "${FILTER_DIR}/confirmed_vs2_contigs.txt" > "${FILTER_DIR}/confirmed_vs2_contigs_regex.txt"

###############################################################################
# 2. Per-batch: filter for-dramv fasta + tab to the confirmed contig list,
#    then run DRAM-v annotate + distill on the filtered subset.
###############################################################################
run_one_batch () {
    local batch_id="$1"
    local fordramv_dir="${VS2_OUT}/${batch_id}/for-dramv"
    local fa_in="${fordramv_dir}/final-viral-combined-for-dramv.fa"
    local tab_in="${fordramv_dir}/viral-affi-contigs-for-dramv.tab"
    local out_dir="${DRAMV_OUT}/${batch_id}"

    if [[ ! -s "${fa_in}" || ! -s "${tab_in}" ]]; then
        echo "[skip] ${batch_id} -- no for-dramv output (did 01_virsorter2.sh finish for this batch?)"
        return 0
    fi
    if [[ -f "${out_dir}/annotate/annotations.tsv" ]]; then
        echo "[skip] ${batch_id} already done"
        return 0
    fi

    local fa_filt="${FILTER_DIR}/${batch_id}_filtered.fa"
    local tab_filt="${FILTER_DIR}/${batch_id}_filtered.tab"

    # Prefix-match filter the fasta (same technique verified in 03_merge_checkv.sh)
    seqkit grep -n -r -f "${FILTER_DIR}/confirmed_vs2_contigs_regex.txt" "${fa_in}" > "${fa_filt}" 2>/dev/null

    if [[ ! -s "${fa_filt}" ]]; then
        echo "[skip] ${batch_id} -- no confirmed VS2 contigs in this batch after filtering"
        rm -f "${fa_filt}"
        return 0
    fi

    # Filter the affi-contigs tab to the same retained headers (prefix match
    # on column 1, since gene rows are typically "<contig>_gene_N" style).
    grep '^>' "${fa_filt}" | sed 's/^>//' > "${FILTER_DIR}/${batch_id}_retained_headers.txt"
    awk -F'\t' '
        NR==FNR { keep[$1]=1; next }
        {
            matched = 0
            for (h in keep) { if (index($1, h) == 1) { matched = 1; break } }
            if (matched) print
        }
    ' "${FILTER_DIR}/${batch_id}_retained_headers.txt" "${tab_in}" > "${tab_filt}"

    mkdir -p "${out_dir}"
    DRAM-v.py annotate \
        -i "${fa_filt}" \
        -v "${tab_filt}" \
        -o "${out_dir}/annotate" \
        --threads "${BATCH_THREADS}" \
        > "${LOG_DIR}/dramv_annotate_${batch_id}.log" 2>&1

    if [[ ! -f "${out_dir}/annotate/annotations.tsv" ]]; then
        echo "[FAIL] ${batch_id} annotate -- see ${LOG_DIR}/dramv_annotate_${batch_id}.log"
        return 0
    fi

    DRAM-v.py distill \
        -i "${out_dir}/annotate/annotations.tsv" \
        -o "${out_dir}/distill" \
        > "${LOG_DIR}/dramv_distill_${batch_id}.log" 2>&1

    rm -f "${fa_filt}" "${tab_filt}" "${FILTER_DIR}/${batch_id}_retained_headers.txt"
    echo "[done] ${batch_id}"
}
export -f run_one_batch
export VS2_OUT DRAMV_OUT LOG_DIR FILTER_DIR BATCH_THREADS

find "${BATCH_DIR}" -maxdepth 1 -name "batch_*.fa" -printf "%f\n" | sed 's/\.fa$//' | sort | \
    parallel -j "${BATCH_JOBS}" run_one_batch {}

###############################################################################
# 3. Combine per-batch outputs
###############################################################################
log "Combining per-batch DRAM-v annotations..."
COMBINED_ANNOT="${DRAMV_OUT}/combined_annotations.tsv"
first=1
for f in "${DRAMV_OUT}"/batch_*/annotate/annotations.tsv; do
    [[ -f "${f}" ]] || continue
    if [[ "${first}" -eq 1 ]]; then
        cp "${f}" "${COMBINED_ANNOT}"
        first=0
    else
        tail -n +2 "${f}" >> "${COMBINED_ANNOT}"
    fi
done

n_annotated=$(($(wc -l < "${COMBINED_ANNOT}" 2>/dev/null || echo 1) - 1))
log "Step 4 finished: ${n_annotated} genes annotated across all batches -> ${COMBINED_ANNOT}"
log "VIBRANT-only confirmed sequences were NOT annotated (see script header + METHODS_NOTES.md)"
log "Next: 06_votu_clustering.sh"