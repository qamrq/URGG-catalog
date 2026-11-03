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
#  SCRIPT  : 9p_virus_annotation.sh  -  Viral gene annotation (DRAM-v) and prophage flagging
#  AUTHOR  : QRQ
#  INPUT   : confirmed viral sequences, CheckV quality_summary.tsv
#            $URGG_DB/dramv_config/database_download/vog.annotations.tsv.gz
#  OUTPUT  : DRAM-v annotations
#            phase_5_virome/08_results/prophage_flags.tsv
#  USAGE   : bash 9p_virus_annotation.sh [dramv prophage]   (no argument: all steps in this order)
# ==========================================================

export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

step_dramv() (
source "$(dirname "$0")/config_virome.sh"

log "Step 4: DRAM-v annotation (VS2-associated confirmed viral sequences only)"

CONFIRMED_IDS="${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt"
SUPPORT_TSV="${CONFIRMED_VIRAL_OUT}/caller_support.tsv"
FILTER_DIR="${DRAMV_OUT}/_filtered_input"
mkdir -p "${FILTER_DIR}"

if [[ ! -f "${CONFIRMED_IDS}" || ! -f "${SUPPORT_TSV}" ]]; then
    log "ERROR: missing ${CONFIRMED_IDS} or ${SUPPORT_TSV} -- run 8p_virus_detection.sh checkv first."
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    log "ERROR: seqkit required (exact-match extraction has been shown to silently"
    log "       drop sequences -- see 8p_virus_detection.sh checkv)."
    exit 1
fi

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
    log "ERROR: no VS2-associated confirmed contigs found. Check 8p_virus_detection.sh checkv output."
    exit 1
fi

sed 's/^/^/' "${FILTER_DIR}/confirmed_vs2_contigs.txt" > "${FILTER_DIR}/confirmed_vs2_contigs_regex.txt"

run_one_batch () {
    local batch_id="$1"
    local fordramv_dir="${VS2_OUT}/${batch_id}/for-dramv"
    local fa_in="${fordramv_dir}/final-viral-combined-for-dramv.fa"
    local tab_in="${fordramv_dir}/viral-affi-contigs-for-dramv.tab"
    local out_dir="${DRAMV_OUT}/${batch_id}"

    if [[ ! -s "${fa_in}" || ! -s "${tab_in}" ]]; then
        echo "[skip] ${batch_id} -- no for-dramv output (did 8p_virus_detection.sh virsorter2 finish for this batch?)"
        return 0
    fi
    if [[ -f "${out_dir}/annotate/annotations.tsv" ]]; then
        echo "[skip] ${batch_id} already done"
        return 0
    fi

    local fa_filt="${FILTER_DIR}/${batch_id}_filtered.fa"
    local tab_filt="${FILTER_DIR}/${batch_id}_filtered.tab"

    seqkit grep -n -r -f "${FILTER_DIR}/confirmed_vs2_contigs_regex.txt" "${fa_in}" > "${fa_filt}" 2>/dev/null

    if [[ ! -s "${fa_filt}" ]]; then
        echo "[skip] ${batch_id} -- no confirmed VS2 contigs in this batch after filtering"
        rm -f "${fa_filt}"
        return 0
    fi

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
log "Next: 10p_votus.sh cluster"
)

step_prophage() (
source "$(dirname "$0")/config_virome.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
mkdir -p "${RESULTS_DIR}"
log "Step 5: prophage flagging (CheckV provirus + VOGDB gene hits)"

QUALITY_SUMMARY="${CHECKV_OUT}/quality_summary.tsv"
CONFIRMED_IDS="${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt"
VOTU_INDEX="${VOTU_OUT}/vOTU_index.tsv"
VOG_ANNOT_GZ="${URGG_DB}/dramv_config/database_download/vog.annotations.tsv.gz"

for f in "${QUALITY_SUMMARY}" "${CONFIRMED_IDS}" "${VOTU_INDEX}" "${VOG_ANNOT_GZ}"; do
    if [[ ! -s "${f}" ]]; then
        log "ERROR: missing required input ${f}"
        exit 1
    fi
done

MARKER_VOG_LIST="${RESULTS_DIR}/prophage_marker_vogs_discovered.tsv"
log "Discovering prophage-marker VOGs from ${VOG_ANNOT_GZ}..."
zcat "${VOG_ANNOT_GZ}" | awk -F'\t' '
    NR==1 { next }
    tolower($5) ~ /integrase|excisionase|recombinase|repressor|lysogen/ { print $1"\t"$4"\t"$5 }
' > "${MARKER_VOG_LIST}"
n_marker_vogs=$(wc -l < "${MARKER_VOG_LIST}")
log "  ${n_marker_vogs} marker VOGs discovered"

python3 << PYEOF
import csv, glob, re
from collections import defaultdict

confirmed_ids = set()
with open("${CONFIRMED_IDS}") as f:
    for line in f:
        line = line.strip()
        if line:
            confirmed_ids.add(line)

checkv_provirus = {}
with open("${QUALITY_SUMMARY}") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        cid = row["contig_id"]
        if cid in confirmed_ids:
            checkv_provirus[cid] = (row.get("provirus", "No") == "Yes")
print(f"CheckV quality_summary entries matched to confirmed set: {len(checkv_provirus)}")

marker_vogs = {}
with open("${MARKER_VOG_LIST}") as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 3:
            marker_vogs[parts[0]] = parts[2]
print(f"Using {len(marker_vogs)} discovered marker VOGs")

def reduce_to_base(s):
    for sep in ("__vs2", "__vib"):
        idx = s.find(sep)
        if idx != -1:
            return s[:idx]
    return s

confirmed_roots = {cid: reduce_to_base(cid) for cid in confirmed_ids}

# Scan every per-batch annotations.tsv (verified: real leading tab -> blank
# first header field is the true index/gene-id; "scaffold" column holds the
# real parent-contig identifier -- see script header comments).
vogdb_by_confirmed_id = defaultdict(set)
batch_files = sorted(glob.glob("${DRAMV_OUT}/batch_*/annotate/annotations.tsv"))
vib_file = "${DRAMV_OUT}/vibrant_only/annotate/annotations.tsv"
import os
all_files = batch_files + ([vib_file] if os.path.exists(vib_file) else [])
print(f"Scanning {len(all_files)} annotation files...")

n_gene_hits = 0
for fpath in all_files:
    with open(fpath, newline="") as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        try:
            scaffold_idx = header.index("scaffold")
            vogdb_idx = header.index("vogdb_id")
        except ValueError:
            continue
        for row in reader:
            if len(row) <= max(scaffold_idx, vogdb_idx):
                continue
            vogdb_field = row[vogdb_idx]
            if not vogdb_field:
                continue
            hit_vogs = {v for v in marker_vogs if v in vogdb_field}
            if not hit_vogs:
                continue
            n_gene_hits += 1
            scaffold_val = row[scaffold_idx]
            for cid, root in confirmed_roots.items():
                if scaffold_val.startswith(root):
                    vogdb_by_confirmed_id[cid] |= hit_vogs
                    break

print(f"Genes with a prophage-marker VOGDB hit: {n_gene_hits}")
print(f"Confirmed sequences reached via VOGDB marker: {len(vogdb_by_confirmed_id)}")

out_rows = []
n_checkv_only = n_vogdb_only = n_both = n_prophage_total = 0
for cid in confirmed_ids:
    is_provirus = checkv_provirus.get(cid, False)
    vog_hits = vogdb_by_confirmed_id.get(cid, set())
    has_vog = len(vog_hits) > 0
    is_prophage = is_provirus or has_vog

    if is_prophage:
        n_prophage_total += 1
        if is_provirus and has_vog: n_both += 1
        elif is_provirus: n_checkv_only += 1
        else: n_vogdb_only += 1

    hit_desc = ";".join(sorted(marker_vogs.get(v, v) for v in vog_hits)) if vog_hits else "NA"
    out_rows.append([cid, cid[:7], is_provirus, has_vog, hit_desc, is_prophage])

flags_path = "${RESULTS_DIR}/prophage_flags.tsv"
with open(flags_path, "w", newline="") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(["contig_id", "mag_index", "is_provirus_checkv", "vogdb_prophage_gene",
                "vogdb_hit_description", "is_prophage"])
    w.writerows(out_rows)

n_total_confirmed = len(confirmed_ids)
pct = 100 * n_prophage_total / n_total_confirmed if n_total_confirmed else 0

summary_path = "${RESULTS_DIR}/prophage_summary.tsv"
with open(summary_path, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"n_confirmed_viral_sequences\t{n_total_confirmed}\n")
    f.write(f"n_prophage_flagged\t{n_prophage_total}\n")
    f.write(f"pct_prophage_of_confirmed\t{pct:.1f}\n")
    f.write(f"n_flagged_by_checkv_provirus_only\t{n_checkv_only}\n")
    f.write(f"n_flagged_by_vogdb_gene_only\t{n_vogdb_only}\n")
    f.write(f"n_flagged_by_both_criteria\t{n_both}\n")

print(f"\n=== PROPHAGE SUMMARY (comparable to Yan & Yu 2024 Fig. 1a) ===")
print(f"Confirmed viral sequences evaluated : {n_total_confirmed}")
print(f"Prophage-flagged                    : {n_prophage_total} ({pct:.1f}%)")
print(f"  via CheckV provirus call only     : {n_checkv_only}")
print(f"  via VOGDB marker gene only        : {n_vogdb_only}")
print(f"  via both criteria                 : {n_both}")

prophage_by_contig = {row[0]: row[5] for row in out_rows}
votu_prophage_frac = []
with open("${VOTU_INDEX}") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        rep = row["representative"]
        votu_prophage_frac.append([row["votu_id"], rep, prophage_by_contig.get(rep, None)])

votu_path = "${RESULTS_DIR}/votu_prophage_summary.tsv"
with open(votu_path, "w", newline="") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(["votu_id", "representative", "representative_is_prophage"])
    w.writerows(votu_prophage_frac)

n_votu_prophage = sum(1 for r in votu_prophage_frac if r[2] is True)
n_votu_total = len(votu_prophage_frac)
print(f"\nvOTUs whose representative is prophage-flagged: {n_votu_prophage}/{n_votu_total} ({100*n_votu_prophage/n_votu_total:.1f}%)")
PYEOF

log "Step 5 finished."
)

for s in ${@:-dramv prophage}; do
  echo "==> 9p_virus_annotation.sh: $s"
  "step_$s"
done
