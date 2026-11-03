#!/usr/bin/env bash
###############################################################################
# 05_prophage_flagging.sh — Phase 5: Virome pipeline, Step 5
#
# Flags confirmed viral sequences as PROPHAGES using two criteria:
#   (a) CheckV's "provirus" call (Yes/No), OR
#   (b) a gene hit to a VOGDB entry whose functional description indicates
#       mobile lysogeny machinery (integrase, excisionase, recombinase,
#       repressor, lysogen).
#
# ROOT CAUSE HISTORY (extensive debugging, kept for the record):
#   DRAM's per-batch annotations.tsv has a REAL, intentional unnamed leading
#   column (the pandas row index -- the true gene id) -- confirmed via raw
#   hex byte inspection (the header genuinely starts with a literal tab
#   byte). The column literally named "fasta" does NOT hold the gene's own
#   identifier -- it holds the name of the INPUT FASTA FILE that gene came
#   from (e.g. "batch_0001_filtered"). The correct column for parent-contig
#   identity is "scaffold", which already holds the contig-level identifier
#   directly (no per-gene suffix to strip). Every earlier attempt in this
#   script's history joined on "fasta" (wrong column, wrong semantics) or on
#   a mis-corrected header alignment (wrong theory) -- both produced zero or
#   garbage matches despite passing every intermediate sanity check, because
#   the actual bug was a misunderstanding of DRAM's own column semantics,
#   not a parsing or alignment bug. Verified against real data before
#   deploying this version.
#
# Uses the per-batch annotations.tsv files directly (verified clean) rather
# than any previously-built combined file, to avoid any residual risk from
# earlier merge steps.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
mkdir -p "${RESULTS_DIR}"
log "Step 5: prophage flagging (CheckV provirus + VOGDB gene hits)"

QUALITY_SUMMARY="${CHECKV_OUT}/quality_summary.tsv"
CONFIRMED_IDS="${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt"
VOTU_INDEX="${VOTU_OUT}/vOTU_index.tsv"
VOG_ANNOT_GZ="/disk252/qamar/databases/dramv_config/database_download/vog.annotations.tsv.gz"

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