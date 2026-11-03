#!/usr/bin/env bash
###############################################################################
# 08a_generate_results_tables.sh — Phase 5: manuscript results tables
#
# Extracts summary tables from existing pipeline output (steps 1-6, already
# complete). Pure extraction/aggregation, no tool reruns, safe to rerun any
# time.
#
# FIX (this version): the cross-host analysis (item 4) previously used the
# RAW per-BioProject "animal" label directly (mag_id's second __-delimited
# field), which is NOT a normalized species label -- different BioProjects
# used different conventions for the SAME species (e.g. "Dairy_cattle",
# "Beef_cattle", "Nellore_cattle", "Boran_cattle", "Cattle", "Cow" are all
# Bos taurus). Using raw labels fragmented one species into 7+ separate
# groups, making "core" (present in ALL groups) essentially unreachable and
# artificially inflating "host_specific" -- a labeling artifact that
# produced core=0 on a real run, not a genuine biological result. Fixed by
# consolidating through host_species_map.py before classification.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
mkdir -p "${RESULTS_DIR}"
log "Generating manuscript results tables -> ${RESULTS_DIR}"

###############################################################################
# 1. Viral discovery summary table
###############################################################################
log "[1/6] Discovery summary table..."
n_mags=$(find "${MAGS_99ANI_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l)
n_vs2=$(cat "${VS2_OUT}"/*/final-viral-score.tsv 2>/dev/null | grep -v "^seqname" | wc -l)
n_vib=$(cat "${VIBRANT_OUT}"/*/VIBRANT_batch_*/VIBRANT_phages_batch_*/*.phages_combined.fna 2>/dev/null | grep -c '^>')
n_union=$(grep -c '^>' "${CONFIRMED_VIRAL_OUT}/all_viral_union.fa" 2>/dev/null || echo NA)
n_confirmed=$(grep -c '^>' "${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa" 2>/dev/null || echo NA)
n_votus=$(($(wc -l < "${VOTU_OUT}/vOTU_index.tsv") - 1))

{
  printf 'stage\tcount\n'
  printf 'MAGs_screened\t%s\n' "${n_mags}"
  printf 'VirSorter2_calls\t%s\n' "${n_vs2}"
  printf 'VIBRANT_calls\t%s\n' "${n_vib}"
  printf 'Union_pre_filter\t%s\n' "${n_union}"
  printf 'CheckV_confirmed\t%s\n' "${n_confirmed}"
  printf 'Nonredundant_vOTUs\t%s\n' "${n_votus}"
} > "${RESULTS_DIR}/table1_discovery_summary.tsv"
cat "${RESULTS_DIR}/table1_discovery_summary.tsv"

###############################################################################
# 2. CheckV quality tier distribution
###############################################################################
log "[2/6] CheckV quality tier distribution..."
awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="checkv_quality") q=i; next }
    { print $q }
' "${CHECKV_OUT}/quality_summary.tsv" | sort | uniq -c | sort -rn \
  | awk '{print $2"\t"$1}' > "${RESULTS_DIR}/checkv_quality_tiers.tsv"
sed -i '1i quality_tier\tcount' "${RESULTS_DIR}/checkv_quality_tiers.tsv"
cat "${RESULTS_DIR}/checkv_quality_tiers.tsv"

###############################################################################
# 3. vOTU size distribution
###############################################################################
log "[3/6] vOTU size distribution..."
awk -F'\t' 'NR>1 {print $3}' "${VOTU_OUT}/vOTU_index.tsv" \
  | sort -n | uniq -c | awk '{print $2"\t"$1}' \
  > "${RESULTS_DIR}/votu_size_distribution.tsv"
sed -i '1i cluster_size\tn_votus' "${RESULTS_DIR}/votu_size_distribution.tsv"
n_singletons=$(awk -F'\t' 'NR>1 && $3==1' "${VOTU_OUT}/vOTU_index.tsv" | wc -l)
log "  ${n_singletons} / $(($(wc -l < "${VOTU_OUT}/vOTU_index.tsv")-1)) vOTUs are singletons"

###############################################################################
# 4. Cross-host vOTU distribution -- FIXED to consolidate host labels to
#    canonical species before classification (see host_species_map.py)
###############################################################################
log "[4/6] Cross-host vOTU distribution (key comparative result)..."
python3 << PYEOF
import csv, sys
from collections import defaultdict

sys.path.insert(0, "/home/claude/phase5_scripts")
try:
    from host_species_map import normalize_host
except ImportError:
    # fallback inline copy if the module file isn't found on this path --
    # keep in sync with host_species_map.py
    RAW_LABEL_TO_SPECIES = {
        "Cattle": "Bos_taurus", "Cow": "Bos_taurus", "Beef_cattle": "Bos_taurus",
        "Beef_Heifers": "Bos_taurus", "Dairy_cattle": "Bos_taurus",
        "Nellore_cattle": "Bos_taurus", "Boran_cattle": "Bos_taurus",
        "Yak": "Bos_grunniens", "Gaur": "Bos_frontalis", "Buffalo": "Bubalus_bubalis",
        "Sheep": "Ovis_aries", "Goat": "Capra_hircus",
        "Antelope": "Pantholops_hodgsonii", "Gazelle": "Procapra_picticaudata",
        "Camel": "Camelus_dromedarius", "Giraffe": "Giraffa_camelopardalis",
        "Moose": "Alces_alces", "North_American_Moose": "Alces_americanus",
        "Red_deer": "Cervus_elaphus", "Reindeer": "Rangifer_tarandus",
        "Water_deer": "Hydropotes_inermis",
    }
    def normalize_host(raw_label):
        return RAW_LABEL_TO_SPECIES.get(raw_label, f"UNMAPPED_{raw_label}")

mag_to_host = {}
unmapped_seen = set()
with open("${MAG_INDEX_TSV}") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        mag_id = row["mag_id"]
        parts = mag_id.split("__")
        raw_label = parts[1] if len(parts) > 1 else "UNKNOWN"
        host = normalize_host(raw_label)
        if host.startswith("UNMAPPED_"):
            unmapped_seen.add(raw_label)
        mag_to_host[row["mag_index"]] = host

if unmapped_seen:
    print(f"WARNING: {len(unmapped_seen)} raw label(s) have NO species mapping "
          f"and were left as UNMAPPED_<label> (treated as their own group -- "
          f"add these to host_species_map.py for a correct result): {sorted(unmapped_seen)}")

all_hosts = sorted(set(mag_to_host.values()))
print(f"Detected {len(all_hosts)} CONSOLIDATED host species: {all_hosts}")

votu_hosts = defaultdict(set)
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        votu_id = row["votu_id"]
        members = row["member_ids"].split(",")
        for m in members:
            mag_idx = m.strip()[:7]
            host = mag_to_host.get(mag_idx, "UNKNOWN")
            votu_hosts[votu_id].add(host)

matrix_path = "${RESULTS_DIR}/cross_host_votu_matrix.tsv"
with open(matrix_path, "w", newline="") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(["votu_id"] + all_hosts + ["n_hosts", "category"])
    n_core = n_shared = n_specific = 0
    for votu_id, hosts in votu_hosts.items():
        presence = [1 if h in hosts else 0 for h in all_hosts]
        n_h = len(hosts)
        if n_h == len(all_hosts):
            cat = "core"; n_core += 1
        elif n_h > 1:
            cat = "shared"; n_shared += 1
        else:
            cat = "host_specific"; n_specific += 1
        w.writerow([votu_id] + presence + [n_h, cat])

summary_path = "${RESULTS_DIR}/cross_host_summary.tsv"
with open(summary_path, "w") as f:
    f.write("category\tn_votus\n")
    f.write(f"core\t{n_core}\n")
    f.write(f"shared\t{n_shared}\n")
    f.write(f"host_specific\t{n_specific}\n")

print(f"core={n_core} shared={n_shared} host_specific={n_specific}")
print(f"-> {matrix_path}")
print(f"-> {summary_path}")
PYEOF

###############################################################################
# 5. Genome length / completeness
###############################################################################
log "[5/6] Length/completeness table..."
awk -F'\t' '
    NR==1 {
        for(i=1;i<=NF;i++){
            if($i=="contig_id") id=i
            if($i=="contig_length") len=i
            if($i=="completeness") comp=i
            if($i=="checkv_quality") q=i
        }
        print "contig_id\tlength\tcompleteness\tquality"; next
    }
    { print $id"\t"$len"\t"$comp"\t"$q }
' "${CHECKV_OUT}/quality_summary.tsv" > "${RESULTS_DIR}/length_completeness_all.tsv"

awk -F'\t' 'NR==FNR{keep[$1]=1;next} FNR==1{print;next} ($1 in keep)' \
    "${CONFIRMED_VIRAL_OUT}/confirmed_ids.txt" "${RESULTS_DIR}/length_completeness_all.tsv" \
    > "${RESULTS_DIR}/length_completeness.tsv"
log "  $(($(wc -l < "${RESULTS_DIR}/length_completeness.tsv")-1)) confirmed entries with length/completeness"

###############################################################################
# 6. AMG summary
#    KNOWN LIMITATION (confirmed on a real run): auxiliary_score was
#    uniformly 5 across all annotated genes -- DRAM-v's AMG scoring needs
#    non-viral flanking gene context to assess genome-boundary position,
#    which the boundary-trimmed confirmed-viral-only input used for
#    annotation does not retain. amg_summary.tsv is expected to be near-
#    empty for this reason; this is a stated scope limitation (deferred to
#    follow-up work with untrimmed contig context), not a bug to chase here.
###############################################################################
log "[6/6] AMG summary across all batches..."
AMG_COMBINED="${RESULTS_DIR}/amg_summary_combined.tsv"
first=1
for f in "${DRAMV_OUT}"/batch_*/distill/amg_summary.tsv; do
    [[ -f "${f}" ]] || continue
    if [[ "${first}" -eq 1 ]]; then
        cp "${f}" "${AMG_COMBINED}"; first=0
    else
        tail -n +2 "${f}" >> "${AMG_COMBINED}"
    fi
done
[[ -f "${DRAMV_OUT}/vibrant_only/distill/amg_summary.tsv" ]] && \
    tail -n +2 "${DRAMV_OUT}/vibrant_only/distill/amg_summary.tsv" >> "${AMG_COMBINED}" 2>/dev/null

n_amgs=$(($(wc -l < "${AMG_COMBINED}" 2>/dev/null || echo 1) - 1))
log "  ${n_amgs} total AMG entries -> ${AMG_COMBINED}"
if [[ "${n_amgs}" -le 0 ]]; then
    log "  EXPECTED (see script comments): AMG scoring lost flanking context during"
    log "  boundary trimming. Treat as a stated methods limitation, not a bug."
fi

log "All 6 result tables written to ${RESULTS_DIR}"