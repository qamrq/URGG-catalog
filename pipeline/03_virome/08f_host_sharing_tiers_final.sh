#!/usr/bin/env bash
###############################################################################
# 08f_host_sharing_tiers_final.sh — Phase 5: standalone vOTU host-sharing
# tier file, using the CORRECTED fixed-threshold definition
#
# cross_host_votu_matrix.tsv's own "category" column uses the ORIGINAL
# strict "present in every host group" definition (core = all 15 species),
# which produced core=0 -- an artifact of requiring universal presence
# across a species panel this size, not a real biological null. The
# manuscript instead uses the same fixed-threshold scheme as the bacterial
# catalog (1 species = host-specific; 2-3 = shared; >3 = core/broad),
# recomputed inline wherever needed since (08c, 08d, manuscript tables).
# This script makes that corrected classification a real, standalone file
# instead of something recomputed ad hoc each time.
#
# Output: RESULTS_DIR/votu_host_sharing_tiers.tsv
#         votu_id | representative | n_hosts | host_sharing_tier | host_species_list
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
MATRIX="${RESULTS_DIR}/cross_host_votu_matrix.tsv"

if [[ ! -s "${MATRIX}" ]]; then
    log "ERROR: ${MATRIX} not found -- run 08a_generate_results_tables.sh first."
    exit 1
fi

log "Building standalone vOTU host-sharing tier file (corrected fixed-threshold scheme)..."

python3 << PYEOF
import csv

votu_id_to_rep = {}
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        votu_id_to_rep[row["votu_id"]] = row["representative"]

out_path = "${RESULTS_DIR}/votu_host_sharing_tiers.tsv"
n_specific = n_shared = n_core = 0

with open("${MATRIX}") as fin, open(out_path, "w", newline="") as fout:
    r = csv.DictReader(fin, delimiter='\t')
    host_cols = [c for c in r.fieldnames if c not in ("votu_id", "n_hosts", "category")]
    w = csv.writer(fout, delimiter='\t')
    w.writerow(["votu_id", "representative", "n_hosts", "host_sharing_tier", "host_species_list"])

    for row in r:
        n_hosts = int(row["n_hosts"])
        if n_hosts == 1:
            tier = "host_specific"; n_specific += 1
        elif n_hosts <= 3:
            tier = "shared"; n_shared += 1
        else:
            tier = "core_broad"; n_core += 1

        present_species = [h for h in host_cols if row[h] == "1"]
        rep = votu_id_to_rep.get(row["votu_id"], "NA")

        w.writerow([row["votu_id"], rep, n_hosts, tier, ";".join(present_species)])

total = n_specific + n_shared + n_core
print(f"=== vOTU HOST-SHARING TIERS (corrected fixed-threshold scheme) ===")
print(f"Host-specific (1 species) : {n_specific} ({100*n_specific/total:.1f}%)")
print(f"Shared (2-3 species)      : {n_shared} ({100*n_shared/total:.1f}%)")
print(f"Core/broad (>3 species)   : {n_core} ({100*n_core/total:.1f}%)")
print(f"Total                     : {total}")
print(f"-> {out_path}")
PYEOF

log "Step 08f finished."
