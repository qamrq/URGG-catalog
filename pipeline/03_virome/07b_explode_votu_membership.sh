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
#  SCRIPT  : 07b_explode_votu_membership.sh  -  vOTU membership, one row per sequence
#  AUTHOR  : QRQ
#  INPUT   : vOTU_index.tsv
#  OUTPUT  : phase_5_virome/08_results/votu_membership_full.tsv
# ==========================================================

source "$(dirname "$0")/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
VOTU_INDEX="${VOTU_OUT}/vOTU_index.tsv"

if [[ ! -s "${VOTU_INDEX}" ]]; then
    log "ERROR: ${VOTU_INDEX} not found."
    exit 1
fi

log "Exploding vOTU_index.tsv into one row per sequence..."

python3 << PYEOF
import csv

out_path = "${RESULTS_DIR}/votu_membership_full.tsv"
n_rows = 0
n_clusters = 0

with open("${VOTU_INDEX}") as fin, open(out_path, "w", newline="") as fout:
    r = csv.DictReader(fin, delimiter='\t')
    w = csv.writer(fout, delimiter='\t')
    w.writerow(["sequence_id", "votu_id", "representative", "is_representative", "n_members"])

    for row in r:
        n_clusters += 1
        votu_id = row["votu_id"]
        rep = row["representative"]
        n_members = row["n_members"]
        members = [m.strip() for m in row["member_ids"].split(",") if m.strip()]

        for member in members:
            is_rep = (member == rep)
            w.writerow([member, votu_id, rep, is_rep, n_members])
            n_rows += 1

print(f"Exploded {n_clusters} clusters -> {n_rows} sequence-level rows")
print(f"-> {out_path}")
PYEOF

n_check=$(($(wc -l < "${RESULTS_DIR}/votu_membership_full.tsv") - 1))
log "Sanity check: ${n_check} sequence-level rows written (expect close to 15,301, the confirmed viral sequence count)"
log "Step 07b finished."
