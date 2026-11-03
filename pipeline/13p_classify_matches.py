#!/usr/bin/env python3
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
#  SCRIPT  : 13p_classify_matches.py  -  Classify vOTU matches to a reference virome (called by 13p)
#  AUTHOR  : QRQ
#  INPUT   : ANI table and vOTU metadata (arguments)
#  OUTPUT  : match classification table
# ==========================================================

import argparse
import csv

def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def tier_for(ani, cov, species_ani, species_cov, genus_ani, genus_cov):
    if ani >= species_ani and cov >= species_cov:
        return "species_matched"
    if ani >= genus_ani and cov >= genus_cov:
        return "genus_related"
    return "distant_hit"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query-fa", required=True)
    ap.add_argument("--ani", required=True, help="anicalc output tsv")
    ap.add_argument("--species-ani", type=float, default=95.0)
    ap.add_argument("--species-cov", type=float, default=85.0)
    ap.add_argument("--genus-ani", type=float, default=80.0)
    ap.add_argument("--genus-cov", type=float, default=50.0)
    ap.add_argument("--out-table", required=True, help="full per-vOTU table")
    ap.add_argument("--out-plot-csv", required=True, help="slim CSV ready for plotting")
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    all_query_ids = read_fasta_ids(args.query_fa)

    tier_rank = {"species_matched": 0, "genus_related": 1, "distant_hit": 2}

    # track the best-tier hit per query (tie-broken by ani, then coverage)
    best_hit = {}

    with open(args.ani) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if not row:
                continue
            qid, tid = row[0], row[1]
            ani = float(row[3])
            qcov = float(row[4])
            tcov = float(row[5])
            cov_shorter = max(qcov, tcov)

            tier = tier_for(ani, cov_shorter, args.species_ani, args.species_cov,
                             args.genus_ani, args.genus_cov)
            candidate = (tier_rank[tier], ani, cov_shorter, tid, qcov, tcov, tier)
            current = best_hit.get(qid)
            # lower tier_rank is better; within same tier, higher ani/cov wins
            if current is None or (candidate[0], -candidate[1], -candidate[2]) < (current[0], -current[1], -current[2]):
                best_hit[qid] = candidate

    counts = {"species_matched": 0, "genus_related": 0, "distant_hit": 0, "no_hit": 0}

    with open(args.out_table, "w", newline="") as tab_fh, \
         open(args.out_plot_csv, "w", newline="") as plot_fh:

        writer = csv.writer(tab_fh, delimiter="\t")
        writer.writerow([
            "query_votu_id", "tier", "best_rvd_hit",
            "ani", "query_coverage", "target_coverage", "coverage_of_shorter"
        ])

        plot_writer = csv.writer(plot_fh)
        plot_writer.writerow(["query_votu_id", "ani", "coverage_of_shorter", "tier"])

        for qid in all_query_ids:
            hit = best_hit.get(qid)
            if hit is None:
                tier = "no_hit"
                writer.writerow([qid, tier, "NA", "NA", "NA", "NA", "NA"])
            else:
                _, ani, cov_shorter, tid, qcov, tcov, tier = hit
                writer.writerow([qid, tier, tid, f"{ani:.2f}", f"{qcov:.2f}", f"{tcov:.2f}", f"{cov_shorter:.2f}"])
                plot_writer.writerow([qid, f"{ani:.2f}", f"{cov_shorter:.2f}", tier])

            counts[tier] += 1

    total = sum(counts.values())

    with open(args.out_summary, "w") as out_fh:
        out_fh.write("RVD comparison summary (tiered)\n")
        out_fh.write("================================\n")
        out_fh.write(f"Species tier:  ANI >= {args.species_ani}%, coverage(shorter) >= {args.species_cov}%\n")
        out_fh.write(f"Genus tier:    ANI >= {args.genus_ani}%, coverage(shorter) >= {args.genus_cov}%\n\n")
        out_fh.write(f"Total query vOTUs:                  {total}\n")
        for tier in ["species_matched", "genus_related", "distant_hit", "no_hit"]:
            n = counts[tier]
            pct = 100 * n / total if total else 0.0
            out_fh.write(f"  {tier:<18} {n:>6}  ({pct:5.1f}%)\n")

    print(", ".join(f"{k}: {v} ({100*v/total:.1f}%)" for k, v in counts.items()))

if __name__ == "__main__":
    main()
