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
#  SCRIPT  : 4p_genome_ids.py  -  Assign URGG identifiers to species-level genomes
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_95_species_final/
#            annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#  OUTPUT  : drep_analysis/drep_95_species_final/genome_id_mapping.tsv
# ==========================================================

URGG_PROJECT = __import__("os").environ.get("URGG_PROJECT", "/path/to/metagenome_analysis")

import pandas as pd
import re
from pathlib import Path

ANNOT   = Path((URGG_PROJECT + "/annotation"))
DREP    = Path((URGG_PROJECT + "/drep_analysis/drep_95_species_final"))
GTDBTK_SUMMARY = ANNOT / "gtdbtk_output_v232" / "gtdbtk.all.summary.tsv"
OUT_FILE = DREP / "genome_id_mapping.tsv"

N_DIGITS = 5
PREFIX   = "URGG"

strip_ext = lambda s: str(s).removesuffix(".fa")

def parse_rank(classification: str, code: str) -> str:
    m = re.search(rf"{code}__([^;]*)", str(classification))
    val = m.group(1) if m else ""
    return val if val else "zzz_unclassified"

def main():
    gtdb = pd.read_csv(GTDBTK_SUMMARY, sep="\t", dtype=str)
    gcol = next(c for c in gtdb.columns if "genome" in c.lower())
    ccol = next(c for c in gtdb.columns if c.lower() == "classification")
    gtdb["genome"] = gtdb[gcol].apply(strip_ext)

    for rank, code in [("domain", "d"), ("phylum", "p"), ("class", "c"),
                       ("order", "o"), ("family", "f"), ("genus", "g"), ("species", "s")]:
        gtdb[rank] = gtdb[ccol].apply(lambda x: parse_rank(x, code))

    sort_cols = ["domain", "phylum", "class", "order", "family", "genus", "species", "genome"]
    gtdb_sorted = gtdb.sort_values(sort_cols).reset_index(drop=True)
    gtdb_sorted["public_id"] = [f"{PREFIX}{str(i+1).zfill(N_DIGITS)}" for i in range(len(gtdb_sorted))]

    out = gtdb_sorted[["genome", "public_id", ccol]].rename(
        columns={"genome": "original_genome_id", ccol: "gtdb_classification"})
    out.to_csv(OUT_FILE, sep="\t", index=False)

    print(f"Assigned {len(out)} URGG IDs -> {OUT_FILE}")
    print(f"Range: {out['public_id'].iloc[0]} to {out['public_id'].iloc[-1]}")
    print(out.head(5).to_string(index=False))

if __name__ == "__main__":
    main()
