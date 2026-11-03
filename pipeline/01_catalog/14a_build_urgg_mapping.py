#!/usr/bin/env python3
# ============================================================
# 14a_build_urgg_mapping.py
#
# Assigns clean public accession IDs (URGG00001 ... URGG25451) to the final
# species catalog. Does NOT touch any files — only writes a mapping table.
#
# Saved at: drep_95_species_final/genome_id_mapping.tsv (sibling to
# dereplicated_genomes/, so it doesn't get picked up by any `find *.fa`
# call elsewhere in the pipeline).
#
# ORDERING: taxonomic order (domain->phylum->class->order->family->genus->
# species, alphabetical within each rank), unclassified genomes sorted last
# within their nearest resolved rank. Deterministic — re-running on the same
# GTDB-Tk output always reproduces the identical mapping.
#
# Host is deliberately NOT part of the ID: many species are shared across
# multiple hosts, and the representative genome's originating host is
# essentially arbitrary (whichever genome won dRep's scoring) — encoding
# host into the ID would misrepresent shared species as host-exclusive.
# Host stays as metadata, joined separately (already done via the website
# per-host JSON work).
# ============================================================
import pandas as pd
import re
from pathlib import Path

ANNOT   = Path("/disk252/qamar/metagenome_analysis/annotation")
DREP    = Path("/disk252/qamar/metagenome_analysis/drep_analysis/drep_95_species_final")
GTDBTK_SUMMARY = ANNOT / "gtdbtk_output_v232" / "gtdbtk.all.summary.tsv"
OUT_FILE = DREP / "genome_id_mapping.tsv"

N_DIGITS = 5     # 25,451 fits in 5 digits
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
    print("\nFirst 5:")
    print(out.head(5).to_string(index=False))
    print("\nOriginal genome files are UNCHANGED. Next: 14b_create_release_copies.sh")


if __name__ == "__main__":
    main()
