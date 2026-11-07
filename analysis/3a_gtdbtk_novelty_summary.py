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
#  SCRIPT  : 3a_gtdbtk_novelty_summary.py  -  Species- and genus-level novelty from GTDB-Tk
#  AUTHOR  : QRQ
#  INPUT   : annotation/gtdbtk_output/gtdbtk.all.summary.tsv
#  OUTPUT  : annotation/gtdbtk_output/novelty_summary/
# ==========================================================

URGG_PROJECT = __import__("os").environ.get("URGG_PROJECT", "/path/to/metagenome_analysis")

import pandas as pd
from pathlib import Path

ANNOT = Path((URGG_PROJECT + "/annotation"))
GTDBTK_SUMMARY = ANNOT / "gtdbtk_output" / "gtdbtk.all.summary.tsv"
OUT_DIR = ANNOT / "gtdbtk_output" / "novelty_summary"
OUT_DIR.mkdir(exist_ok=True)

def parse_classification(classification: str) -> dict:
    """
    GTDB-Tk classification strings look like:
      d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;
      f__Lachnospiraceae;g__Blautia;s__Blautia sp001517625
    An empty rank (e.g. 's__' with nothing after it) means GTDB-Tk did not
    confidently assign that level - that is the novelty signal.
    """
    ranks = {"d": None, "p": None, "c": None, "o": None, "f": None, "g": None, "s": None}
    if not isinstance(classification, str):
        return ranks
    for part in classification.split(";"):
        part = part.strip()
        if len(part) < 3:
            continue
        prefix, _, value = part.partition("__")
        value = value.strip()
        if prefix in ranks:
            ranks[prefix] = value if value else None
    return ranks

def novelty_tier(ranks: dict) -> str:
    if ranks["s"]:
        return "known_species"
    elif ranks["g"]:
        return "novel_species_known_genus"
    elif ranks["f"]:
        return "novel_genus_known_family"
    elif ranks["o"]:
        return "novel_family_known_order"
    elif ranks["c"]:
        return "novel_order_known_class"
    elif ranks["p"]:
        return "novel_class_known_phylum"
    elif ranks["d"]:
        return "novel_phylum_known_domain"
    else:
        return "unclassified"

def main():
    if not GTDBTK_SUMMARY.exists():
        print(f"ERROR: {GTDBTK_SUMMARY} not found.")
        return

    df = pd.read_csv(GTDBTK_SUMMARY, sep="\t")
    print(f"Loaded {len(df)} genome classifications from GTDB-Tk output")

    parsed = df["classification"].apply(parse_classification)
    for rank in ["d", "p", "c", "o", "f", "g", "s"]:
        df[f"rank_{rank}"] = parsed.apply(lambda r: r[rank])

    df["novelty_tier"] = parsed.apply(novelty_tier)

    keep_cols = ["user_genome", "classification", "novelty_tier",
                 "rank_g", "rank_f"]
    optional_cols = ["fastani_reference", "fastani_ani", "fastani_af",
                      "closest_placement_reference", "closest_placement_ani",
                      "closest_placement_af", "classification_method",
                      "note", "warnings"]
    for c in optional_cols:
        if c in df.columns:
            keep_cols.append(c)

    result = df[keep_cols].copy()
    result.to_csv(OUT_DIR / "taxonomic_novelty_per_genome.tsv", sep="\t", index=False)

    # ---- Summary counts ----
    counts = df["novelty_tier"].value_counts()
    total = len(df)

    print("\n" + "=" * 60)
    print("=" * 60)
    print(f"Total genomes classified: {total}\n")
    for tier, n in counts.items():
        pct = 100 * n / total
        print(f"  {tier:35s} {n:6d}  ({pct:5.1f}%)")

    n_novel_species_or_higher = total - counts.get("known_species", 0)
    n_novel_genus_or_higher = sum(
        counts.get(t, 0) for t in [
            "novel_genus_known_family", "novel_family_known_order",
            "novel_order_known_class", "novel_class_known_phylum",
            "novel_phylum_known_domain", "unclassified"
        ]
    )
    print(f"\nTotal novel at species level or above: {n_novel_species_or_higher} "
          f"({100*n_novel_species_or_higher/total:.1f}%)")
    print(f"Total novel at genus level or above:   {n_novel_genus_or_higher} "
          f"({100*n_novel_genus_or_higher/total:.1f}%)  <- your genus-novelty headline number")

    genus_novel = result[result["novelty_tier"].isin([
        "novel_genus_known_family", "novel_family_known_order",
        "novel_order_known_class", "novel_class_known_phylum",
        "novel_phylum_known_domain"
    ])].sort_values("novelty_tier")
    genus_novel_path = OUT_DIR / "genus_novel_candidates.tsv"
    genus_novel.to_csv(genus_novel_path, sep="\t", index=False)
    print(f"\n{len(genus_novel)} genus-or-higher novel candidates written to:")
    print(f"  {genus_novel_path}")

    print(f"\nFull per-genome table: {OUT_DIR / 'taxonomic_novelty_per_genome.tsv'}")

if __name__ == "__main__":
    main()
