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
#  SCRIPT  : 13c_remap_hosts_to_species.py  -  Collapse host labels to the 17 host species
#  AUTHOR  : QRQ
#  INPUT   : annotation/website_per_host/per_host_species_profile.tsv
#  OUTPUT  : annotation/website_per_host_by_species/per_species_profile.tsv
# ==========================================================

URGG_PROJECT = __import__("os").environ.get("URGG_PROJECT", "/path/to/metagenome_analysis")

import pandas as pd
from pathlib import Path

IN_TSV  = Path((URGG_PROJECT + "/annotation/website_per_host/per_host_species_profile.tsv"))
OUT_DIR = Path((URGG_PROJECT + "/annotation/website_per_host_by_species"))
OUT_DIR.mkdir(exist_ok=True)

HOST_MAP = {
    "Cattle":               ("Cattle",              "Bos taurus"),
    "Beef_cattle":          ("Cattle",              "Bos taurus"),
    "Beef_Heifers":         ("Cattle",              "Bos taurus"),
    "Dairy_cattle":         ("Cattle",              "Bos taurus"),
    "Cow":                  ("Cattle",              "Bos taurus"),
    "Boran_cattle":         ("Zebu cattle",         "Bos indicus"),
    "Nellore_cattle":       ("Zebu cattle",         "Bos indicus"),
    "Yak":                  ("Yak",                 "Bos grunniens"),
    "Gaur":                 ("Gaur",                "Bos frontalis"),
    "Sheep":                ("Sheep",               "Ovis aries"),
    "Goat":                 ("Goat",                "Capra hircus"),
    "Buffalo":              ("Buffalo",             "Bubalus bubalis"),
    "Antelope":             ("Tibetan antelope",    "Pantholops hodgsonii"),
    "Gazelle":              ("Tibetan gazelle",     "Procapra picticaudata"),
    "Camel":                ("Camel",               "Camelus sp."),
    "Giraffe":              ("Giraffe",             "Giraffa camelopardalis"),
    "Water_deer":           ("Water deer",          "Hydropotes inermis"),
    "Moose":                ("Moose",               "Alces alces"),
    "North_American_Moose": ("Moose (N. American)", "Alces americanus"),
    "Red_deer":             ("Red deer",            "Cervus elaphus"),
    "Reindeer":             ("Reindeer",            "Rangifer tarandus"),
}

df = pd.read_csv(IN_TSV, sep="\t", dtype=str)
df["n_contributing_genomes"] = pd.to_numeric(df["n_contributing_genomes"], errors="coerce")

unmapped = set(df["host"].unique()) - set(HOST_MAP.keys())
if unmapped:
    print(f"WARNING: {len(unmapped)} raw host label(s) have no mapping and will "
          f"be dropped: {sorted(unmapped)}")
    print("  Add them to HOST_MAP above if they should be included.")

df["common_name"] = df["host"].map(lambda h: HOST_MAP.get(h, (None, None))[0])
df["species_binomial"] = df["host"].map(lambda h: HOST_MAP.get(h, (None, None))[1])
df = df.dropna(subset=["common_name"])

n_species_final = df["common_name"].nunique()
print(f"Collapsed {df['host'].nunique()} raw host labels into {n_species_final} "
      f"common-name categories.")
if n_species_final != 17:
    print(f"  NOTE: this is {n_species_final}, not 17 — check HOST_MAP against "
          f"your actual raw host values and the manuscript's species table.")

profile_cols = [c for c in df.columns
                if c not in ("host", "common_name", "species_binomial",
                             "final_genome", "n_contributing_genomes")]

agg = (df.groupby(["common_name", "species_binomial", "final_genome"])
       .agg({**{"n_contributing_genomes": "sum"},
             **{c: "first" for c in profile_cols}})
       .reset_index())

out_tsv = OUT_DIR / "per_species_profile.tsv"
agg.to_csv(out_tsv, sep="\t", index=False)
print(f"Wrote {len(agg)} rows to {out_tsv}")

for name, sub in agg.groupby("common_name"):
    safe = name.replace(" ", "_").replace("(", "").replace(")", "").replace(".", "")
    (OUT_DIR / f"{safe}.json").write_text(
        sub.drop(columns=["common_name"]).to_json(orient="records", indent=2)
    )

print(f"Wrote {agg['common_name'].nunique()} per-species JSON files to {OUT_DIR}")
print("\nFinal common-name categories:")
for name in sorted(agg["common_name"].unique()):
    n = agg[agg["common_name"] == name]["final_genome"].nunique()
    print(f"  {name}: {n} species clusters present")
