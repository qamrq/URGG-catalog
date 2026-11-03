#!/usr/bin/env python3
# ============================================================
# 13c_remap_hosts_to_species.py
#
# Collapses the 21 raw `host` metadata values (breed/production-system
# granularity, e.g. Beef_cattle/Dairy_cattle/Boran_cattle/Cow all being
# Bos taurus or Bos indicus) into the 17 true host species used in the
# manuscript's Family/Species table, using COMMON NAMES for the website
# rather than Latin binomials.
#
# Re-aggregates the ALREADY-COMPUTED per_host_species_profile.tsv — does
# NOT re-run Steps 1-4, since nothing about the underlying joins changes,
# only how rows get grouped afterward.
#
# *** FLAGGED FOR YOUR CONFIRMATION BEFORE TRUSTING THIS ***
#   - "Gaur" (raw label) mapped to Bos frontalis per your table, but Gaur
#     colloquially = wild Bos gaurus, a DIFFERENT species from the
#     domesticated Bos frontalis (gayal/mithun) your table lists. Confirm
#     which is actually correct for your samples.
#   - Bactrian vs Dromedary camel cannot be split — raw metadata only has
#     one generic "Camel" label. Both your table's Camelidae rows collapse
#     into one "Camel" category unless you have species-level info elsewhere.
# ============================================================
import pandas as pd
from pathlib import Path

IN_TSV  = Path("/disk252/qamar/metagenome_analysis/annotation/website_per_host/per_host_species_profile.tsv")
OUT_DIR = Path("/disk252/qamar/metagenome_analysis/annotation/website_per_host_by_species")
OUT_DIR.mkdir(exist_ok=True)

# raw host label (as it appears in metadata.tsv) -> (common name, Latin binomial)
HOST_MAP = {
    "Cattle":               ("Cattle",              "Bos taurus"),
    "Beef_cattle":          ("Cattle",              "Bos taurus"),
    "Beef_Heifers":         ("Cattle",              "Bos taurus"),
    "Dairy_cattle":         ("Cattle",              "Bos taurus"),
    "Cow":                  ("Cattle",              "Bos taurus"),
    "Boran_cattle":         ("Zebu cattle",         "Bos indicus"),
    "Nellore_cattle":       ("Zebu cattle",         "Bos indicus"),
    "Yak":                  ("Yak",                 "Bos grunniens"),
    "Gaur":                 ("Gaur",                "Bos frontalis"),   # CONFIRM — see caveat above
    "Sheep":                ("Sheep",               "Ovis aries"),
    "Goat":                 ("Goat",                "Capra hircus"),
    "Buffalo":              ("Buffalo",             "Bubalus bubalis"),
    "Antelope":             ("Tibetan antelope",    "Pantholops hodgsonii"),
    "Gazelle":              ("Tibetan gazelle",     "Procapra picticaudata"),
    "Camel":                ("Camel",               "Camelus sp."),     # unresolved — see caveat above
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

# ---- Re-aggregate: same final_genome can now appear under a merged
# common_name if it was contributed by e.g. both "Cattle" and "Beef_cattle"
# rows for the same species cluster — sum genome counts, keep one profile row.
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
