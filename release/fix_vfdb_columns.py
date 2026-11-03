#!/usr/bin/env python3
"""Replace the VFDB column(s) of web_data/web_genome_function.tsv with one column per VFDB level
(genome, level, trait in 64_vf_presence_propagated.tsv). Other columns and files are untouched."""
import re, pandas as pd
B  = "/disk252/qamar/metagenome_analysis"
WR = f"{B}/website_release"
VF = f"{B}/annotation/results_phase3/s64_vfdb_propagated/64_vf_presence_propagated.tsv"
GF = f"{WR}/web_data/web_genome_function.tsv"
tq = pd.read_csv(f"{WR}/urgg_tax_QC_table.tsv", sep="\t", usecols=["original_genome_id", "public_id"])
o2u = dict(zip(tq.original_genome_id, tq.public_id))
gf = pd.read_csv(GF, sep="\t", dtype=str, keep_default_na=False)
gf = gf[[c for c in gf.columns if not c.startswith("vfdb_")]]
v = pd.read_csv(VF, sep="\t", dtype=str); v["u"] = v.genome.map(o2u); v = v[v.u.notna()]
amr_pos = gf.columns.get_loc("amr_genes")
for lv, sub in v.groupby("level"):
    key = "vfdb_" + re.sub(r"[^A-Za-z0-9]+", "_", lv).strip("_").lower()
    d = sub.groupby("u").trait.apply(lambda s: ";".join(sorted(set(s)))).to_dict()
    gf.insert(amr_pos, key, gf.urgg_id.map(d).fillna("")); amr_pos += 1
    print(f"{key}: {sub.trait.nunique():,} distinct values, {len(d):,} genomes")
gf.to_csv(GF, sep="\t", index=False)
print("updated", GF, "columns:", list(gf.columns))
