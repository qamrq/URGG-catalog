#!/usr/bin/env python3
"""
86a_build_genome_to_species.py — build the per-genome input table for
86b_species_accumulation.py.

Joins the full-scale genome -> final species cluster mapping (117,145 genomes)
to metadata.tsv.  Fixes the earlier join bug: the mapping file stores genome
names WITH ".fa", metadata.tsv stores them WITHOUT, so a direct join gave NA
hosts for every genome.

    python3 86a_build_genome_to_species.py \
        --mapping  checks/fullscale_genome_cluster_mapping.csv \
        --metadata metadata.tsv \
        --out      genome_to_species.tsv

Output columns:
  genome, species_rep, host (common-name label), host_species (Latin, from the
  metadata 'species' column), bioproject, sample_site, sample (SRR/ERR/DRR/CRR
  run accession parsed from the genome name; empty if none)

Also writes <out stem>_sample_id_coverage.tsv: per BioProject, how many genomes
carry a run accession, so you can see which projects drop out of the
sample-based curves.
"""
import argparse, os, re, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--mapping", required=True, help="fullscale_genome_cluster_mapping.csv")
ap.add_argument("--mapping-sep", default=",")
ap.add_argument("--mapping-genome-col", default="genome")
ap.add_argument("--mapping-cluster-col", default="cluster_primary",
                help="cluster_primary = two-stage dRep species; cluster_galah for the Galah check")
ap.add_argument("--metadata", required=True, help="metadata.tsv")
ap.add_argument("--out", default="genome_to_species.tsv")
ap.add_argument("--expect-genomes", type=int, default=117145)
ap.add_argument("--expect-species", type=int, default=25451)
ap.add_argument("--no-strict", action="store_true", help="warn instead of stopping on count mismatches")
a = ap.parse_args()


def fail(msg):
    if a.no_strict:
        print("WARNING: " + msg)
    else:
        sys.exit("ERROR: " + msg + "  (re-run with --no-strict to continue anyway)")


def clean_name(s):
    """basename, strip CR/whitespace, drop a trailing .fa/.fasta/.fna"""
    s = s.astype(str).str.replace("\r", "", regex=False).str.strip()
    s = s.map(os.path.basename)
    return s.str.replace(r"\.(fa|fasta|fna)$", "", regex=True)


# ---- mapping
m = pd.read_csv(a.mapping, sep=a.mapping_sep, dtype=str)
for c in (a.mapping_genome_col, a.mapping_cluster_col):
    if c not in m.columns:
        sys.exit(f"ERROR: column '{c}' not in {a.mapping}. Found: {list(m.columns)}")
m = pd.DataFrame({
    "genome": clean_name(m[a.mapping_genome_col]),
    "species_rep": m[a.mapping_cluster_col].astype(str).str.strip()
                   .str.replace(r"^[A-Za-z_]+::", "", regex=True),   # drop "primary::"
})
if m["genome"].duplicated().any():
    fail(f"{m['genome'].duplicated().sum()} duplicated genome names in mapping")
print(f"mapping : {len(m):,} genomes, {m['species_rep'].nunique():,} species clusters")

# ---- metadata
md = pd.read_csv(a.metadata, sep="\t", dtype=str)
need = ["genome", "bioproject", "host", "species", "sample_site"]
miss = [c for c in need if c not in md.columns]
if miss:
    sys.exit(f"ERROR: metadata missing {miss}. Found: {list(md.columns)}")
md = md[need].copy()
for c in need:
    md[c] = md[c].astype(str).str.replace("\r", "", regex=False).str.strip()
md["genome"] = clean_name(md["genome"])
if md["genome"].duplicated().any():
    fail(f"{md['genome'].duplicated().sum()} duplicated genome names in metadata")
md = md.drop_duplicates("genome").rename(columns={"species": "host_species"})

# ---- join
g = m.merge(md, on="genome", how="left", indicator=True)
n_unmatched = int((g["_merge"] != "both").sum())
print(f"joined  : {len(g):,} rows; unmatched to metadata: {n_unmatched:,}")
if n_unmatched:
    print("  examples:", g.loc[g["_merge"] != "both", "genome"].head(5).tolist())
    fail(f"{n_unmatched} genomes in mapping have no metadata row")
g = g.drop(columns="_merge")

# ---- run accession (sample unit). NA when the name carries no accession;
#      those genomes stay in the genome-based curves but not the sample-based ones.
acc = re.compile(r"(?<![A-Za-z0-9])((?:SRR|ERR|DRR|CRR)\d{5,})")
g["sample"] = g["genome"].str.extract(acc, expand=False)

# ---- checks
if len(g) != a.expect_genomes:
    fail(f"expected {a.expect_genomes:,} genomes, got {len(g):,}")
if g["species_rep"].nunique() != a.expect_species:
    fail(f"expected {a.expect_species:,} species clusters, got {g['species_rep'].nunique():,}")
for c in ("host", "host_species"):
    bad = g[c].isna() | g[c].isin(["", "nan", "NA", "Unknown", "unknown"])
    if bad.any():
        print(f"WARNING: {bad.sum():,} genomes with missing/unknown '{c}'")

print(f"\nhost labels  : {g['host'].nunique()}  |  host species: {g['host_species'].nunique()}")
print(g.groupby("host_species")["host"].agg(lambda v: ", ".join(sorted(set(v)))).to_string())
n_acc = g["sample"].notna().sum()
print(f"\ngenomes with a run accession: {n_acc:,} / {len(g):,} ({100*n_acc/len(g):.1f}%), "
      f"{g['sample'].nunique():,} distinct runs")

cov = (g.assign(has=g["sample"].notna())
         .groupby(["bioproject", "host_species"])
         .agg(genomes=("genome", "size"), with_accession=("has", "sum"),
              runs=("sample", "nunique"))
         .reset_index())
cov["pct_with_accession"] = (100 * cov["with_accession"] / cov["genomes"]).round(1)

out_cols = ["genome", "species_rep", "host", "host_species", "bioproject", "sample_site", "sample"]
g[out_cols].to_csv(a.out, sep="\t", index=False)
stem = re.sub(r"\.(tsv|txt|csv)$", "", a.out)
cov.sort_values("pct_with_accession").to_csv(f"{stem}_sample_id_coverage.tsv", sep="\t", index=False)
print(f"\nWritten: {a.out}\n         {stem}_sample_id_coverage.tsv")
