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
#  SCRIPT  : 13b_build_per_host_tables.py  -  Per-host functional and taxonomic profiles
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/metadata.tsv, stage1/stage2 Cdb.csv
#            annotation/gene_to_genome.tsv, gene catalog clusters
#            annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#  OUTPUT  : annotation/website_per_host/per_host_species_profile.tsv
# ==========================================================

URGG_PROJECT = __import__("os").environ.get("URGG_PROJECT", "/path/to/metagenome_analysis")

import pandas as pd
from pathlib import Path

WORKDIR = Path(URGG_PROJECT)
DREP    = WORKDIR / "drep_analysis"
ANNOT   = WORKDIR / "annotation"
OUT_DIR = ANNOT / "website_per_host"
OUT_DIR.mkdir(exist_ok=True)

METADATA_FILE   = DREP / "metadata.tsv"
STAGE1_DIR      = DREP / "stage1_per_host"
STAGE2_CDB      = DREP / "drep_95_species_final" / "data_tables" / "Cdb.csv"
GENE2GENOME     = ANNOT / "gene_to_genome.tsv"
CLUSTER_TABLE   = ANNOT / "gene_catalog" / "nr_gene_catalog_cluster.tsv"
GTDBTK_SUMMARY  = ANNOT / "gtdbtk_output_v232" / "gtdbtk.all.summary.tsv"
EGGNOG_ANNOT    = ANNOT / "eggnog_output" / "nr_catalog.emapper.annotations"
DBCAN_OVERVIEW  = ANNOT / "dbcan_output" / "overview.tsv"
VFDB_HITS       = ANNOT / "diamond_vfdb" / "diamond_vfdb_filtered.tsv"
AMR_MERGED      = ANNOT / "amrfinderplus_merged.tsv"

strip_ext = lambda s: str(s).removesuffix(".fa")

def build_host_attribution() -> pd.DataFrame:
    """genome (any host, any subgroup) -> final_genome (25,451-catalog rep)."""
    rows = []
    for host_dir in sorted(STAGE1_DIR.glob("*_drep")):
        cdb = host_dir / "data_tables" / "Cdb.csv"
        wdb = host_dir / "data_tables" / "Wdb.csv"
        if not cdb.exists() or not wdb.exists():
            print(f"WARNING: missing tables for {host_dir.name}, skipping")
            continue
        cdb_df = pd.read_csv(cdb)
        wdb_df = pd.read_csv(wdb)
        s1_col = "secondary_cluster" if "secondary_cluster" in cdb_df.columns else "cluster"
        win_col = "cluster" if "cluster" in wdb_df.columns else s1_col

        cdb_df["genome"] = cdb_df["genome"].apply(strip_ext)
        wdb_df["genome"] = wdb_df["genome"].apply(strip_ext)
        winner_map = dict(zip(wdb_df[win_col], wdb_df["genome"]))

        for _, row in cdb_df.iterrows():
            rows.append({
                "genome": row["genome"],
                "stage1_winner": winner_map.get(row[s1_col]),
            })
    stage1_df = pd.DataFrame(rows)

    stage2_df = pd.read_csv(STAGE2_CDB)
    stage2_df["genome"] = stage2_df["genome"].apply(strip_ext)
    s2_col = "secondary_cluster" if "secondary_cluster" in stage2_df.columns else "cluster"

    stage2_wdb_path = STAGE2_CDB.parent / "Wdb.csv"
    if stage2_wdb_path.exists():
        stage2_wdb = pd.read_csv(stage2_wdb_path)
        stage2_wdb["genome"] = stage2_wdb["genome"].apply(strip_ext)
        win_col2 = "cluster" if "cluster" in stage2_wdb.columns else s2_col
        stage2_winner_map = dict(zip(stage2_wdb[win_col2], stage2_wdb["genome"]))
    else:
        stage2_winner_map = None

    stage2_cluster_map = dict(zip(stage2_df["genome"], stage2_df[s2_col]))

    def resolve_final_genome(stage1_winner):
        if stage1_winner is None or stage1_winner not in stage2_cluster_map:
            return None
        s2_cluster = stage2_cluster_map[stage1_winner]
        if stage2_winner_map is not None:
            return stage2_winner_map.get(s2_cluster)
        return stage1_winner

    stage1_df["final_genome"] = stage1_df["stage1_winner"].apply(resolve_final_genome)
    missing = stage1_df["final_genome"].isna().sum()
    if missing:
        print(f"WARNING: {missing} genomes could not be traced to a final genome "
              f"— check Stage1/Stage2 table naming assumptions")
    return stage1_df.dropna(subset=["final_genome"])[["genome", "final_genome"]]

def build_gene_function_map() -> pd.DataFrame:
    """final_genome -> {KOs, CAZy families, VF genes} via cluster membership."""
    print("Loading gene_to_genome.tsv (this is large, ~47M rows)...")
    g2g = pd.read_csv(GENE2GENOME, sep="\t", header=None, names=["gene_id", "genome"],
                       dtype=str)
    g2g["genome"] = g2g["genome"].apply(strip_ext)

    g2g["gene_id_prefixed"] = g2g["genome"] + "~" + g2g["gene_id"]

    print("Loading MMseqs2 cluster membership table...")
    clu = pd.read_csv(CLUSTER_TABLE, sep="\t", header=None,
                       names=["cluster_rep", "member_gene"], dtype=str)

    merged = g2g.merge(clu, left_on="gene_id_prefixed", right_on="member_gene", how="left")
    match_rate = merged["cluster_rep"].notna().mean()
    print(f"gene_to_genome <-> cluster_table join match rate: {match_rate:.1%}")
    if match_rate < 0.95:
        print("WARNING: match rate below 95% even with genome~header key — "
              "inspect a few unmatched rows manually before trusting downstream "
              "annotation joins. Possible cause: Prodigal metadata (the ' # ... '"
              "suffix) wasn't stripped consistently somewhere in the pipeline.")

    return merged[["genome", "cluster_rep"]].dropna()

def load_annotation_tables():
    ann = {}

    if GTDBTK_SUMMARY.exists():
        gtdb = pd.read_csv(GTDBTK_SUMMARY, sep="\t", dtype=str)
        gcol = next(c for c in gtdb.columns if "genome" in c.lower())
        gtdb["genome"] = gtdb[gcol].apply(strip_ext)
        ann["taxonomy"] = gtdb.set_index("genome")
    else:
        print(f"NOTE: {GTDBTK_SUMMARY} not found — taxonomy will be blank. "
              f"Re-run this script once GTDB-Tk completes.")
        ann["taxonomy"] = None

    if EGGNOG_ANNOT.exists():
        header_line_idx = None
        with open(EGGNOG_ANNOT) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#query"):
                    header_line_idx = i
                    break
                if i > 20:
                    break

        if header_line_idx is None:
            print("WARNING: could not locate the '#query' header line in "
                  "eggNOG output within the first 20 lines — check the file "
                  "manually. Skipping eggNOG annotation.")
            ann["eggnog"] = None
        else:
            egg = pd.read_csv(EGGNOG_ANNOT, sep="\t", dtype=str,
                               skiprows=header_line_idx, header=0)
            egg = egg[~egg.iloc[:, 0].astype(str).str.startswith("#")]
            egg.columns = [str(c).lstrip("#") for c in egg.columns]
            qcol = egg.columns[0]
            ko_col = next((c for c in egg.columns if "KEGG_ko" in c), None)
            ann["eggnog"] = egg.set_index(qcol)[[ko_col]] if ko_col else None
            if ko_col is None:
                print("WARNING: could not find a KEGG_ko column in eggNOG "
                      f"output. Actual columns: {list(egg.columns)}")
    else:
        ann["eggnog"] = None

    if DBCAN_OVERVIEW.exists():
        dbcan_raw = pd.read_csv(DBCAN_OVERVIEW, sep="\t", dtype=str)
        ann["dbcan"] = dbcan_raw.set_index(dbcan_raw.columns[0])
    else:
        print(f"NOTE: {DBCAN_OVERVIEW} not found — check the real dbCAN output "
              f"filename/path and update DBCAN_OVERVIEW above.")
        ann["dbcan"] = None

    if VFDB_HITS.exists():
        vfdb = pd.read_csv(VFDB_HITS, sep="\t", header=None, dtype=str)
        vfdb.columns = ["qseqid"] + [f"col{i}" for i in range(1, vfdb.shape[1])]
        ann["vfdb"] = vfdb.groupby("qseqid").size().rename("n_vf_hits")
    else:
        ann["vfdb"] = None

    if AMR_MERGED.exists():
        amr = pd.read_csv(AMR_MERGED, sep="\t", dtype=str)
        amr["genome"] = amr["genome"].apply(strip_ext)
        ann["amr"] = amr
    else:
        ann["amr"] = None

    return ann

def main():
    print("=== Step 1: host attribution via Stage1->Stage2 chaining ===")
    attribution = build_host_attribution()

    meta = pd.read_csv(METADATA_FILE, sep="\t", dtype=str)
    meta["genome"] = meta["genome"].apply(strip_ext)
    attribution = attribution.merge(meta[["genome", "host"]], on="genome", how="left")

    host_species = (attribution.groupby(["host", "final_genome"])
                     .size().rename("n_contributing_genomes").reset_index())
    print(f"{host_species['host'].nunique()} hosts, "
          f"{host_species['final_genome'].nunique()} distinct final species referenced")

    print("\n=== Step 2: gene -> genome -> cluster_rep map ===")
    gene_func_map = build_gene_function_map()

    print("\n=== Step 3: loading annotation tables ===")
    ann = load_annotation_tables()

    print("\n=== Step 4: assembling per-genome functional summaries ===")
    final_genomes = host_species["final_genome"].unique()

    print("Pre-grouping gene_func_map by genome (one-time cost, a few minutes)...")
    reps_by_genome = gene_func_map.groupby("genome")["cluster_rep"].apply(
        lambda s: s.unique()
    ).to_dict()
    print(f"Grouped {len(reps_by_genome)} genomes.")

    genome_profiles = []
    for i, fg in enumerate(final_genomes):
        if i % 2000 == 0:
            print(f"  {i}/{len(final_genomes)} genomes processed...")
        profile = {"final_genome": fg}

        if ann["taxonomy"] is not None and fg in ann["taxonomy"].index:
            row = ann["taxonomy"].loc[fg]
            profile["classification"] = row.get("classification", None)

        my_reps = reps_by_genome.get(fg, [])

        if ann["eggnog"] is not None:
            kos = ann["eggnog"].reindex(my_reps).iloc[:, 0].dropna()
            kos = kos[kos != "-"]
            profile["KOs"] = ";".join(sorted(set(
                k for cell in kos for k in str(cell).split(",")
            )))

        if ann["dbcan"] is not None:
            fams = ann["dbcan"].reindex(my_reps)
            import re
            FAMILY_COLS = [c for c in ["dbCAN_hmm", "dbCAN_sub", "DIAMOND", "Recommend Results"]
                           if c in fams.columns]
            all_fams = set()
            for col in FAMILY_COLS:
                for cell in fams[col].dropna():
                    cell = str(cell)
                    if cell == "-" or not cell:
                        continue
                    for piece in cell.split("+"):
                        fam = re.sub(r"\(.*?\)", "", piece).strip()
                        if fam and fam != "-":
                            all_fams.add(fam)
            profile["CAZy_families"] = ";".join(sorted(all_fams))

            if "Substrate" in fams.columns:
                subs = set(s for s in fams["Substrate"].dropna() if s not in ("-", ""))
                profile["CAZy_substrates"] = ";".join(sorted(subs))

        if ann["vfdb"] is not None:
            profile["n_VF_genes"] = int(ann["vfdb"].reindex(my_reps).fillna(0).sum())

        if ann["amr"] is not None:
            amr_rows = ann["amr"][ann["amr"]["genome"] == fg]
            symcol = next((c for c in amr_rows.columns if "symbol" in c.lower()), None)
            profile["AMR_genes"] = ";".join(sorted(amr_rows[symcol].dropna().unique())) if symcol else ""

        genome_profiles.append(profile)

    profiles_df = pd.DataFrame(genome_profiles)
    final_table = host_species.merge(profiles_df, on="final_genome", how="left")

    out_tsv = OUT_DIR / "per_host_species_profile.tsv"
    final_table.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote {len(final_table)} rows to {out_tsv}")

    for host, sub in final_table.groupby("host"):
        safe_name = str(host).replace(" ", "_").replace("/", "_")
        (OUT_DIR / f"{safe_name}.json").write_text(
            sub.drop(columns=["host"]).to_json(orient="records", indent=2)
        )
    print(f"Wrote one JSON file per host to {OUT_DIR}")

if __name__ == "__main__":
    main()
