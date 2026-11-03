#!/usr/bin/env python3
"""
build_web_tables.py — build the website data tables for the URGG catalog portal
(see URGG_website_build_spec.md). Reads existing pipeline outputs only; writes into OUT.

Usage (on the server):
    python3 build_web_tables.py                 # everything
    python3 build_web_tables.py --skip-fasta    # tables only (no gzip of 25,451 genomes)
    python3 build_web_tables.py --skip-ko       # skip the KO index (streams the 8.5 GB cluster file)
    python3 build_web_tables.py --threads 64    # default 64

Parallel layout: the FASTA step (25,451 genomes -> pigz) runs in the background on
--threads cores (worker processes x pigz threads) while the tables and the KO index
(streaming the 8.5 GB cluster file) are built at the same time in the main process.
Uses pigz if found on PATH (falls back to Python gzip, much slower).

Requires: python3 with pandas (the drep conda env is fine).
The script STOPS (exit 1) if any headline total differs from the paper.
"""
import argparse, collections, csv, gzip, json, os, re, shutil, subprocess, sys, threading, time
from multiprocessing import Pool
import pandas as pd

# ------------------------------------------------------------------ paths (edit if moved)
B   = "/disk252/qamar/metagenome_analysis"
WR  = f"{B}/website_release"
P = dict(
    tax_qc      = f"{WR}/urgg_tax_QC_table.tsv",
    id_map      = f"{WR}/genome_id_mapping.tsv",
    genomes_dir = f"{WR}/genomes",
    metadata    = f"{B}/drep_analysis/metadata.tsv",
    g2s         = f"{B}/drep_analysis/genome_to_species.tsv",
    host_range  = f"{B}/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv",
    strain_map  = f"{B}/drep_analysis/strain_id_map.tsv",
    thiamine    = f"{B}/annotation/results_phase3/s63_thiamine_guild_propagated/63_genome_thiamine_guilds.tsv",
    cazy        = f"{B}/annotation/results_phase3/s54_cazyme_propagated/54_cazy_family_presence_propagated.tsv",
    subpath50   = f"{B}/annotation/results_phase3/s50_subpathway_propagated/subpathway_presence_propagated_50.tsv",
    coop        = f"{B}/annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_category_presence_per_genome_WIDE.csv",
    vf          = f"{B}/annotation/results_phase3/s64_vfdb_propagated/64_vf_presence_propagated.tsv",
    amr         = f"{B}/annotation/amr_results_with_hits/all_mags_amr_combined.tsv",
    reads       = f"{B}/annotation/results_phase3/s79_readbased_host_range/79_TableS_species_host_range_assembly_vs_reads.csv",
    emapper     = f"{B}/annotation/eggnog_output/nr_catalog.emapper.annotations",
    clusters    = f"{B}/annotation/gene_catalog/nr_gene_catalog_cluster.tsv",
)
OUT = f"{WR}/web_data"

# ------------------------------------------------------------------ paper values (hard checks)
EXPECT = dict(species=25451, qc_genomes=117145, strains=23999, near_complete=10546, novel=16548,
              described=8903, novel_species=16174, novel_genus=367, novel_family=6, novel_order=1,
              host_specific=20186, shared=4811, core=454, singletons=12571, ruc=76)
QC_BY_HOST = {"Ovis aries":37202,"Bos taurus":23730,"Bos grunniens":20470,"Pantholops hodgsonii":20026,
    "Capra hircus":6849,"Bubalus bubalis":4449,"Bos indicus":1979,"Hydropotes inermis":1254,"Alces alces":403,
    "Procapra picticaudata":179,"Cervus elaphus":175,"Giraffa camelopardalis":159,"Rangifer tarandus":125,
    "Alces americanus":82,"Camelus bactrianus":49,"Camelus dromedarius":10,"Bos frontalis":4}
SPECIES_BY_HOST = {"Ovis aries":6476,"Pantholops hodgsonii":4214,"Bos taurus":3755,"Bos grunniens":3605,
    "Capra hircus":2778,"Bubalus bubalis":2456,"Bos indicus":950,"Hydropotes inermis":374,"Alces alces":296,
    "Procapra picticaudata":150,"Cervus elaphus":125,"Rangifer tarandus":100,"Alces americanus":80,
    "Giraffa camelopardalis":52,"Camelus bactrianus":34,"Camelus dromedarius":3,"Bos frontalis":3}
COMMON = {"Ovis aries":"sheep","Pantholops hodgsonii":"Tibetan antelope","Bos taurus":"cattle","Bos grunniens":"yak",
    "Capra hircus":"goat","Bubalus bubalis":"buffalo","Bos indicus":"zebu cattle","Hydropotes inermis":"water deer",
    "Alces alces":"moose","Procapra picticaudata":"Tibetan gazelle","Cervus elaphus":"red deer","Rangifer tarandus":"reindeer",
    "Alces americanus":"North American moose","Giraffa camelopardalis":"giraffe","Camelus bactrianus":"Bactrian camel",
    "Camelus dromedarius":"dromedary","Bos frontalis":"gaur"}
FAMILY = {h:"Bovidae" for h in ["Ovis aries","Pantholops hodgsonii","Bos taurus","Bos grunniens","Capra hircus",
    "Bubalus bubalis","Bos indicus","Procapra picticaudata","Bos frontalis"]}
FAMILY.update({h:"Cervidae" for h in ["Hydropotes inermis","Alces alces","Cervus elaphus","Rangifer tarandus","Alces americanus"]})
FAMILY.update({"Giraffa camelopardalis":"Giraffidae","Camelus bactrianus":"Camelidae","Camelus dromedarius":"Camelidae"})
READ_BIOPROJECTS = set("""PRJNA60251 PRJEB10338 PRJEB8939 PRJEB12797 PRJNA301235 PRJNA362629 PRJNA512570 PRJNA214227
PRJNA485488 PRJNA486341 PRJEB14312 PRJNA419351 PRJNA448333 PRJNA657455 PRJNA526070 PRJNA298982 PRJCA008889 PRJNA972320
PRJNA983076 PRJEB38078 PRJNA260105 PRJNA603266 PRJNA723218 PRJNA1069655 PRJNA1034944""".split())

FAIL = []
def check(name, got, want):
    ok = got == want
    print(f"  [{'OK' if ok else 'FAIL'}] {name}: {got:,} (paper {want:,})")
    if not ok: FAIL.append(name)
def log(m): print(time.strftime("%H:%M:%S"), m, flush=True)
def host_name(x): return str(x).replace("_", " ").strip()
def tier_label(x):
    x = str(x).lower()
    return "Core" if x.startswith(("core", "broad")) else "Shared" if x.startswith("shared") else "Host-specific" if x.startswith("host") else ""

def novelty(cls):
    r = dict(p.split("__", 1) for p in str(cls).split(";") if "__" in p)
    if not r.get("o"): return "Novel order"
    if not r.get("f"): return "Novel family"
    if not r.get("g"): return "Novel genus"
    if not r.get("s"): return "Novel species"
    return "Described"
def rank(cls, k):
    r = dict(p.split("__", 1) for p in str(cls).split(";") if "__" in p)
    return r.get(k, "")

def read_any(path):
    sep = "," if path.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)

def per_genome_lists(path, orig2urgg, label):
    """Return {URGG: [features]} from a wide (genome x feature 0/1) or long (genome, feature[, present]) table."""
    if not os.path.exists(path):
        log(f"  WARNING {label}: {path} not found, skipped"); return {}
    df = read_any(path)
    gcol = max(df.columns, key=lambda c: df[c].astype(str).isin(orig2urgg.keys()).mean() if not pd.api.types.is_numeric_dtype(df[c]) else 0)
    hit = df[gcol].astype(str).isin(orig2urgg.keys()).mean()
    if hit < 0.5:
        log(f"  WARNING {label}: no genome column matched URGG originals (best {gcol}, {hit:.0%}); skipped"); return {}
    other = [c for c in df.columns if c != gcol]
    def is_binary(c):
        v = set(df[c].dropna().astype(str).str.upper().unique())
        return len(v) > 0 and v <= {"0", "1", "0.0", "1.0", "TRUE", "FALSE"}
    num = [c for c in other if is_binary(c)]
    out = collections.defaultdict(list)
    if len(num) >= 2 and len(num) >= len(other) - 3:                     # wide 0/1
        vals = df[num].replace({"True":1,"False":0,"TRUE":1,"FALSE":0}).apply(pd.to_numeric, errors="coerce").fillna(0)
        for g, row in zip(df[gcol].astype(str), vals.values):
            u = orig2urgg.get(g)
            if u: out[u] = [c for c, v in zip(num, row) if v > 0]
        log(f"  {label}: wide table, {len(num)} features, {len(out):,} genomes")
    else:                                                                 # long
        objs = [c for c in other if not pd.api.types.is_numeric_dtype(df[c]) and c not in num]
        fcol = max(objs, key=lambda c: df[c].nunique())          # the feature column has the most distinct values
        pres = [c for c in other if c in num]
        sub = df if not pres else df[pd.to_numeric(df[pres[0]].replace({"True":1,"False":0,"TRUE":1,"FALSE":0}), errors="coerce").fillna(0) > 0]
        for g, f in zip(sub[gcol].astype(str), sub[fcol].astype(str)):
            u = orig2urgg.get(g)
            if u: out[u].append(f)
        out = {k: sorted(set(v)) for k, v in out.items()}
        log(f"  {label}: long table (feature column '{fcol}'), {len(out):,} genomes")
    return out

def vfdb_by_level(path, orig2urgg):
    """64_vf_presence_propagated.tsv is long: genome, level, trait. Split by level -> {level: {urgg: [traits]}}."""
    if not os.path.exists(path): log(f"  WARNING VFDB: {path} not found, skipped"); return {}
    v = pd.read_csv(path, sep="\t", dtype=str)
    if "level" not in v.columns:
        return {"homologs": per_genome_lists(path, orig2urgg, "VFDB homologs")}
    v["u"] = v.genome.map(orig2urgg); v = v[v.u.notna()]
    out = {}
    for lv, sub in v.groupby("level"):
        key = re.sub(r"[^A-Za-z0-9]+", "_", str(lv)).strip("_").lower()
        out[key] = sub.groupby("u").trait.apply(lambda s: sorted(set(s))).to_dict()
        log(f"  VFDB level '{lv}': {sub.trait.nunique():,} distinct, {len(out[key]):,} genomes")
    return out

# ------------------------------------------------------------------ FASTA: gzip + header fix + stats
PIGZ = shutil.which("pigz")

def gz_bytes(data, dst, threads):
    """Compress bytes to dst with pigz (multi-threaded) or Python gzip as fallback."""
    if PIGZ:
        with open(dst, "wb") as fo:
            subprocess.run([PIGZ, "-6", "-p", str(threads), "-c"], input=data, stdout=fo, check=True)
    else:
        with gzip.open(dst, "wb", compresslevel=6) as fo: fo.write(data)

def gz_write_text(dst, lines, threads=8):
    gz_bytes("".join(lines).encode(), dst, threads)

def fasta_job(args):
    src, dst, pz = args
    with open(src, "rb") as fi: lines = fi.read().split(b"\n")
    contigs, seq = 0, []
    for i, l in enumerate(lines):
        if l.startswith(b">"):
            contigs += 1; lines[i] = l.replace(b"|", b"_")
        else:
            seq.append(l.strip())
    s = b"".join(seq).upper()
    gc = s.count(b"G") + s.count(b"C"); at = s.count(b"A") + s.count(b"T")
    gz_bytes(b"\n".join(lines), dst, pz)
    return os.path.basename(src).split(".")[0], len(s), contigs, round(100 * gc / max(1, gc + at), 2)

def run_fasta(jobs, workers, stats, err):
    """Background thread: drives a process pool; each worker pipes one genome through pigz."""
    try:
        with Pool(workers) as pool:
            for i, (u, n, c, g) in enumerate(pool.imap_unordered(fasta_job, jobs, chunksize=20), 1):
                stats[u] = (n, c, g)
                if i % 2500 == 0: log(f"  [FASTA] {i:,} / {len(jobs):,}")
        log(f"  [FASTA] done: {len(stats):,} genomes gzipped")
    except Exception as e:
        err.append(repr(e))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fasta", action="store_true")
    ap.add_argument("--skip-ko", action="store_true")
    ap.add_argument("--threads", type=int, default=64, help="threads for the FASTA/pigz step (default 64)")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    # ---------- core tables
    log("Reading genome tables")
    tq = pd.read_csv(P["tax_qc"], sep="\t")
    o2u = dict(zip(tq.original_genome_id, tq.public_id))
    meta = pd.read_csv(P["metadata"], sep="\t", low_memory=False)
    g2s = pd.read_csv(P["g2s"], sep="\t", low_memory=False)
    hr = pd.read_csv(P["host_range"], low_memory=False)
    sm = pd.read_csv(P["strain_map"], sep="\t", low_memory=False)

    print("Checks:")
    check("species-level genomes", len(tq), EXPECT["species"])
    check("QC-passing genomes (metadata)", len(meta), EXPECT["qc_genomes"])
    check("genome->species map", len(g2s), EXPECT["qc_genomes"])
    check("strain-level genomes", len(sm), EXPECT["strains"])
    check("near-complete, all QC genomes", int(((meta.completeness > 90) & (meta.contamination <= 5)).sum()), 38408)
    check("near-complete species-level", int(((tq.completeness > 90) & (tq.contamination <= 5)).sum()), EXPECT["near_complete"])

    tq["novelty"] = tq.gtdb_classification.map(novelty)
    nv = tq.novelty.value_counts()
    check("described", int(nv.get("Described", 0)), EXPECT["described"])
    check("novel species", int(nv.get("Novel species", 0)), EXPECT["novel_species"])
    check("novel genus", int(nv.get("Novel genus", 0)), EXPECT["novel_genus"])
    check("novel family", int(nv.get("Novel family", 0)), EXPECT["novel_family"])
    check("novel order", int(nv.get("Novel order", 0)), EXPECT["novel_order"])

    hr["urgg_id"] = hr.representative.map(o2u)
    check("host-range rows matched to URGG IDs", int(hr.urgg_id.notna().sum()), EXPECT["species"])
    hr["tier"] = hr.sharing_category.map(tier_label)
    tc = hr.tier.value_counts()
    check("host-specific", int(tc.get("Host-specific", 0)), EXPECT["host_specific"])
    check("shared", int(tc.get("Shared", 0)), EXPECT["shared"])
    check("core", int(tc.get("Core", 0)), EXPECT["core"])
    check("singletons", int((hr.n_member_genomes == 1).sum()), EXPECT["singletons"])

    meta["host_species"] = meta.species.map(host_name)
    qc_host = meta.host_species.value_counts().to_dict()
    for h, v in QC_BY_HOST.items(): check(f"QC genomes, {h}", int(qc_host.get(h, 0)), v)
    tq["rep_host"] = tq.species.map(host_name)
    sp_host = tq.rep_host.value_counts().to_dict()
    for h, v in SPECIES_BY_HOST.items(): check(f"species-level, {h}", int(sp_host.get(h, 0)), v)
    tq["in_ruc"] = tq.gtdb_classification.str.contains(";o__UMGS1840;", regex=False)
    check("RUC lineage (order UMGS1840)", int(tq.in_ruc.sum()), EXPECT["ruc"])
    if FAIL:
        print(f"\nSTOP: {len(FAIL)} check(s) failed: {FAIL}. Nothing written."); sys.exit(1)

    # ---------- cluster -> URGG, members
    cl2u = dict(zip(hr.cluster.astype(str), hr.urgg_id))
    g2s["urgg_id"] = g2s.species_rep.astype(str).map(cl2u)
    if g2s.urgg_id.isna().any():
        # species_rep may hold the representative genome name instead of a cluster id
        g2s.loc[g2s.urgg_id.isna(), "urgg_id"] = g2s.loc[g2s.urgg_id.isna(), "species_rep"].map(o2u)
    log(f"members mapped to a species: {g2s.urgg_id.notna().sum():,} / {len(g2s):,}")
    if g2s.urgg_id.isna().any():
        print("STOP: some genomes could not be mapped to a species URGG ID"); sys.exit(1)

    s_orig2id = dict(zip(sm.original_genome_id, sm.strain_id))
    mem = g2s[["genome", "urgg_id"]].merge(meta[["genome", "bioproject", "host_species", "sample_site",
                                                   "completeness", "contamination", "N50"]], on="genome", how="left")
    mem["strain_id"] = mem.genome.map(s_orig2id).fillna("")
    mem["is_representative"] = mem.genome.isin(o2u.keys())
    mem["origin"] = mem.bioproject.astype(str).map(lambda b: "Reprocessed reads" if b in READ_BIOPROJECTS else "Public MAG")
    mem["quality"] = ((mem.completeness > 90) & (mem.contamination <= 5)).map({True: "Near-complete", False: "Other"})
    mem = mem.rename(columns={"genome": "original_genome_id", "sample_site": "site"})
    mem.to_csv(f"{OUT}/web_members.tsv", sep="\t", index=False)
    log(f"web_members.tsv: {len(mem):,} rows")

    # ---------- strains
    st = sm.merge(meta[["genome", "host_species", "sample_site", "completeness", "contamination"]],
                  left_on="original_genome_id", right_on="genome", how="left").drop(columns="genome")
    st["species_urgg_id"] = st.original_genome_id.map(dict(zip(mem.original_genome_id, mem.urgg_id)))
    st.to_csv(f"{OUT}/web_strains.tsv", sep="\t", index=False)
    strains_by_sp = st.groupby("species_urgg_id").strain_id.apply(lambda s: ";".join(sorted(s))).to_dict()
    log(f"web_strains.tsv: {len(st):,} rows")

    # ---------- FASTA
    # Started in the background; the tables + KO index below run at the same time.
    stats, fasta_err, fasta_thread = {}, [], None
    if not a.skip_fasta:
        fdir = f"{OUT}/fasta"; os.makedirs(fdir, exist_ok=True)
        pz = 4 if PIGZ else 1                              # pigz threads per genome
        workers = max(1, a.threads // pz)                  # e.g. 64 threads = 16 workers x pigz -p 4
        jobs = [(f"{P['genomes_dir']}/{u}.fa", f"{fdir}/{u}.fa.gz", pz) for u in tq.public_id]
        log(f"[FASTA] started in background: {len(jobs):,} genomes, {workers} workers x "
            + (f"pigz -p {pz}" if PIGZ else "python gzip (pigz not found!)"))
        fasta_thread = threading.Thread(target=run_fasta, args=(jobs, workers, stats, fasta_err), daemon=False)
        fasta_thread.start()

    # ---------- genomes
    g = tq.rename(columns={"public_id": "urgg_id", "host": "rep_host_common_old", "species": "rep_host_species_raw",
                           "sample_site": "rep_site", "N50": "n50"})
    g["rep_host_species"] = g.rep_host_species_raw.map(host_name)
    g["rep_host_common"] = g.rep_host_species.map(COMMON)
    g["rep_host_family"] = g.rep_host_species.map(FAMILY)
    g["quality"] = ((g.completeness > 90) & (g.contamination <= 5)).map({True: "Near-complete", False: "Other"})
    for k, col in [("d", "domain"), ("p", "phylum"), ("c", "class"), ("o", "order"), ("f", "family"), ("g", "genus"), ("s", "gtdb_species")]:
        g[col] = g.gtdb_classification.map(lambda x, k=k: rank(x, k))
    h = hr.set_index("urgg_id")
    g["hosts_n"] = g.urgg_id.map(h.n_host_species)
    g["hosts_list"] = g.urgg_id.map(h.host_species_list)
    g["sharing_tier"] = g.urgg_id.map(h.tier)
    g["cluster_size"] = g.urgg_id.map(h.n_member_genomes)
    g["origin"] = g.bioproject.astype(str).map(lambda b: "Reprocessed reads" if b in READ_BIOPROJECTS else "Public MAG")
    g["strain_ids"] = g.urgg_id.map(strains_by_sp).fillna("")
    def write_genomes():
      if stats:
        g["genome_size_bp"] = g.urgg_id.map(lambda u: stats.get(u, ("", "", ""))[0])
        g["contigs"] = g.urgg_id.map(lambda u: stats.get(u, ("", "", ""))[1])
        g["gc_pct"] = g.urgg_id.map(lambda u: stats.get(u, ("", "", ""))[2])
      g["in_ruc_lineage"] = g.in_ruc
      g["fasta"] = "fasta/" + g.urgg_id + ".fa.gz"
      keep = ["urgg_id", "original_genome_id", "origin", "bioproject", "rep_host_species", "rep_host_common", "rep_host_family",
              "rep_site", "hosts_n", "hosts_list", "sharing_tier", "cluster_size", "completeness", "contamination", "quality",
              "n50"] + (["genome_size_bp", "contigs", "gc_pct"] if stats else []) + \
             ["domain", "phylum", "class", "order", "family", "genus", "gtdb_species", "novelty", "gtdb_classification",
              "strain_ids", "in_ruc_lineage", "fasta"]
      g[keep].to_csv(f"{OUT}/web_genomes.tsv", sep="\t", index=False)
      g[keep].drop(columns="fasta").to_csv(f"{WR}/URGG_genome_metadata.tsv", sep="\t", index=False)   # OMIX upload copy
      log(f"web_genomes.tsv: {len(g):,} rows (also URGG_genome_metadata.tsv for OMIX)")

    # ---------- per-genome function
    log("Per-genome function tables")
    thi = pd.read_csv(P["thiamine"], sep="\t", low_memory=False)
    thi_g = dict(zip(thi.genome.map(o2u), thi.guild))
    cazy = per_genome_lists(P["cazy"], o2u, "CAZy families")
    sub = per_genome_lists(P["subpath50"], o2u, "metabolic sub-pathways (50%)")
    coop = per_genome_lists(P["coop"], o2u, "cooperative categories")
    vf_levels = vfdb_by_level(P["vf"], o2u)        # one column per level (vf_gene, factor, category ...)
    amr = collections.defaultdict(set)
    with open(P["amr"]) as f:                       # AMRFinderPlus rows, genome name prepended; Type is column 10
        for line in f:
            x = line.rstrip("\n").split("\t")
            if len(x) > 9 and x[9] == "AMR" and x[0] in o2u: amr[o2u[x[0]]].add(x[6])
    gf = pd.DataFrame({"urgg_id": g.urgg_id})
    gf["thiamine_guild"] = gf.urgg_id.map(thi_g).fillna("")
    for name, d in [("cazy_families", cazy), ("subpathways_50pct", sub), ("cooperative_categories", coop)] + \
                   [(f"vfdb_{lv}", d_) for lv, d_ in vf_levels.items()]:
        gf[name] = gf.urgg_id.map(lambda u, d=d: ";".join(d.get(u, [])))
    gf["amr_genes"] = gf.urgg_id.map(lambda u: ";".join(sorted(amr.get(u, []))))
    gf.to_csv(f"{OUT}/web_genome_function.tsv", sep="\t", index=False)
    log(f"web_genome_function.tsv: AMR-positive {sum(1 for v in amr.values() if v):,} (paper 294)")

    # ---------- family index helper (counts by host and tier, top genera)
    hosts_of = {}
    sep_guess = lambda s: ";" if ";" in s else ","
    for u, l in zip(h.index, h.host_species_list):
        hosts_of[u] = [x.strip() for x in str(l).split(sep_guess(str(l))) if x.strip()]
    tier_of = dict(zip(g.urgg_id, g.sharing_tier)); genus_of = dict(zip(g.urgg_id, g.genus))
    tier_n = collections.Counter(tier_of.values())
    def index_rows(fam2genomes, names=None):
        rows, lists = [], []
        for fam, gs in sorted(fam2genomes.items()):
            gs = sorted(gs); hc = collections.Counter(hh for u in gs for hh in hosts_of.get(u, []))
            tcount = collections.Counter(tier_of.get(u) for u in gs)
            gen = collections.Counter(genus_of.get(u) or "unassigned" for u in gs).most_common(10)
            rows.append(dict(family=fam, name=(names or {}).get(fam, ""), genomes_n=len(gs),
                             genomes_pct=round(100 * len(gs) / EXPECT["species"], 2),
                             **{f"core_pct": round(100 * tcount["Core"] / max(tier_n["Core"], 1), 1),
                                f"shared_pct": round(100 * tcount["Shared"] / max(tier_n["Shared"], 1), 1),
                                f"hostspecific_pct": round(100 * tcount["Host-specific"] / max(tier_n["Host-specific"], 1), 1)},
                             by_host=";".join(f"{k}:{v}" for k, v in hc.most_common()),
                             top_genera=";".join(f"{k}:{v}" for k, v in gen)))
            lists.append((fam, ";".join(gs)))
        return rows, lists

    fam2g = collections.defaultdict(set)
    for u, fams in cazy.items():
        for fa in fams: fam2g[fa].add(u)
    rows, lists = index_rows(fam2g)
    pd.DataFrame(rows).to_csv(f"{OUT}/web_cazy_index.tsv", sep="\t", index=False)
    gz_write_text(f"{OUT}/web_cazy_genomes.tsv.gz", ["family\turgg_ids\n"] + [f"{a_}\t{b_}\n" for a_, b_ in lists])
    log(f"web_cazy_index.tsv: {len(rows)} families")

    # ---------- KO index (propagated: representative-gene KO -> all cluster members -> genome)
    if not a.skip_ko:
        log("Reading eggNOG KO calls of representative genes")
        rep_ko, ko_id, ko_name = {}, {}, collections.defaultdict(collections.Counter)
        with open(P["emapper"]) as f:
            hdr = None
            for line in f:
                if line.startswith("##"): continue
                x = line.rstrip("\n").split("\t")
                if hdr is None: hdr = {k.lstrip("#"): i for i, k in enumerate(x)}; continue
                kos = x[hdr["KEGG_ko"]]
                if kos in ("-", ""): continue
                ids = []
                for k in kos.split(","):
                    k = k.replace("ko:", "")
                    if k not in ko_id: ko_id[k] = len(ko_id)
                    ids.append(ko_id[k])
                    pn = x[hdr["Preferred_name"]]
                    if pn not in ("-", ""): ko_name[k][pn] += 1
                rep_ko[x[hdr["query"]]] = tuple(ids)
        log(f"  {len(rep_ko):,} representative genes with a KO; {len(ko_id):,} KOs")
        uidx = {u: i for i, u in enumerate(g.urgg_id)}; nG = len(uidx)
        bits = [bytearray((nG + 7) // 8) for _ in range(len(ko_id))]
        log("Streaming gene-catalog clusters (this takes a while)")
        with open(P["clusters"]) as f:
            for i, line in enumerate(f, 1):
                rep, member = line.rstrip("\n").split("\t", 1)
                ids = rep_ko.get(rep)
                if not ids: continue
                u = o2u.get(member.split("~", 1)[0])
                if u is None: continue
                j = uidx[u]
                for k in ids: bits[k][j >> 3] |= 1 << (j & 7)
                if i % 10_000_000 == 0: log(f"  {i:,} lines")
        import numpy as np
        ulist = np.array(g.urgg_id)
        fam2g = {}
        for k, i in ko_id.items():
            on = np.unpackbits(np.frombuffer(bytes(bits[i]), dtype=np.uint8), bitorder="little")[:nG].nonzero()[0]
            if len(on): fam2g[k] = set(ulist[on])
        names = {k: c.most_common(1)[0][0] for k, c in ko_name.items()}
        rows, lists = index_rows(fam2g, names)
        pd.DataFrame(rows).rename(columns={"family": "ko"}).to_csv(f"{OUT}/web_ko_index.tsv", sep="\t", index=False)
        gz_write_text(f"{OUT}/web_ko_genomes.tsv.gz", ["ko\turgg_ids\n"] + [f"{a_}\t{b_}\n" for a_, b_ in lists])
        thic = len(fam2g.get("K03147", ()))
        log(f"web_ko_index.tsv: {len(rows):,} KOs; thiC (K03147) carriers {thic:,} (paper, Suppl. Table 10: 3,671)")

    # ---------- read-based detection (Suppl. Table 3 rule), original names -> URGG IDs
    if os.path.exists(P["reads"]):
        rd = read_any(P["reads"])
        for c in rd.columns:
            if not pd.api.types.is_numeric_dtype(rd[c]) and rd[c].astype(str).isin(o2u.keys()).mean() > 0.5:
                rd.insert(0, "urgg_id", rd[c].map(o2u)); break
        rd.to_csv(f"{OUT}/web_read_detection.tsv", sep="\t", index=False)
        log(f"web_read_detection.tsv: {len(rd):,} rows")

    # ---------- overview JSON
    site_counts = meta.sample_site.value_counts().to_dict()
    phy = g.groupby("phylum").agg(n=("urgg_id", "size"), novel=("novelty", lambda s: (s != "Described").sum()))
    ov = dict(
        at_a_glance={"studies": 39, "bioprojects": 54, "reprocessed_bioprojects": 25, "public_mag_studies": 14,
                     "public_mag_bioprojects": 29, "quality_filtered_metagenomes": 1447, "host_species": 17,
                     "ruminant": 15, "pseudo_ruminant": 2, "sampling_sites": 13, "bins_from_assembled_samples": 207159,
                     "qc_genomes": len(meta), "near_complete_all": 38408, "species_level": len(g),
                     "near_complete_species": int((g.quality == "Near-complete").sum()), "novel_species_level": int((g.novelty != "Described").sum()),
                     "strain_level": len(st), "median_completeness": 87.5, "median_contamination": 1.57,
                     "predicted_genes": 46918358, "nonredundant_genes": 37129697, "votus": 8980, "votus_viral": 6619},
        hosts=[{"host": hh, "common": COMMON[hh], "family": FAMILY[hh], "qc_genomes": int(qc_host.get(hh, 0)),
                "species_level": int(sp_host.get(hh, 0))} for hh in QC_BY_HOST],
        sites_qc_genomes=site_counts,
        novelty=g.novelty.value_counts().to_dict(),
        tiers={k: int(v) for k, v in tc.items()}, singletons=int((hr.n_member_genomes == 1).sum()),
        phyla={p_: {"n": int(r.n), "novel_pct": round(100 * r.novel / r.n, 1)} for p_, r in phy.sort_values("n", ascending=False).iterrows()},
        reference_matches_pct={"RGMC": 48.7, "GTDB R232": 34.7, "RGMGC": 13.5, "MGnify cow rumen": 10.4, "MGnify sheep rumen": 5.1,
                               "Human gut": 1.9, "Pig gut": 1.5, "Mouse gut": 0.6, "Hungate1000": 0.4, "Chicken gut": 0.3,
                               "No species-level match": 41.1},
        virome={"sequences": 15301, "votus": 8980, "viral": 6619, "temperate": 4246, "virulent": 2196, "prophage_flagged": 4190,
                "host_specific": 8755, "shared": 219, "core": 6, "rvd_species_matched_pct": 19.8, "rgvc_species_matched_pct": 10.5})
    with open(f"{OUT}/web_overview_stats.json", "w") as f: json.dump(ov, f, indent=1)
    log(f"web_overview_stats.json written; sites found: {len(site_counts)} (paper 13)")

    # ---------- wait for the background FASTA step, then write the genome table (needs size/contigs/GC)
    if fasta_thread is not None:
        if fasta_thread.is_alive(): log("Tables done; waiting for the FASTA step to finish")
        fasta_thread.join()
        if fasta_err: print(f"STOP: FASTA step failed: {fasta_err[0]}"); sys.exit(1)
        missing = len(tq) - len(stats)
        if missing: print(f"WARNING: {missing} genomes were not gzipped")
    write_genomes()
    log(f"Done. Output in {OUT}")

if __name__ == "__main__":
    main()
