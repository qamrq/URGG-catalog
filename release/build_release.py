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
#  SCRIPT  : build_release.py  -  Download files for the website and OMIX
#  AUTHOR  : QRQ
#  INPUT   : website_release/, drep_analysis/, annotation/gene_catalog/
#  OUTPUT  : website_release/release/ (archives, gene catalog, matrices, MD5SUMS.txt)
# ==========================================================

URGG_PROJECT = __import__("os").environ.get("URGG_PROJECT", "/path/to/metagenome_analysis")

import argparse, collections, glob, gzip, io, json, os, re, shutil, subprocess, sys, tarfile, time
from multiprocessing import Process
import pandas as pd

B   = URGG_PROJECT
WR  = f"{B}/website_release"
WEB = f"{WR}/web_data"
OUT = f"{WR}/release"
P = dict(
    tax_qc      = f"{WR}/urgg_tax_QC_table.tsv",
    species_dir = f"{WR}/genomes",
    strain_dir  = f"{B}/drep_analysis/drep_99_strain/dereplicated_genomes",
    strain_map  = f"{B}/drep_analysis/strain_id_map.tsv",
    gene_dir    = f"{B}/annotation/gene_catalog",
    clusters    = f"{B}/annotation/gene_catalog/nr_gene_catalog_cluster.tsv",
    virome_dir  = f"{B}/phase_5_virome",
    thiamine    = f"{B}/annotation/results_phase3/s63_thiamine_guild_propagated/63_genome_thiamine_guilds.tsv",
    amr         = f"{B}/annotation/amr_results_with_hits/all_mags_amr_combined.tsv",
    vf          = f"{B}/annotation/results_phase3/s64_vfdb_propagated/64_vf_presence_propagated.tsv",
    members     = f"{WEB}/web_members.tsv",
    strains_web = f"{WEB}/web_strains.tsv",
    votus_web   = f"{WEB}/web_votus.tsv",
    ko_web      = f"{WEB}/web_ko_genomes.tsv.gz",
    cazy_web    = f"{WEB}/web_cazy_genomes.tsv.gz",
    genomes_web = f"{WEB}/web_genomes.tsv",
)
EXPECT = dict(species=25451, strains=23999, genes=37129697, votus=8980, predicted_genes=46918358)
AMR_COLS = ["genome", "protein_id", "contig_id", "start", "stop", "strand", "gene_symbol", "sequence_name", "scope",
            "element_type", "element_subtype", "class", "subclass", "method", "target_length", "ref_length",
            "pct_coverage", "pct_identity", "alignment_length", "closest_accession", "closest_name", "hmm_id", "hmm_description"]

PIGZ = shutil.which("pigz")
def log(tag, m): print(time.strftime("%H:%M:%S"), f"[{tag}]", m, flush=True)

class GzOut:
    """Binary writer that compresses through pigz into <dst>.part, renamed to <dst> on close()."""
    def __init__(self, dst, threads):
        self.dst, self.part = dst, dst + ".part"
        self.fo = open(self.part, "wb")
        if PIGZ:
            self.p = subprocess.Popen([PIGZ, "-6", "-p", str(max(1, threads)), "-c"], stdin=subprocess.PIPE,
                                      stdout=self.fo, bufsize=8 << 20)
            self.w = self.p.stdin
        else:
            self.p, self.w = None, gzip.GzipFile(fileobj=self.fo, mode="wb", compresslevel=6)
    def write(self, b): self.w.write(b)
    def close(self):
        self.w.close()
        if self.p and self.p.wait() != 0: raise RuntimeError(f"pigz failed for {self.dst}")
        self.fo.close(); os.replace(self.part, self.dst)

def done(path, force):
    return os.path.exists(path) and not force

def write_plain(dst, df):
    df.to_csv(dst + ".part", sep="\t", index=False); os.replace(dst + ".part", dst)

def fasta_records(path):
    """Yield (header_without_>, sequence_bytes) from a plain or .gz FASTA."""
    op = gzip.open if path.endswith(".gz") else open
    h, seq = None, []
    with op(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if h is not None: yield h, b"".join(seq)
                h, seq = line[1:].rstrip(b"\r\n"), []
            else:
                seq.append(line.rstrip(b"\r\n"))
        if h is not None: yield h, b"".join(seq)

def wrap(seq, w=60):
    return b"\n".join(seq[i:i + w] for i in range(0, len(seq), w)) + b"\n"

def find_one(patterns, label, exclude=()):
    hits = sorted({p for pat in patterns for p in glob.glob(pat, recursive=True)
                   if os.path.isfile(p) and not any(e in os.path.basename(p) for e in exclude)})
    return hits

def tar_add(tf, name, data, mtime):
    ti = tarfile.TarInfo(name); ti.size = len(data); ti.mtime = mtime; ti.mode = 0o644
    tf.addfile(ti, io.BytesIO(data))

def strain_file(d, orig):
    o = str(orig)
    for cand in (o, o + ".fa", o + ".fasta", o + ".fna", re.sub(r"\.(fa|fasta|fna)$", "", o) + ".fa"):
        p = os.path.join(d, cand)
        if os.path.isfile(p): return p
    return None

def job_species(threads, force):
    T = "species"; dst = f"{OUT}/URGG_species_genomes.tar.gz"
    if done(dst, force): log(T, "exists, skipped"); return
    ids = pd.read_csv(P["tax_qc"], sep="\t", usecols=["public_id"]).public_id.tolist()
    miss = [u for u in ids if not os.path.isfile(f"{P['species_dir']}/{u}.fa")]
    if len(ids) != EXPECT["species"] or miss:
        raise SystemExit(f"species: {len(ids)} IDs, {len(miss)} FASTA missing (e.g. {miss[:3]})")
    log(T, f"packing {len(ids):,} genomes, pigz -p {threads}")
    out = GzOut(dst, threads); now = int(time.time())
    with tarfile.open(fileobj=out.w, mode="w|") as tf:
        for i, u in enumerate(ids, 1):
            with open(f"{P['species_dir']}/{u}.fa", "rb") as f:
                data = b"\n".join(l.replace(b"|", b"_") if l.startswith(b">") else l for l in f.read().split(b"\n"))
            tar_add(tf, f"URGG_species_genomes/{u}.fa", data, now)
            if i % 5000 == 0: log(T, f"{i:,}")
    out.close(); log(T, f"done -> {os.path.basename(dst)}")

def job_strains(threads, force):
    T = "strains"; dst = f"{OUT}/URGG_strain_genomes.tar.gz"; cmap = f"{OUT}/URGG_strain_contig_ids.tsv.gz"
    if done(dst, force) and done(cmap, force): log(T, "exists, skipped"); return
    sm = pd.read_csv(P["strain_map"], sep="\t", low_memory=False)
    files = [(s, strain_file(P["strain_dir"], o)) for s, o in zip(sm.strain_id, sm.original_genome_id)]
    miss = [s for s, f in files if f is None]
    if len(sm) != EXPECT["strains"] or miss:
        raise SystemExit(f"strains: {len(sm)} rows, {len(miss)} FASTA not found in strain_dir (e.g. {miss[:3]})")
    log(T, f"packing {len(files):,} genomes (contigs -> <URGGS>_<n>), pigz -p {threads}")
    out = GzOut(dst, max(1, threads - 2)); cm = GzOut(cmap, 2); now = int(time.time())
    cm.write(b"contig_id\tstrain_id\toriginal_header\n")
    with tarfile.open(fileobj=out.w, mode="w|") as tf:
        for i, (s, path) in enumerate(files, 1):
            buf, k = [], 0
            for h, seq in fasta_records(path):
                k += 1; new = f"{s}_{k}".encode()
                buf.append(b">" + new + b"\n" + wrap(seq))
                cm.write(new + b"\t" + s.encode() + b"\t" + h + b"\n")
            tar_add(tf, f"URGG_strain_genomes/{s}.fa", b"".join(buf), now)
            if i % 5000 == 0: log(T, f"{i:,}")
    out.close(); cm.close(); log(T, f"done -> {os.path.basename(dst)}, {os.path.basename(cmap)}")

def job_genes(threads, force, faa, fna):
    T = "genes"
    f_faa, f_fna = f"{OUT}/URGG_gene_catalog.faa.gz", f"{OUT}/URGG_gene_catalog.fna.gz"
    f_map, f_cl = f"{OUT}/URGG_gene_id_map.tsv.gz", f"{OUT}/URGG_gene_clusters.tsv.gz"
    if all(done(x, force) for x in (f_faa, f_map, f_cl)) and (not fna or done(f_fna, force)): log(T, "exists, skipped"); return
    mem = pd.read_csv(P["members"], sep="\t", usecols=["original_genome_id", "urgg_id", "strain_id"], dtype=str).fillna("")
    g2sp = dict(zip(mem.original_genome_id, mem.urgg_id)); g2st = dict(zip(mem.original_genome_id, mem.strain_id))
    del mem

    log(T, "pass 0: reading cluster representatives")
    reps = set()
    with open(P["clusters"], "rb") as f:
        for line in f: reps.add(line.split(b"\t", 1)[0])
    log(T, f"  {len(reps):,} representatives (paper {EXPECT['genes']:,}){'' if len(reps) == EXPECT['genes'] else '  <-- CHECK'}")

    log(T, f"pass 1: {os.path.basename(faa)} -> {os.path.basename(f_faa)} + ID map (pigz -p {threads})")
    gid = {}
    o_faa, o_map = GzOut(f_faa, threads - 4), GzOut(f_map, 4)
    o_map.write(b"gene_id\toriginal_gene_id\tsource_genome\tspecies_urgg_id\tstrain_id\n")
    unmapped = n = scanned = 0
    for scanned, (h, seq) in enumerate(fasta_records(faa), 1):
        orig = h.split(None, 1)[0]
        if orig not in reps or orig in gid: continue
        n += 1
        new = b"URGGG%08d" % n; gid[orig] = n
        genome = orig.split(b"~", 1)[0].decode()
        sp = g2sp.get(genome, "");  unmapped += sp == ""
        o_faa.write(b">" + new + b"\n" + wrap(seq.rstrip(b"*")))
        o_map.write(new + b"\t" + orig + b"\t" + genome.encode() + b"\t" + sp.encode() + b"\t" + g2st.get(genome, "").encode() + b"\n")
        if n % 5_000_000 == 0: log(T, f"  faa {n:,}")
    o_faa.close(); o_map.close()
    log(T, f"  scanned {scanned:,} proteins; kept {n:,} representatives (paper {EXPECT['genes']:,})"
           f"{'' if n == EXPECT['genes'] else '  <-- CHECK'}; representatives with no protein: {len(reps) - n:,}; "
           f"genes whose genome has no species ID: {unmapped:,}")
    del reps

    # pass 2: nucleotide FASTA, same IDs
    if not fna:
        log(T, "pass 2 skipped: no nucleotide FASTA given (--fna); protein catalog only")
    else:
      log(T, f"pass 2: {os.path.basename(fna)} -> {os.path.basename(f_fna)}")
      recs = []; missing = 0
      o_fna = GzOut(f_fna, threads)
      for m, (h, seq) in enumerate(fasta_records(fna), 1):
          k = gid.get(h.split(None, 1)[0])
          if k is None: missing += 1; continue
          o_fna.write(b">" + (b"URGGG%08d" % k) + b"\n" + wrap(seq))
          if m % 5_000_000 == 0: log(T, f"  fna {m:,}")
      o_fna.close()
      log(T, f"  {m:,} nucleotide records scanned; written {m - missing:,} (expected {len(gid):,}); non-representative skipped: {missing:,}")

    # pass 3: clusters (rep -> member)
    log(T, f"pass 3: clusters -> {os.path.basename(f_cl)}")
    o_cl = GzOut(f_cl, threads); o_cl.write(b"gene_id\tmember_original_gene_id\tmember_species_urgg_id\n")
    nl = norep = 0
    with open(P["clusters"], "rb") as f:
        for nl, line in enumerate(f, 1):
            rep, member = line.rstrip(b"\r\n").split(b"\t", 1)
            k = gid.get(rep)
            if k is None: norep += 1; continue
            o_cl.write((b"URGGG%08d" % k) + b"\t" + member + b"\t" + g2sp.get(member.split(b"~", 1)[0].decode(), "").encode() + b"\n")
            if nl % 10_000_000 == 0: log(T, f"  clusters {nl:,}")
    o_cl.close()
    log(T, f"  {nl:,} member genes (paper predicted genes {EXPECT['predicted_genes']:,}); rep not in catalog: {norep:,}")
    log(T, "done")

def job_small(threads, force, votu_fna, wait_hours):
    T = "tables"
    tq = pd.read_csv(P["tax_qc"], sep="\t", usecols=["public_id", "original_genome_id"])
    o2u = dict(zip(tq.original_genome_id, tq.public_id))

    # vOTU sequences
    dst = f"{OUT}/URGG_vOTUs.fna.gz"
    if not done(dst, force):
        if not (votu_fna and os.path.exists(P["votus_web"])):
            log(T, "WARNING vOTU FASTA or web_votus.tsv missing -> URGG_vOTUs.fna.gz skipped")
        else:
            vt = pd.read_csv(P["votus_web"], sep="\t", dtype=str)
            r2v = dict(zip(vt.representative, vt.votu_id))
            sm = pd.read_csv(P["strain_map"], sep="\t", dtype=str)
            p2s = dict(zip(sm.pipeline_index, sm.strain_id))
            p2s.update({str(int(re.sub(r"\D", "", k))): s for k, s in p2s.items() if re.sub(r"\D", "", str(k))})
            strip = lambda name: name.rsplit("__", 1)[0] if name.count("__") >= 2 else name
            short = collections.Counter(strip(r) for r in r2v)
            r2v_short = {strip(r): v for r, v in r2v.items() if short[strip(r)] == 1}
            def lookup(hdr):
                if hdr in r2v: return r2v[hdr]
                pre, _, rest = hdr.partition("__")
                s = p2s.get(pre) or (p2s.get(str(int(re.sub(r"\D", "", pre)))) if re.sub(r"\D", "", pre) else None)
                if not s: return None
                cand = f"{s}__{rest}"
                return r2v.get(cand) or r2v_short.get(strip(cand))
            out = GzOut(dst, threads); hm = GzOut(f"{OUT}/URGG_vOTU_header_map.tsv.gz", 1)
            hm.write(b"votu_id\toriginal_header\n"); seen, unmatched = set(), []
            for h, seq in fasta_records(votu_fna):
                hdr = h.split(None, 1)[0].decode()
                v = lookup(hdr)
                if v is None: unmatched.append(hdr); continue
                if v in seen: continue
                seen.add(v); out.write(b">" + v.encode() + b"\n" + wrap(seq)); hm.write(f"{v}\t{hdr}\n".encode())
            out.close(); hm.close()
            log(T, f"URGG_vOTUs.fna.gz: {len(seen):,} vOTUs (paper {EXPECT['votus']:,}){'' if len(seen) == EXPECT['votus'] else '  <-- CHECK'}; "
                   f"unmatched headers: {len(unmatched):,}" + (f" e.g. {unmatched[:3]}" if unmatched else ""))
            if len(seen) != EXPECT["votus"]:
                log(T, f"  example Suppl. Table 15 names: {list(r2v)[:3]}")
    if os.path.exists(P["votus_web"]) and not done(f"{OUT}/URGG_vOTU_table.tsv", force):
        shutil.copy(P["votus_web"], f"{OUT}/URGG_vOTU_table.tsv")

    if not done(f"{OUT}/URGG_thiamine_guilds.tsv", force):
        t = pd.read_csv(P["thiamine"], sep="\t", low_memory=False)
        t.insert(0, "urgg_id", t.genome.map(o2u)); t = t[t.urgg_id.notna()].rename(columns={"genome": "original_genome_id"})
        write_plain(f"{OUT}/URGG_thiamine_guilds.tsv", t.sort_values("urgg_id")); log(T, f"URGG_thiamine_guilds.tsv: {len(t):,} rows")
    if not done(f"{OUT}/URGG_AMR.tsv", force):
        a = pd.read_csv(P["amr"], sep="\t", header=None, dtype=str, keep_default_na=False)
        a.columns = AMR_COLS[:a.shape[1]] + [f"col{i}" for i in range(len(AMR_COLS), a.shape[1])]
        a.insert(0, "urgg_id", a.genome.map(o2u)); a = a[a.urgg_id.notna()].rename(columns={"genome": "original_genome_id"})
        write_plain(f"{OUT}/URGG_AMR.tsv", a.sort_values("urgg_id"))
        log(T, f"URGG_AMR.tsv: {len(a):,} hits in {a.urgg_id.nunique():,} species genomes (paper 294 AMR-positive)")
    if not done(f"{OUT}/URGG_VFDB.tsv", force):
        v = pd.read_csv(P["vf"], sep="\t", low_memory=False)
        gcol = max(v.columns, key=lambda c: v[c].astype(str).isin(o2u.keys()).mean())
        v.insert(0, "urgg_id", v[gcol].map(o2u)); v = v[v.urgg_id.notna()].rename(columns={gcol: "original_genome_id"})
        write_plain(f"{OUT}/URGG_VFDB.tsv", v.sort_values("urgg_id")); log(T, f"URGG_VFDB.tsv: {len(v):,} rows")
    if os.path.exists(P["strains_web"]) and not done(f"{OUT}/URGG_strain_metadata.tsv", force):
        shutil.copy(P["strains_web"], f"{OUT}/URGG_strain_metadata.tsv")

    def web_finished():
        return (os.path.exists(P["genomes_web"]) and os.path.exists(P["members"])
                and os.path.getmtime(P["genomes_web"]) >= os.path.getmtime(P["members"]))
    t0 = time.time()
    if not web_finished(): log(T, "waiting for build_web_tables.py to finish (KO/CAZy matrices) ...")
    while not web_finished():
        if time.time() - t0 > wait_hours * 3600:
            log(T, "WARNING gave up waiting; rerun later with the same command to add the matrices"); return
        time.sleep(60)
    for src, dst, col in [(P["ko_web"], f"{OUT}/URGG_KO_matrix.tsv.gz", "ko"), (P["cazy_web"], f"{OUT}/URGG_CAZy_matrix.tsv.gz", "cazy_family")]:
        if done(dst, force): continue
        if not os.path.exists(src): log(T, f"WARNING {src} missing, {os.path.basename(dst)} skipped"); continue
        by_g = collections.defaultdict(list)
        with gzip.open(src, "rt") as f:
            next(f)
            for line in f:
                feat, ids = line.rstrip("\n").split("\t")
                for u in ids.split(";"):
                    if u: by_g[u].append(feat)
        out = GzOut(dst, threads); out.write(f"urgg_id\t{col}\n".encode()); n = 0
        for u in sorted(by_g):
            for feat in sorted(by_g[u]): out.write(f"{u}\t{feat}\n".encode()); n += 1
        out.close(); log(T, f"{os.path.basename(dst)}: {n:,} rows, {len(by_g):,} genomes")
    log(T, "done")

def run(fn, *args):
    try: fn(*args)
    except BaseException as e:
        print(time.strftime("%H:%M:%S"), f"[{fn.__name__}] FAILED: {e!r}", flush=True); sys.exit(1)

def md5_all():
    files = sorted(f for f in os.listdir(OUT) if not f.endswith(".part") and f != "MD5SUMS.txt")
    procs = [(f, subprocess.Popen(["md5sum", f], cwd=OUT, stdout=subprocess.PIPE, text=True)) for f in files]
    lines = [p.communicate()[0] for _, p in procs]
    with open(f"{OUT}/MD5SUMS.txt", "w") as fo: fo.writelines(lines)
    log("md5", f"MD5SUMS.txt: {len(files)} files")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--faa"); ap.add_argument("--fna"); ap.add_argument("--votu-fna")
    ap.add_argument("--only", choices=["species", "strains", "genes", "tables"], action="append")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--wait-hours", type=float, default=24)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    # ---- locate inputs
    gd = P["gene_dir"]
    faa = a.faa or find_one([f"{gd}/*.faa", f"{gd}/*.faa.gz", f"{gd}/*protein*.fa*"], "faa")
    fna = a.fna or find_one([f"{gd}/*.fna", f"{gd}/*.fna.gz", f"{gd}/*.ffn", f"{gd}/*.ffn.gz", f"{gd}/*nucl*.fa*"], "fna")
    vfa = a.votu_fna or find_one([f"{P['virome_dir']}/**/*[vV][oO][tT][uU]*.f*a", f"{P['virome_dir']}/**/*[vV][oO][tT][uU]*.f*a.gz",
                                  f"{P['virome_dir']}/**/*[vV][oO][tT][uU]*.fasta"], "votu", exclude=("faa",))
    print("Inputs:")
    problems = []
    for lab, v in [("gene protein FASTA (--faa)", faa), ("gene nucleotide FASTA (--fna)", fna), ("vOTU FASTA (--votu-fna)", vfa)]:
        if isinstance(v, list):
            if len(v) == 1: v = v[0]
            else:
                problems.append(lab); print(f"  {lab}: {'NOT FOUND' if not v else 'several candidates, pass one:'}")
                for x in v: print(f"      {x}  ({os.path.getsize(x)/1e9:.1f} GB)")
                continue
        print(f"  {lab}: {v}  ({os.path.getsize(v)/1e9:.1f} GB)" if os.path.exists(v) else f"  {lab}: {v}  NOT FOUND")
        if lab.startswith("gene protein"): faa = v
        elif lab.startswith("gene nucleotide"): fna = v
        else: vfa = v
    sm = pd.read_csv(P["strain_map"], sep="\t", low_memory=False)
    nf = sum(strain_file(P["strain_dir"], o) is not None for o in sm.original_genome_id)
    print(f"  strain FASTAs found: {nf:,} / {len(sm):,}  ({P['strain_dir']})")
    ns = len(glob.glob(f"{P['species_dir']}/URGG*.fa"))
    print(f"  species FASTAs found: {ns:,} / {EXPECT['species']:,}")
    for k in ("members", "votus_web"):
        print(f"  {os.path.basename(P[k])}: {'ok' if os.path.exists(P[k]) else 'MISSING'}")
    try:
        memgb = int(open("/proc/meminfo").read().split("MemAvailable:")[1].split()[0]) / 1e6
        print(f"  available RAM: {memgb:.0f} GB (gene ID table needs ~10-15 GB)")
    except Exception: pass
    print(f"  pigz: {PIGZ or 'NOT FOUND (python gzip fallback, slow)'};  output: {OUT}")
    if a.dry_run: return

    jobs = a.only or ["species", "strains", "genes", "tables"]
    if "genes" in jobs and isinstance(fna, list): fna = None
    if "genes" in jobs and (not isinstance(faa, str) or not os.path.exists(P["members"])):
        sys.exit("STOP: protein FASTA not resolved (see above); pass --faa")
    os.makedirs(OUT, exist_ok=True)
    T = a.threads
    share = {"species": max(2, T // 4), "strains": max(2, T // 4), "genes": max(6, T * 3 // 8), "tables": max(2, T // 8)}
    target = {"species": (job_species, share["species"], a.force),
              "strains": (job_strains, share["strains"], a.force),
              "genes":   (job_genes, share["genes"], a.force, faa, fna),
              "tables":  (job_small, share["tables"], a.force, vfa if isinstance(vfa, str) else None, a.wait_hours)}
    log("main", "starting " + ", ".join(f"{j} (pigz {share[j]} thr)" for j in jobs))
    procs = {j: Process(target=run, args=target[j], name=j) for j in jobs}
    for p in procs.values(): p.start()
    for p in procs.values(): p.join()
    bad = [j for j, p in procs.items() if p.exitcode != 0]
    md5_all()
    if bad: print(f"\nFINISHED WITH ERRORS in: {bad}. Fix and rerun the same command (finished files are kept)."); sys.exit(1)
    log("main", f"all done -> {OUT}")

if __name__ == "__main__":
    main()
