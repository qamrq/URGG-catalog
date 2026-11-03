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
#  SCRIPT  : 86b_species_accumulation.py  -  Species accumulation curves and Chao1
#  AUTHOR  : QRQ
#  INPUT   : genome_to_species.tsv (arguments)
#  OUTPUT  : 86_curves_by_genome.tsv, 86_curves_by_sample.tsv, 86_host_summary.tsv, 86_accumulation_*.pdf
# ==========================================================

import argparse, os, sys, warnings, logging
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--membership", required=True)
ap.add_argument("--sep", default="\t")
ap.add_argument("--genome-col", default="genome")
ap.add_argument("--cluster-col", default="species_rep")
ap.add_argument("--host-col", default="host_species")
ap.add_argument("--sample-col", default=None)
ap.add_argument("--perm", type=int, default=100)
ap.add_argument("--min-genomes", type=int, default=100, help="minimum genomes for a figure panel")
ap.add_argument("--min-chao", type=int, default=1000,
                help="minimum genomes for Chao1 to be reported as reliable (flag column + panel label)")
ap.add_argument("--italic", choices=["auto", "yes", "no"], default="auto",
                help="italic host titles; auto = italic only when --host-col is host_species")
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out", default="s86_accumulation")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng(a.seed)

d0 = pd.read_csv(a.membership, sep=a.sep, dtype=str, keep_default_na=True)
need = [a.genome_col, a.cluster_col, a.host_col] + ([a.sample_col] if a.sample_col else [])
miss = [c for c in need if c not in d0.columns]
if miss:
    sys.exit(f"Missing columns {miss}. Found: {list(d0.columns)}")
d = d0[need].copy()
for c in need:
    d[c] = d[c].str.replace("\r", "", regex=False).str.strip().replace({"": np.nan, "nan": np.nan, "NA": np.nan})

core = [a.genome_col, a.cluster_col, a.host_col]
n_in = len(d)
d = d.dropna(subset=core)
n_na = n_in - len(d)
dup = d[a.genome_col].duplicated(keep=False)
if dup.any():
    conflict = d[dup].groupby(a.genome_col)[a.cluster_col].nunique().gt(1).sum()
    if conflict:
        sys.exit(f"{conflict} genomes appear more than once with DIFFERENT clusters — fix the input.")
d = d.drop_duplicates(a.genome_col)
n_dup = n_in - n_na - len(d)
d[a.host_col] = d[a.host_col].str.replace("_", " ", regex=False)
print(f"{n_in:,} rows read; dropped {n_na:,} with missing genome/cluster/host and {n_dup:,} duplicates")
print(f"{len(d):,} genomes, {d[a.cluster_col].nunique():,} species clusters, {d[a.host_col].nunique()} hosts")
if not 110000 <= len(d) <= 120000:
    print("WARNING: expected ~117,145 genomes; the table should hold all quality-filtered genomes, "
          "not only representatives.")
if a.sample_col:
    n_s = d[a.sample_col].notna().sum()
    print(f"{n_s:,} genomes ({100*n_s/len(d):.1f}%) carry a sample ID; the rest are excluded from sample-based curves only")

def accumulate(unit_idx, clus_idx, n_units, n_perm):
    """Exact permutation accumulation.
    unit_idx[i], clus_idx[i]: the i-th (unit, cluster) incidence pair (unique pairs).
    For a random order of units, the cluster count after k units = number of
    clusters whose earliest unit has rank < k. Returns (x, mean, lo, hi, all)."""
    _, clus_idx = np.unique(clus_idx, return_inverse=True)
    n_clus = clus_idx.max() + 1
    out = np.empty((n_perm, n_units), dtype=np.int32)
    for p in range(n_perm):
        rank = rng.permutation(n_units)[unit_idx]
        first = np.full(n_clus, n_units, dtype=np.int64)
        np.minimum.at(first, clus_idx, rank)
        out[p] = np.cumsum(np.bincount(first, minlength=n_units)[:n_units])
    return (np.arange(1, n_units + 1), out.mean(0),
            np.percentile(out, 2.5, 0), np.percentile(out, 97.5, 0), out)

def chao1_bc(S, f1, f2, n):
    """bias-corrected Chao1 (Chao 2005; vegan::estimateR): defined for f2 = 0 too."""
    return S + ((n - 1) / n) * f1 * (f1 - 1) / (2 * (f2 + 1))

def tail_slope(mean, per):
    """new species per `per` units over the last 10% of the curve"""
    n = len(mean)
    if n < 10:
        return np.nan
    k = int(np.floor(0.9 * n)) - 1
    return per * (mean[-1] - mean[k]) / (n - 1 - k)

def thin(n, target=400):
    step = max(1, n // target)
    return np.unique(np.r_[np.arange(0, n, step), n - 1])

hosts = d[a.host_col].value_counts()
groups = [("All hosts", d)] + [(h, d[d[a.host_col] == h]) for h in hosts.index]

rows_g, rows_s, summ = [], [], []
for h, g in groups:
    n = len(g)
    c = g[a.cluster_col].to_numpy()
    # genome-based: each genome is its own unit
    x, m, lo, hi, _ = accumulate(np.arange(n), c, n, a.perm)
    rows_g += [dict(host=h, n_genomes=int(x[i]), species_mean=round(float(m[i]), 2),
                    species_lo=float(lo[i]), species_hi=float(hi[i])) for i in thin(n)]
    cnt = pd.Series(c).value_counts()
    S, f1, f2 = len(cnt), int((cnt == 1).sum()), int((cnt == 2).sum())
    ch1 = chao1_bc(S, f1, f2, n)
    s = dict(host=h, genomes=n, species_observed=S, singletons_f1=f1, doubletons_f2=f2,
             chao1=round(ch1, 1), observed_of_chao1_pct=round(100 * S / ch1, 1),
             goods_coverage_pct=round(100 * (1 - f1 / n), 1),
             new_species_per_100_genomes_last10pct=round(tail_slope(m, 100), 2),
             chao1_reliable=bool(n >= a.min_chao and f2 > 0))

    if a.sample_col:
        gs = g.dropna(subset=[a.sample_col])
        pairs = gs[[a.sample_col, a.cluster_col]].drop_duplicates()
        s["genomes_with_sample_id"] = len(gs)
        s["samples"] = pairs[a.sample_col].nunique()
        if s["samples"] >= 2:
            u_codes, u_idx = np.unique(pairs[a.sample_col].to_numpy(), return_inverse=True)
            mu = len(u_codes)
            xs, ms, los, his, _ = accumulate(u_idx, pairs[a.cluster_col].to_numpy(), mu, a.perm)
            rows_s += [dict(host=h, n_samples=int(xs[i]), species_mean=round(float(ms[i]), 2),
                            species_lo=float(los[i]), species_hi=float(his[i])) for i in thin(mu)]
            inc = pairs[a.cluster_col].value_counts()
            Ss, q1, q2 = len(inc), int((inc == 1).sum()), int((inc == 2).sum())
            ch2 = chao1_bc(Ss, q1, q2, mu)
            s.update(species_observed_in_sampled_runs=Ss, uniques_Q1=q1, duplicates_Q2=q2,
                     chao2=round(ch2, 1), observed_of_chao2_pct=round(100 * Ss / ch2, 1),
                     new_species_per_10_samples_last10pct=round(tail_slope(ms, 10), 2))
    summ.append(s)
    extra = f"  runs {s.get('samples', 0):>6,}" if a.sample_col else ""
    print(f"{h:28s} genomes {n:>7,}  species {S:>6,}  Chao1 {ch1:>9,.0f}  "
          f"obs/Chao1 {s['observed_of_chao1_pct']:>5}%  Good's {s['goods_coverage_pct']:5.1f}%{extra}")

cg = pd.DataFrame(rows_g); cs = pd.DataFrame(rows_s); sm = pd.DataFrame(summ)
cg.to_csv(f"{a.out}/86_curves_by_genome.tsv", sep="\t", index=False)
if len(cs):
    cs.to_csv(f"{a.out}/86_curves_by_sample.tsv", sep="\t", index=False)
sm.to_csv(f"{a.out}/86_host_summary.tsv", sep="\t", index=False)
pd.DataFrame([dict(membership=a.membership, rows_read=n_in, dropped_missing=n_na, dropped_duplicates=n_dup,
                   genomes_used=len(d), species=d[a.cluster_col].nunique(), hosts=d[a.host_col].nunique(),
                   host_col=a.host_col, sample_col=a.sample_col or "", permutations=a.perm, seed=a.seed,
                   chao_estimator="bias-corrected: S + ((n-1)/n) f1(f1-1) / (2(f2+1))")]
             ).T.to_csv(f"{a.out}/86_run_info.tsv", sep="\t", header=False)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message="findfont")
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
italic = a.italic == "yes" or (a.italic == "auto" and a.host_col == "host_species")

def panel_fig(curves, xcol, show, label_fn, xlab, fname):
    nc = 4; nr = int(np.ceil(len(show) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(7.1, 1.65 * nr), squeeze=False, constrained_layout=True)
    for ax, h in zip(axs.flat, show):
        c = curves[curves["host"] == h]
        ax.fill_between(c[xcol], c.species_lo, c.species_hi, color="#9ecae1", lw=0)
        ax.plot(c[xcol], c.species_mean, color="#08519c", lw=1)
        ax.set_title(h, fontsize=7, pad=2, style="italic" if (italic and h != "All hosts") else "normal")
        ax.text(0.03, 0.97, label_fn(sm.loc[sm.host == h].iloc[0]), transform=ax.transAxes, va="top", fontsize=6)
        ax.tick_params(length=2, width=0.5, labelsize=6)
        kfmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:g}k" if v >= 1000 else f"{v:g}")
        ax.xaxis.set_major_formatter(kfmt); ax.yaxis.set_major_formatter(kfmt)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.set_ylim(bottom=0); ax.set_xlim(left=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axs.flat[len(show):]:
        ax.axis("off")
    fig.supxlabel(xlab, fontsize=7); fig.supylabel("Species-level clusters", fontsize=7)
    fig.savefig(f"{a.out}/{fname}.pdf"); fig.savefig(f"{a.out}/{fname}.png", dpi=300)
    plt.close(fig)

show_g = [h for h, g in groups if len(g) >= a.min_genomes]
panel_fig(cg, "n_genomes", show_g,
          lambda r: f"{int(r.species_observed):,} species\n" +
                    (f"obs/Chao1 {r.observed_of_chao1_pct}%" if r.chao1_reliable else f"n = {int(r.genomes):,} genomes"),
          "Genomes sampled", "86_accumulation_by_genome")
if len(cs):
    show_s = [h for h in show_g if h in set(cs.host)]
    panel_fig(cs, "n_samples", show_s,
              lambda r: f"{int(r.samples):,} runs\nobs/Chao2 {r.observed_of_chao2_pct}%",
              "Sequencing runs sampled", "86_accumulation_by_sample")
print(f"\nWritten to {a.out}/  (permutations = {a.perm} for every host, seed = {a.seed})")
