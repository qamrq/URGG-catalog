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
#  SCRIPT  : 2f_fig1e_quality_density.py  -  Fig. 1e, completeness and contamination by tier
#  AUTHOR  : QRQ
#  INPUT   : catalog_tables/web_genomes.tsv
#  OUTPUT  : fig1/Fig1e_quality_density.{png,pdf,svg}, Fig1e_tier_shares.csv
# ==========================================================

import os
from pathlib import Path

URGG_PROJECT = Path(os.environ.get("URGG_PROJECT", "data"))
URGG_FIGURES = Path(os.environ.get("URGG_FIGURES", "figures_out"))


def urgg_in(rel):
    p = URGG_PROJECT / rel
    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p} (set URGG_PROJECT)")
    return p


def urgg_out(fig):
    d = URGG_FIGURES / fig
    d.mkdir(parents=True, exist_ok=True)
    return d

import matplotlib as mpl
from matplotlib import font_manager

BLUE, YELLOW, ORANGE = "#2a78d6", "#eda100", "#eb6834"
TIER_ORDER = ["Host-specific", "Shared", "Core"]
TIER_COL = {"Host-specific": BLUE, "Shared": YELLOW, "Core": ORANGE}
TIER_LABEL = {"Host-specific": "Host-specific (1 host species)",
              "Shared": "Shared (2-3 host species)",
              "Core": "Core (>3 host species)"}
PALETTE8 = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
FAMILY_COL = {"Bovidae": "#2a78d6", "Cervidae": "#eb6834", "Giraffidae": "#1baf7a", "Camelidae": "#e87ba4"}

def novelty_tier(v):
    """web_genomes.tsv 'novelty' -> the three tiers used in Fig. 1e (genus or higher pooled)."""
    if v == "Described":
        return "Described species"
    if v == "Novel species":
        return "Novel species"
    return "Novel genus or higher"

def set_font(prefer=("Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"), size=7):
    import os
    avail = {f.name for f in font_manager.fontManager.ttflist}
    want = os.environ.get("URGG_FONT")
    fam = want if want in avail else next((f for f in prefer if f in avail), "serif")
    mpl.rcParams.update({"font.family": fam, "font.size": size, "axes.linewidth": 0.6,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    return fam

def save(fig, outdir, stem, dpi=600, tight=True):
    kw = dict(bbox_inches="tight", pad_inches=0.05) if tight else {}
    for ext in ("png", "pdf", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=dpi if ext == "png" else None, **kw)
    print("Written:", outdir / f"{stem}.png|pdf|svg")

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


NOV_ORDER = ["Described species", "Novel species", "Novel genus or higher"]
NOV_STYLE = {"Described species": "-", "Novel species": "--", "Novel genus or higher": ":"}
SPECS = [("completeness", np.linspace(70, 100, 300), 0.9, "Completeness (%)"),
         ("contamination", np.linspace(0, 5, 300), 0.15, "Contamination (%)")]

g = pd.read_csv(urgg_in("catalog_tables/web_genomes.tsv"), sep="\t",
                usecols=["urgg_id", "sharing_tier", "novelty", "completeness", "contamination"])
g["nov"] = g.novelty.map(novelty_tier)
N = len(g)
share = (100 * g.sharing_tier.value_counts() / N).reindex(TIER_ORDER)
nov = (100 * g.nov.value_counts() / N).reindex(NOV_ORDER)
out = urgg_out("fig1")
pd.concat([share.rename("pct").to_frame().assign(kind="sharing tier"),
           nov.rename("pct").to_frame().assign(kind="novelty tier")]).to_csv(out / "Fig1e_tier_shares.csv")

def kde(v, grid, bw):
    v = np.asarray(v, float)
    return np.exp(-0.5 * ((grid[:, None] - v[None, :]) / bw) ** 2).sum(axis=1) if v.size else np.zeros_like(grid)

set_font(size=9)
fig, axes = plt.subplots(2, 1, figsize=(3.6, 5.6))
fig.subplots_adjust(hspace=0.32)
for ax, (col, grid, bw, xlab) in zip(axes, SPECS):
    bands = {t: kde(g.loc[g.sharing_tier == t, col], grid, bw) for t in TIER_ORDER}
    total = sum(bands.values()); peak = total.max(); base = np.zeros_like(grid)
    for t in TIER_ORDER:
        top = base + bands[t] / peak
        ax.fill_between(grid, base, top, color=TIER_COL[t], lw=0, zorder=2); base = top
    ax.plot(grid, total / peak, color="black", lw=0.8, zorder=3)
    for t in NOV_ORDER:
        ax.plot(grid, kde(g.loc[g.nov == t, col], grid, bw) / peak, color="black",
                ls=NOV_STYLE[t], lw=1.0, zorder=4)
    ax.set_xlim(grid[0], grid[-1]); ax.set_ylim(0, 1.05)
    ax.set_xlabel(xlab); ax.set_ylabel("Relative density")
    ax.spines[["top", "right"]].set_visible(False)

h = [plt.Rectangle((0, 0), 1, 1, fc=TIER_COL[t], ec="none") for t in TIER_ORDER]
h += [plt.Line2D([0], [0], color="black", ls=NOV_STYLE[t], lw=1.0) for t in NOV_ORDER]
lab = [f"{TIER_LABEL[t]}: {share[t]:.1f}%" for t in TIER_ORDER] + [f"{t}: {nov[t]:.1f}%" for t in NOV_ORDER]
axes[1].legend(h, lab, loc="upper right", frameon=False, fontsize=7, handlelength=1.8)
save(fig, out, "Fig1e_quality_density")
print(share.round(1).to_dict(), nov.round(1).to_dict())
