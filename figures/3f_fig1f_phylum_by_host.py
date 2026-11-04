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
#  SCRIPT  : 3f_fig1f_phylum_by_host.py  -  Fig. 1f, phylum composition per host species
#  AUTHOR  : QRQ
#  INPUT   : catalog_tables/web_genomes.tsv
#  OUTPUT  : fig1/Fig1f_phylum_by_host.{png,pdf,svg}, Fig1f_phylum_by_host.csv
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

g = pd.read_csv(urgg_in("catalog_tables/web_genomes.tsv"), sep="\t",
                usecols=["urgg_id", "rep_host_species", "phylum"])
g["phylum"] = g.phylum.fillna("NA")
top = g.phylum.value_counts().index[:7].tolist()
cats = top + ["Other"]
g["ph"] = np.where(g.phylum.isin(top), g.phylum, "Other")
vc = g.rep_host_species.value_counts()
hosts = sorted(vc.index, key=lambda h: (-vc[h], h))
tab = pd.crosstab(g.rep_host_species, g.ph).reindex(index=hosts, columns=cats, fill_value=0)
pct = 100 * tab.div(tab.sum(axis=1), axis=0)
out = urgg_out("fig1")
tab.assign(total=tab.sum(axis=1)).join(pct.add_suffix(" (%)").round(2)).to_csv(out / "Fig1f_phylum_by_host.csv")

set_font(size=8)
fig, ax = plt.subplots(figsize=(7.1, 5.2))
fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.36)
x = np.arange(len(hosts)); bottom = np.zeros(len(hosts))
for c, col in zip(cats, PALETTE8):
    ax.bar(x, pct[c].values, bottom=bottom, color=col, edgecolor="white", linewidth=0.4, width=0.7, label=c)
    bottom += pct[c].values
ax.set_xticks(x)
ax.set_xticklabels(hosts, rotation=45, ha="right", rotation_mode="anchor", style="italic", fontsize=8.5)
ax.set_ylim(0, 100); ax.set_ylabel("Representatives by phylum (%)", fontsize=10.5)
ax.set_xlabel("Host species", fontsize=10.5, labelpad=2)
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="x", length=0)
fig.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.54, 0.0), frameon=False, title="Phylum")
save(fig, out, "Fig1f_phylum_by_host", tight=False)
