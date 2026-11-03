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
#  SCRIPT  : Fig1e_quality_density.py  -  Fig. 1e, completeness and contamination by tier
#  AUTHOR  : QRQ
#  INPUT   : website_release/web_data/web_genomes.tsv
#  OUTPUT  : fig1/Fig1e_quality_density.{png,pdf,svg}, Fig1e_tier_shares.csv
# ==========================================================

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_setup"))
from urgg_paths import urgg_in, urgg_out  # noqa: E402
from urgg_style import TIER_COL, TIER_LABEL, TIER_ORDER, novelty_tier, save, set_font  # noqa: E402

NOV_ORDER = ["Described species", "Novel species", "Novel genus or higher"]
NOV_STYLE = {"Described species": "-", "Novel species": "--", "Novel genus or higher": ":"}
SPECS = [("completeness", np.linspace(70, 100, 300), 0.9, "Completeness (%)"),
         ("contamination", np.linspace(0, 5, 300), 0.15, "Contamination (%)")]

g = pd.read_csv(urgg_in("website_release/web_data/web_genomes.tsv"), sep="\t",
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
