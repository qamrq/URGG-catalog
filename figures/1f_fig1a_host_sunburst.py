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
#  SCRIPT  : 1f_fig1a_host_sunburst.py  -  Fig. 1a, genomes by host family and host species
#  AUTHOR  : QRQ
#  INPUT   : catalog_tables/web_genomes.tsv
#  OUTPUT  : fig1/Fig1a_host_sunburst.{png,pdf,svg}, Fig1a_host_counts.csv
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

import colorsys
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FAMILIES = ["Bovidae", "Cervidae", "Giraffidae", "Camelidae"]
N_BIG = 6
BIG_LABEL_ANG = [44, -22, -88, -105, -135, -160]
R_LAB, COL_X, COL_TOP, COL_STEP = 1.5, -1.85, 2.35, 0.24

g = pd.read_csv(urgg_in("catalog_tables/web_genomes.tsv"), sep="\t",
                usecols=["urgg_id", "rep_host_family", "rep_host_species"])
N = len(g)
cnt = (g.groupby(["rep_host_family", "rep_host_species"]).size().rename("n").reset_index())
cnt["fam_rank"] = cnt.rep_host_family.map({f: i for i, f in enumerate(FAMILIES)})
cnt = cnt.sort_values(["fam_rank", "n"], ascending=[True, False]).reset_index(drop=True)
cnt["pct"] = 100 * cnt.n / N
out = urgg_out("fig1")
cnt.drop(columns="fam_rank").to_csv(out / "Fig1a_host_counts.csv", index=False)
assert N == 25451, f"expected 25,451 representatives, found {N:,}"

def shade(hex_color, k, n):
    r, gg, b = mpl.colors.to_rgb(hex_color)
    h, _, s = colorsys.rgb_to_hls(r, gg, b)
    order = [0, 4, 2, 6, 1, 5, 3, 7, 8]
    pos = order[k % len(order)] / max(n - 1, 1)
    return colorsys.hls_to_rgb(h, 0.36 + 0.36 * min(pos, 1.0), s)

cols, seen = [], {f: 0 for f in FAMILIES}
per_fam = cnt.rep_host_family.value_counts()
for f in cnt.rep_host_family:
    cols.append(shade(FAMILY_COL[f], seen[f], per_fam[f]))
    seen[f] += 1
fam_n = cnt.groupby("rep_host_family").n.sum().reindex(FAMILIES)

set_font(prefer=("Arial", "Liberation Sans", "DejaVu Sans"))
fig = plt.figure(figsize=(11, 9))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_aspect("equal"); ax.set_xlim(-3.3, 3.0); ax.set_ylim(-2.3, 2.3); ax.axis("off")
ax.pie(fam_n.values, radius=0.62, startangle=90, counterclock=False,
       colors=[FAMILY_COL[f] for f in FAMILIES], wedgeprops=dict(edgecolor="white", linewidth=2.5))
ax.pie(cnt.n.values, radius=1.0, startangle=90, counterclock=False, colors=cols,
       wedgeprops=dict(width=0.3, edgecolor="white", linewidth=2.5))

angle, n_hosts = 90.0, len(cnt)
for i, r in cnt.iterrows():
    span = 360 * r.n / N
    th_m = np.deg2rad(angle - span / 2)
    angle -= span
    x0, y0 = 1.03 * np.cos(th_m), 1.03 * np.sin(th_m)
    label = f"{r.rep_host_species}\n{r.n:,} ({r.pct:.2f}%)"
    if i < N_BIG:
        th_l = np.deg2rad(BIG_LABEL_ANG[i])
        x1, y1 = 1.12 * np.cos(th_m), 1.12 * np.sin(th_m)
        x2, y2 = R_LAB * np.cos(th_l), R_LAB * np.sin(th_l)
        ax.plot([x0, x1, x2], [y0, y1, y2], color="black", lw=0.9, zorder=4)
        right = np.cos(th_l) >= 0
        s_l = np.sin(th_l)
        ax.text(x2 + (0.04 if right else -0.04), y2, label, ha="left" if right else "right",
                va="bottom" if s_l > 0.6 else ("top" if s_l < -0.6 else "center"), fontsize=15, linespacing=1.2)
    else:
        y2 = COL_TOP - (n_hosts - 1 - i) * COL_STEP
        ax.plot([x0, COL_X + 0.05], [y0, y2], color="black", lw=0.9, zorder=4)
        ax.text(COL_X - 0.04, y2, label, ha="right", va="center", fontsize=13, linespacing=1.2)

ax.text(0, 0, f"Total Genomes: {N:,}\nHost Families: {len(FAMILIES)}\nHost Species: {n_hosts}",
        ha="center", va="center", fontsize=17, linespacing=1.5, zorder=6,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="black", linewidth=1.2))
save(fig, out, "Fig1a_host_sunburst")
print(cnt[["rep_host_family", "rep_host_species", "n"]].to_string(index=False))
