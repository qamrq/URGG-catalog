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
#  SCRIPT  : 5f_fig3b_subpathway_phylum.py  -  Fig. 3b, metabolic sub-pathways: phylum composition and odds ratios
#  AUTHOR  : QRQ
#  INPUT   : figure_inputs/Fig3/supp/subpathway_presence_propagated_50.tsv
#            figure_inputs/Fig2/supp/representative_host_range.csv
#            figure_inputs/Fig3/supp/gtdbtk.{bac120,ar53}.summary.tsv
#  OUTPUT  : fig3/figC_b_subpathway_phylum.{png,pdf}, _composition.csv
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

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Wedge

BASE = urgg_in("figure_inputs")
SHARING_F = BASE / "Fig2/supp/representative_host_range.csv"
PATH_F = BASE / "Fig3/supp/subpathway_presence_propagated_50.tsv"
GTDB_BAC_F = BASE / "Fig3/supp/gtdbtk.bac120.summary.tsv"
GTDB_AR_F = BASE / "Fig3/supp/gtdbtk.ar53.summary.tsv"
OUT_DIR = urgg_out("fig3")

N_PHYLA = 6
X_MAX = 10.0

# tier colours and heatmap-style neutrals, as in Fig. 3a
TIER_COL = {"Core": "#B06040", "Shared": "#C89040", "Host-specific": "#2898B8"}
INK = "#1A1A1A"
NS = "#8C8C8C"        # not significant
NS_FILL = "#E9E6E1"   # "ns" tile
# phylum colours avoid the three tier colours so a pie slice is never read as a tier
PHYLUM_COL = ["#1F5F87", "#6AA6C9", "#3B3B98", "#4FA36B", "#8E5BA8", "#E0B97A"]
OTHER_COL = "#D0D0D0"

MAP = {"butyrate_kinase_route": "Butyrate (kinase)", "propionate_succinate": "Propionate (succinate)",
       "butyrate_transferase_route": "Butyrate (CoA-transferase)", "acetate_production": "Acetate",
       "cobalamin_biosynthesis": "Cobalamin", "urea_utilization": "Urease",
       "ammonia_assimilation": "Ammonia assimilation", "nitrate_reduction": "Nitrate reduction",
       "bcaa_biosynthesis": "Branched-chain amino acids", "lysine_biosynthesis": "Lysine"}

# ---- forest values, from fig_C_template.R (panel b) ------------------------------------
# Same models as Supplementary Table 6 (>=50% markers; adjusted p here is BH over a wider test family,
# so p values differ slightly from ST6; every significance call is identical).
sp = pd.DataFrame({
    "label": ["Butyrate (kinase)", "Propionate (succinate)", "Butyrate (CoA-transferase)", "Acetate", "Cobalamin",
              "Urease", "Ammonia assimilation", "Nitrate reduction", "Branched-chain amino acids", "Lysine"],
    "markers": [3, 5, 3, 2, 4, 3, 3, 2, 2, 2],
    "Core": [32.4, 1.32, 9.25, 89.9, 17.4, 6.61, 74.2, 1.76, 57.7, 88.6],
    "Shared": [21.9, 0.69, 6.80, 91.8, 14.8, 2.60, 69.5, 2.85, 50.7, 88.7],
    "Host-specific": [15.8, 0.54, 5.68, 91.4, 13.9, 1.40, 66.9, 2.81, 46.4, 88.1],
    "or_u": [2.55, 2.47, 1.69, 0.83, 1.30, 4.99, 1.43, 0.62, 1.58, 1.05],
    "lo_u": [2.08, 0.88, 1.20, 0.61, 1.01, 3.27, 1.15, 0.27, 1.30, 0.78],
    "hi_u": [3.13, 5.59, 2.34, 1.16, 1.67, 7.39, 1.78, 1.24, 1.91, 1.43],
    "p_u": [1.1e-16, 0.079, 6.1e-3, 0.35, 0.079, 7.8e-11, 1.9e-3, 0.35, 8.7e-6, 1.0],
    "or_w": [1.24, 5.65, 1.86, 0.82, 1.10, 1.26, 0.73, 0.69, 1.15, 0.39],
    "lo_w": [0.85, 0.56, 0.90, 0.46, 0.73, 0.70, 0.53, 0.25, 0.82, 0.21],
    "hi_w": [1.81, 57.2, 3.83, 1.49, 1.66, 2.28, 1.01, 1.90, 1.63, 0.69],
    "p_w": [0.55, 0.55, 0.55, 0.73, 0.85, 0.73, 0.017, 0.55, 0.73, 0.040],
    "or_75": [1.32, np.nan, 2.09, 1.41, 1.34, 5.87, 1.51, 3.71, 1.76, 1.13],
    "call_75": ["ns", "nt", "core", "core", "ns", "core", "core", "ns", "core", "ns"],
})
sp["or_100"] = np.where(sp.markers <= 3, sp.or_75, np.nan)
sp["call_100"] = np.where(sp.markers <= 3, sp.call_75, "nt")
sp["dagger"] = sp.label == "Lysine"

def call_rule(lo, hi, p):
    """significant = BH p < 0.05 and the CI excludes 1 (same rule as the R template)"""
    return (p < 0.05) & ((lo > 1) | (hi < 1))

sp["sig_u"] = call_rule(sp.lo_u, sp.hi_u, sp.p_u)
sp["sig_w"] = call_rule(sp.lo_w, sp.hi_w, sp.p_w)
assert sp.sig_u.sum() == 5 and list(sp.label[sp.sig_w]) == ["Lysine"], "forest rule check failed"

# group gaps: SCFA | cofactor | nitrogen | amino acids
GAP = [1.5] * 4 + [1.0] + [0.5] * 3 + [0.0] * 2
sp["y"] = np.arange(len(sp), 0, -1) + np.array(GAP)
TOP_Y = sp.y.max()
REF_Y = TOP_Y + 2.0                     # "All genomes" reference row (phylum pie only)

# ---- phylum of every representative -----------------------------------------------------
def read_gtdb(path):
    return pd.read_csv(path, sep="\t", usecols=["user_genome", "classification"], dtype=str)

sh = pd.read_csv(SHARING_F, usecols=["representative"], dtype=str)
tax = pd.concat([read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F)], ignore_index=True)
tax["phylum"] = (tax["classification"].fillna("").str.extract(r"p__([^;]*)", expand=False)
                 .fillna("").replace("", "Unclassified"))
g = sh.merge(tax[["user_genome", "phylum"]], left_on="representative", right_on="user_genome", how="left")
g["phylum"] = g["phylum"].fillna("Unclassified")
counts = g.loc[g.phylum != "Unclassified", "phylum"].value_counts()
top = list(counts.index[:N_PHYLA])
g["ph"] = np.where(g.phylum.isin(top), g.phylum, "Other phyla")
LV = top + ["Other phyla"]
PHC = dict(zip(top, PHYLUM_COL[:len(top)]))
PHC["Other phyla"] = OTHER_COL
print("Phyla shown:", top)

# ---- carriers of each sub-pathway, by phylum -----------------------------------------
pr = pd.read_csv(PATH_F, sep="\t", dtype=str).iloc[:, :2]
pr.columns = ["genome", "subpathway"]
pr = pr[pr.subpathway.isin(MAP)].drop_duplicates()
pr["label"] = pr.subpathway.map(MAP)
m = pr.merge(g[["representative", "ph"]], left_on="genome", right_on="representative", how="left")
matched = m.ph.notna().mean()
print(f"Carrier rows: {len(m)} | matched to a phylum: {100 * matched:.1f}%")
if matched < 0.95:
    raise ValueError("Under 95% of carrier genomes matched the representatives: check genome ID format.")
m = m.dropna(subset=["ph"])

comp = m.groupby(["label", "ph"]).size().rename("N").reset_index()
comp["share"] = 100 * comp.N / comp.groupby("label").N.transform("sum")
ref = g.groupby("ph").size().rename("N").reset_index()
ref["label"] = "All genomes"
ref["share"] = 100 * ref.N / ref.N.sum()
allc = pd.concat([comp, ref], ignore_index=True)
allc.to_csv(OUT_DIR / "figC_b_subpathway_phylum_composition.csv", index=False)

def shares_of(label):
    d = allc[allc.label == label].set_index("ph").share.reindex(LV).fillna(0.0)
    return d.values

# ---- figure -------------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})
FS = 10          # labels
FS_SM = 8.5      # numbers in the forest
FS_TILE = 8

FIG_W, FIG_H = 7.2, 7.0
fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")
AX_BOTTOM, AX_H = 0.17, 0.70
ax_l = fig.add_axes([0.01, AX_BOTTOM, 0.46, AX_H])
ax_r = fig.add_axes([0.52, AX_BOTTOM, 0.46, AX_H])

Y_LO = sp.y.min() - 0.7
Y_HI = 17.0
ax_l.set_ylim(Y_LO, Y_HI)
ax_r.set_ylim(Y_LO, Y_HI)
ax_l.set_axis_off()

# left: pies (drawn in units where one row = one unit in both directions, so they are circles)
unit_in = AX_H * FIG_H / (Y_HI - Y_LO)
W_UNITS = 0.46 * FIG_W / unit_in
ax_l.set_xlim(0, W_UNITS)
BAR_W = 1.45 / unit_in            # stacked 100% bar, right-hand part of the left panel
BAR_X0 = W_UNITS - BAR_W - 0.1
BAR_H = 0.62
DARK_PHYLA = {"Bacillota", "Bacillota_I", "Verrucomicrobiota"}   # white numbers on these

def break_name(s, limit=14):
    """split a long sub-pathway name onto two lines at the space that balances the two halves"""
    if len(s) <= limit or " " not in s:
        return s
    words = s.split(" ")
    best = None
    for i in range(1, len(words)):
        a, b = " ".join(words[:i]), " ".join(words[i:])
        if best is None or max(len(a), len(b)) < best[0]:
            best = (max(len(a), len(b)), a + "\n" + b)
    return best[1]

rows_l = [(lab, y, False) for lab, y in zip(sp.label, sp.y)] + [("All genomes", REF_Y, True)]
for lab, y, is_ref in rows_l:
    text = break_name(lab + (" †" if (lab == "Lysine") else ""))
    # straight labels, right-aligned just left of their bar; long names break onto two lines
    ax_l.text(BAR_X0 - 0.35, y, text, ha="right", va="center", multialignment="right", fontsize=FS,
              color="#555555" if is_ref else INK, linespacing=1.0)
    vals = shares_of(lab)
    if is_ref:   # the all-genome reference is a pie, the sub-pathways are bars
        cx, R = BAR_X0 + BAR_W / 2, 1.1
        start = 90.0
        for i, v in enumerate(vals):
            if v <= 0:
                continue
            sweep = 360.0 * v / 100.0
            ax_l.add_patch(Wedge((cx, y), R, start - sweep, start, facecolor=PHC[LV[i]],
                                 edgecolor="white", linewidth=0.4))
            if v >= 8:
                mid = np.deg2rad(start - sweep / 2)
                ax_l.text(cx + 0.6 * R * np.cos(mid), y + 0.6 * R * np.sin(mid), f"{v:.0f}",
                          ha="center", va="center", fontsize=FS_TILE,
                          color="white" if LV[i] in DARK_PHYLA else INK)
            start -= sweep
        continue
    x_cur = BAR_X0
    for i, v in enumerate(vals):
        if v <= 0:
            continue
        w = BAR_W * v / 100.0
        ax_l.add_patch(Rectangle((x_cur, y - BAR_H / 2), w, BAR_H, facecolor=PHC[LV[i]],
                                 edgecolor="white", linewidth=0.4))
        if v >= 8:   # print the share only where the segment is wide enough to hold it
            ax_l.text(x_cur + w / 2, y, f"{v:.0f}", ha="center", va="center", fontsize=FS_TILE,
                      color="white" if LV[i] in DARK_PHYLA else INK)
        x_cur += w

# right: forest
X_LAB = {"Core": 10 ** -2.55, "Shared": 10 ** -1.92, "Host-specific": 10 ** -1.30}
ax_r.set_xscale("log")
ax_r.set_xlim(10 ** -2.66, 10 ** 2.12)
ax_r.axvline(1, color="#8A8A8A", linestyle=(0, (2, 2)), linewidth=0.6, zorder=0)
for v in [0.2, 0.5, 2, 5]:
    ax_r.axvline(v, color="#E6E6E6", linewidth=0.5, zorder=0)
ax_r.set_xticks([0.2, 0.5, 1, 2, 5, 10])
ax_r.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"], fontsize=FS - 1, color=INK)
ax_r.minorticks_off()
ax_r.set_yticks([])
for s in ["left", "right", "top"]:
    ax_r.spines[s].set_visible(False)
ax_r.spines["bottom"].set_linewidth(0.6)
ax_r.set_xlabel("Odds ratio, core vs host-specific (log scale)", fontsize=FS - 1, color=INK)

# prevalence columns
for tier, x in X_LAB.items():
    ax_r.text(x, Y_HI - 0.6, tier, ha="left", va="bottom",
              fontsize=FS - 2, fontweight="bold", color=TIER_COL[tier], rotation=45, rotation_mode="anchor")
for _, r in sp.iterrows():
    for tier, x in X_LAB.items():
        ax_r.text(x, r.y, f"{r[tier]:.1f}", ha="center", va="center", fontsize=FS_SM,
                  color=INK if not False else INK)

# OR rows: open circle = unadjusted, filled = within lineages
for _, r in sp.iterrows():
    for model, y_off, filled in [("u", 0.17, False), ("w", -0.17, True)]:
        orv, lo, hi = r[f"or_{model}"], r[f"lo_{model}"], r[f"hi_{model}"]
        sig = r[f"sig_{model}"]
        col = INK if sig else NS
        y = r.y + y_off
        hi_d = min(hi, X_MAX)
        ax_r.plot([lo, hi_d], [y, y], color=col, linewidth=0.9, solid_capstyle="butt", zorder=2)
        if hi > X_MAX:
            ax_r.annotate("", xy=(X_MAX, y), xytext=(X_MAX / 1.25, y),
                          arrowprops=dict(arrowstyle="-|>", color=col, lw=0.8, mutation_scale=6), zorder=2)
            ax_r.text(X_MAX, y + (0.33 if model == "u" else -0.33), f"to {hi:g}", ha="right", va="center",
                      fontsize=FS_TILE - 0.5, color=NS)
        ax_r.plot(orv, y, marker="o", markersize=4.6 if not filled else 4.8,
                  markerfacecolor="white" if not filled else col, markeredgecolor=col,
                  markeredgewidth=0.9, zorder=3)

# threshold tiles (>=75% and 100% markers), in log space
TILE_X = {"75": 10 ** 1.36, "100": 10 ** 1.90}
for _, r in sp.iterrows():
    for key, col_key, or_key in [("75", "call_75", "or_75"), ("100", "call_100", "or_100")]:
        call = r[col_key]
        cx = TILE_X[key]
        left, right = 10 ** (np.log10(cx) - 0.15), 10 ** (np.log10(cx) + 0.15)
        if call == "core":
            fc, tc = TIER_COL["Core"], "white"
        elif call == "ns":
            fc, tc = NS_FILL, INK
        else:
            fc, tc = "white", NS
        ax_r.add_patch(Rectangle((left, r.y - 0.41), right - left, 0.82, facecolor=fc,
                                 edgecolor="#C9C9C9", linewidth=0.4, zorder=1))
        txt = "n.t." if pd.isna(r[or_key]) else f"{r[or_key]:.2f}"
        ax_r.text(cx, r.y, txt, ha="center", va="center", fontsize=FS_TILE, color=tc, zorder=4)
for key, cx in TILE_X.items():
    ax_r.text(cx, Y_HI - 0.6, "≥75%" if key == "75" else "100%", ha="center", va="center",
              fontsize=FS - 2, fontweight="bold", color=INK)

# legend
handles = [Patch(facecolor=PHC[p], edgecolor="none", label=p) for p in LV]
handles += [
    Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=INK, markersize=5,
           label="Unadjusted"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markeredgecolor=INK, markersize=5,
           label="Within lineages"),
    Line2D([0], [0], color=INK, linewidth=1, label="BH p < 0.05, CI excludes 1"),
    Line2D([0], [0], color=NS, linewidth=1, label="Not significant"),
]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9.5,
           handlelength=1.2, columnspacing=0.9, handletextpad=0.5, bbox_to_anchor=(0.5, 0.0))
fig.text(0.01, 0.955, "Phylum composition of carrying genomes", fontsize=FS - 0.5, ha="left", color=INK)
fig.text(0.52, 0.955, "Prevalence (%) and odds ratios", fontsize=FS - 0.5, ha="left", color=INK)

out_png = OUT_DIR / "figC_b_subpathway_phylum.png"
fig.savefig(out_png, dpi=600, facecolor="white")
fig.savefig(OUT_DIR / "figC_b_subpathway_phylum.pdf", facecolor="white")
print("Written:", out_png)
