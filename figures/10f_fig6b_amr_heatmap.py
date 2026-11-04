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
#  SCRIPT  : 10f_fig6b_amr_heatmap.py  -  Fig. 6b, AMR gene families across host species
#  AUTHOR  : QRQ
#  INPUT   : figure_inputs/Fig6/result_files/amr_heatmap_matrix_reps.tsv
#            figure_inputs/Fig6/result_files/amr_family_meta.tsv
#  OUTPUT  : fig6/figure6b_amr_heatmap.{png,pdf}
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

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.colors import to_hex, to_rgb
from matplotlib.ticker import MaxNLocator
from scipy.cluster.hierarchy import linkage, dendrogram

TABLE_DIR = urgg_in("figure_inputs/Fig6/result_files")
OUT_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else urgg_out("fig6")
DPI = int(sys.argv[2]) if len(sys.argv) > 2 else 600
OUT_STEM = "figure6b_amr_heatmap"
MAG_SPLIT = 1000

FONT_CANDIDATES = ["Times New Roman", "Liberation Serif", "Nimbus Roman", "Tinos", "DejaVu Serif"]
_available = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in FONT_CANDIDATES if f in _available), "serif")
plt.rcParams["font.family"] = FONT
print("Font used:", FONT)

FS_ROW = 24
FS_COL = 24
FS_AXIS = 27
FS_STRIP = 21
FS_LEG_T = 24
FS_LEG = 21
FS_TICK = 20

mat = pd.read_csv(TABLE_DIR / "amr_heatmap_matrix_reps.tsv", sep="\t")
meta = pd.read_csv(TABLE_DIR / "amr_family_meta.tsv", sep="\t")

RENAME = {
    "AAC(6')-I family aminoglycoside 6'-N-acetyltransferase": "AAC(6')-I acetyltransferase",
    "aminoglycoside 6'-N-acetyltransferase": "AAC(6')-N-acetyltransferase",
    "ANT(9) family aminoglycoside nucleotidyltransferase": "ANT(9) nucleotidyltransferase",
    "CfxA family broad-spectrum class A beta-lactamase": "CfxA family beta-lactamase",
    "tetracycline resistance ribosomal protection protein": "tetracycline ribosomal protection protein",
    "23S rRNA (adenine(2058)-N(6))-methyltransferase": "23S rRNA methyltransferase",
    "trimethoprim-resistant dihydrofolate reductase": "trimethoprim-resistant DHFR",
}
shorten = lambda x: RENAME.get(x, x)
mat["fam"] = mat["family"].map(shorten)
meta["fam"] = meta["family"].map(shorten)

breadth = (mat[mat["n_hits"] > 0].groupby("fam")["host"].nunique()
           .rename("n_species"))
meta = meta.merge(breadth, left_on="fam", right_index=True, how="left")
meta["n_species"] = meta["n_species"].fillna(0).astype(int)
meta["tier"] = np.select([meta["n_species"] <= 1, meta["n_species"] <= 4],
                         ["Host-specific", "Shared"], default="Core/broad")
meta = meta.set_index("fam")

sentence = lambda s: s.str.lower().str.capitalize()
meta["class_lab"] = sentence(meta["class_"])
meta["mech_lab"] = sentence(meta["mechanism"])

BINOM = {
    "Ovis_aries": "Ovis aries", "Bos_taurus": "Bos taurus", "Bos_grunniens": "Bos grunniens",
    "Pantholops_hodgsonii": "Pantholops hodgsonii", "Capra_hircus": "Capra hircus",
    "Bubalus_bubalis": "Bubalus bubalis", "Bos_indicus": "Bos indicus",
    "Hydropotes_inermis": "Hydropotes inermis", "Alces_alces": "Alces alces",
    "Procapra_picticaudata": "Procapra picticaudata", "Cervus_elaphus": "Cervus elaphus",
    "Giraffa_camelopardalis": "Giraffa camelopardalis", "Rangifer_tarandus": "Rangifer tarandus",
    "Alces_americanus": "Alces americanus", "Camelus_bactrianus": "Camelus bactrianus",
    "Camelus_dromedarius": "Camelus dromedarius", "Bos_frontalis": "Bos frontalis",
}

hosts_all = mat[["host", "total_mags"]].drop_duplicates()
hosts_sel = hosts_all.loc[hosts_all["total_mags"] >= MAG_SPLIT, "host"].tolist()
sub = mat[mat["host"].isin(hosts_sel)]

# ---- clustering (on log rate; zero cells = log10(0.01)) ----
fam_lv = sorted(sub["fam"].unique())
W = pd.DataFrame(np.log10(0.01), index=fam_lv, columns=hosts_sel)
for f, h, r in zip(sub["fam"], sub["host"], sub["rate_per_1000"]):
    W.loc[f, h] = np.log10(r + 0.01)

Z_row = linkage(W.values, method="average", metric="euclidean")
Z_col = linkage(W.values.T, method="average", metric="euclidean")
row_dn = dendrogram(Z_row, no_plot=True)
col_dn = dendrogram(Z_col, no_plot=True)
row_order = [fam_lv[i] for i in row_dn["leaves"]]
col_order = [hosts_sel[i] for i in col_dn["leaves"]]
nR, nC = len(row_order), len(col_order)

CLASS_PAL = {"Tetracycline": "#D85A30", "Aminoglycoside": "#0F6E56", "Lincosamide": "#93A8D6",
             "Lincosamide/macrolide/streptogramin": "#3C3489", "Macrolide": "#F0997B",
             "Beta-lactam": "#854F0B", "Nitroimidazole": "#26215C", "Trimethoprim": "#D6923A"}
MECH_PAL = {"Antibiotic efflux": "#1F2A5C", "Antibiotic inactivation": "#D85A30",
            "Antibiotic target alteration": "#5FA83C", "Antibiotic target protection": "#4FA8B4"}
TIER_PAL = {"Host-specific": "#2898B8", "Shared": "#C89040", "Core/broad": "#B06040"}
TIER_LABEL = {"Host-specific": "Host specific", "Core/broad": "Core"}
ZERO_COL = "#B4B2A9"

def ramp_colours(stops, n=100):
    """Linear RGB interpolation, equivalent to R colorRampPalette(stops)(n)."""
    rgb = np.array([to_rgb(s) for s in stops])
    t = np.linspace(0, 1, n)
    seg = np.linspace(0, 1, len(stops))
    return [to_hex(c) for c in np.stack([np.interp(t, seg, rgb[:, k]) for k in range(3)], axis=1)]

RAMP = ramp_colours(["#3B2D96", "#F1EFE8", "#CB519B"])

CAN_W, CAN_H = 2300, 1520
x0, cw = 760, 112
y0, rh = 215, 46
x1, y1 = x0 + nC * cw, y0 + nR * rh
strip_x, strip_w = [24, 60, 96], 30
d_top, d_bot = 60, y0 - 22
lx = 1740

N_HITS = sub.pivot_table(index="fam", columns="host", values="n_hits", aggfunc="sum") \
    .reindex(index=row_order, columns=col_order).fillna(0).values
RATE = sub.pivot_table(index="fam", columns="host", values="rate_per_1000", aggfunc="sum") \
    .reindex(index=row_order, columns=col_order).fillna(0).values
has = N_HITS > 0
logr = np.log10(RATE[has])
r0, r1 = logr.min(), logr.max()
if r1 == r0:
    r0, r1 = r0 - 0.1, r1 + 0.1
idx = np.clip(np.round(1 + 99 * (np.log10(np.maximum(RATE, 1e-9)) - r0) / (r1 - r0)), 1, 100).astype(int)
CELL_COL = np.where(has, np.array(RAMP)[idx - 1], ZERO_COL)

fig = plt.figure(figsize=(CAN_W / 100, CAN_H / 100), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, CAN_W)
ax.set_ylim(CAN_H, 0)
ax.axis("off")

hmax = Z_col[:, 2].max()
for ic, dc in zip(col_dn["icoord"], col_dn["dcoord"]):
    X = x0 + (np.array(ic) - 5) / 10 * cw + cw / 2
    Y = d_bot - np.array(dc) / hmax * (d_bot - d_top)
    ax.plot(X, Y, color="black", lw=2.5, solid_capstyle="butt")

root_k = max(range(len(col_dn["dcoord"])), key=lambda k: col_dn["dcoord"][k][1])
for arm, name, side in [(0, "Clade I", -1), (3, "Clade II", 1)]:
    X = x0 + (col_dn["icoord"][root_k][arm] - 5) / 10 * cw + cw / 2
    Y = d_bot - col_dn["dcoord"][root_k][arm] / hmax * (d_bot - d_top)
    ax.text(X + side * 8, Y - 14, name, ha="right" if side < 0 else "left", va="bottom",
            fontsize=FS_STRIP, fontweight="bold")

for i in range(nR):
    for j in range(nC):
        ax.add_patch(Rectangle((x0 + j * cw, y0 + i * rh), cw, rh, facecolor=CELL_COL[i, j],
                               edgecolor="white", lw=2.5))

mrow = meta.loc[row_order]
for k, (col_map, key) in enumerate([(CLASS_PAL, "class_lab"), (MECH_PAL, "mech_lab"), (TIER_PAL, "tier")]):
    for i, fam in enumerate(row_order):
        ax.add_patch(Rectangle((strip_x[k], y0 + i * rh), strip_w, rh,
                               facecolor=col_map[mrow.loc[fam, key]], edgecolor="white", lw=1.7))
for x, lab in zip(strip_x, ["Drug class", "Mechanism", "Sharing tier"]):
    ax.text(x + strip_w / 2, y0 - 10, lab, rotation=90, ha="left", va="center",
            rotation_mode="anchor", fontsize=FS_STRIP)

for i, fam in enumerate(row_order):
    ax.text(x0 - 14, y0 + (i + 0.5) * rh, fam, ha="right", va="center", fontsize=FS_ROW)

col_label_artists = []
for j, host in enumerate(col_order):
    t = ax.text(x0 + (j + 0.5) * cw + 14, y1 + 16, BINOM[host], rotation=45, ha="right", va="top",
                rotation_mode="anchor", style="italic", fontsize=FS_COL)
    col_label_artists.append(t)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
lowest = max(ax.transData.inverted().transform(t.get_window_extent(renderer))[:, 1].max()
             for t in col_label_artists)
ax.text((x0 + x1) / 2, lowest + 45, "Host species", ha="center", va="center", fontsize=FS_AXIS)

def legend_block(y, title, labels, colours, step=44):
    ax.text(lx, y, title, ha="left", va="center", fontsize=FS_LEG_T, fontweight="bold")
    for i, (lab, col) in enumerate(zip(labels, colours)):
        yy = y + 46 + i * step
        ax.add_patch(Rectangle((lx, yy - 17), 34, 34, facecolor=col, edgecolor="grey", lw=0.8))
        ax.text(lx + 50, yy, lab, ha="left", va="center", fontsize=FS_LEG)
    return y + 46 + (len(labels) - 1) * step + 17

cls_used = sorted(meta["class_lab"].unique())
mech_used = sorted(meta["mech_lab"].unique())
tier_used = [t for t in TIER_PAL if t in set(meta["tier"])]

y_cls = y0 + 8
bottom = legend_block(y_cls, "Drug class", cls_used, [CLASS_PAL[c] for c in cls_used])
y_mech = bottom + 50
bottom = legend_block(y_mech, "Resistance mechanism", mech_used, [MECH_PAL[m] for m in mech_used])
y_tier = bottom + 50
bottom = legend_block(y_tier, "Sharing tier", [TIER_LABEL.get(t, t) for t in tier_used],
                      [TIER_PAL[t] for t in tier_used])

bar_w, bar_h = 340, 26
y_bar = bottom + 110
bar_top = y_bar + 70
bar_bottom = bar_top + bar_h
ax.text(lx, y_bar, "log10(hits per 1,000\nrepresentative genomes)", ha="left", va="bottom",
        fontsize=FS_LEG_T, fontweight="bold", linespacing=0.95)
seg_w = bar_w / 100
for k in range(100):
    ax.add_patch(Rectangle((lx + k * seg_w, bar_top), seg_w + 0.5, bar_h, facecolor=RAMP[k], edgecolor="none"))
ax.add_patch(Rectangle((lx, bar_top), bar_w, bar_h, fill=False, edgecolor="grey", lw=0.8))
ticks = [t for t in MaxNLocator(3).tick_values(r0, r1) if r0 <= t <= r1]
for t in ticks:
    tx = lx + (t - r0) / (r1 - r0) * bar_w
    ax.plot([tx, tx], [bar_bottom, bar_bottom + 7], color="black", lw=1.2)
    ax.text(tx, bar_bottom + 26, f"{t:.1f}", ha="center", va="center", fontsize=FS_TICK)

nh_y = bar_bottom + 66
ax.add_patch(Rectangle((lx, nh_y - 14), 34, 34, facecolor=ZERO_COL, edgecolor="grey", lw=0.8))
ax.text(lx + 50, nh_y, "No hits", ha="left", va="center", fontsize=FS_LEG)

OUT_DIR.mkdir(parents=True, exist_ok=True)
png_path = OUT_DIR / f"{OUT_STEM}.png"
fig.savefig(png_path, dpi=DPI, facecolor="white")
fig.savefig(OUT_DIR / f"{OUT_STEM}.pdf", facecolor="white")
print("Hosts shown (%d): %s" % (nC, ", ".join(BINOM[h] for h in col_order)))
print("Figure written to", png_path, "at", DPI, "dpi")
