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
#  SCRIPT  : 4f_fig3a_cazy_phylum_tier.py  -  Fig. 3a, CAZy family carriage by phylum and tier
#  AUTHOR  : QRQ
#  INPUT   : figure_inputs/Fig2/supp/representative_host_range.csv
#            figure_inputs/Fig3/supp/54_cazy_family_presence_propagated.tsv
#            figure_inputs/Fig3/supp/gtdbtk.{bac120,ar53}.summary.tsv
#  OUTPUT  : fig3/figC_a_phylum_tier_cazy.{png,pdf,csv}
# ==========================================================

import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent / "_setup"))
from urgg_paths import urgg_in, urgg_out  # inputs: <URGG_PROJECT>/figure_inputs/..., outputs: <URGG_FIGURES>/figN

from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

BASE = urgg_in("figure_inputs")
SHARING_F = BASE / "Fig2/supp/representative_host_range.csv"
PRES_F = BASE / "Fig3/supp/54_cazy_family_presence_propagated.tsv"
GTDB_BAC_F = BASE / "Fig3/supp/gtdbtk.bac120.summary.tsv"
GTDB_AR_F = BASE / "Fig3/supp/gtdbtk.ar53.summary.tsv"
OUT_DIR = urgg_out("fig3")

MIN_PHYLUM_N = 300   # phyla with fewer genomes -> "Other phyla"
MIN_CELL_N = 10      # phylum x tier cells with fewer genomes are left blank

TIERS = ["Core/broad", "Shared", "Host-specific"]
TIER_COL = {"Core/broad": "#B06040", "Shared": "#C89040", "Host-specific": "#2898B8"}  # sampled from reference image

FAMILIES = [  # (group, family), in display order
    ("Strongest core", "GT19"), ("Strongest core", "CE11"),
    ("Strongest core", "GT3"), ("Strongest core", "GT117"),
    ("Strongest host-specific", "CE5"), ("Strongest host-specific", "GH112"),
    ("Strongest host-specific", "GH190"),
    ("Within lineages", "CBM98"), ("Within lineages", "CE2"),
    ("Within lineages", "GH123"), ("Within lineages", "GH1"),
    ("Within lineages", "GH38"), ("Within lineages", "GT30"),
]
GROUPS = list(dict.fromkeys(g for g, _ in FAMILIES))
FAM_ORDER = [f for _, f in FAMILIES]
FAM_GROUP = {fam: grp for grp, fam in FAMILIES}

CMAP = LinearSegmentedColormap.from_list("carry", ["#F3F8FC", "#BBD8EA", "#6AA6C9", "#1F5F87"])

# Times New Roman first; fall back to metric-compatible or serif fonts if absent
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})
FS_LABEL = 10     # phylum labels, legend, colour bar
FS_FAM = 10       # bold family names
FS_GROUP = 10     # group labels
FS_CELL = 7       # percentages inside cells (larger values no longer fit the cell width)

def tier_of(cat: str):
    if cat.startswith("Core"):
        return "Core/broad"
    if cat.startswith("Shared"):
        return "Shared"
    if cat.startswith("Host"):
        return "Host-specific"
    return None

def read_gtdb(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=["user_genome", "classification"], dtype=str)

# ---- tier per representative -------------------------------------------------
sh = pd.read_csv(SHARING_F, usecols=["representative", "sharing_category"], dtype=str)
sh["tier"] = sh["sharing_category"].map(tier_of)
if sh["tier"].isna().any():
    raise ValueError("Unrecognised sharing_category values present")
print("Tier counts:\n", sh["tier"].value_counts())

# ---- phylum per genome ---------------------------------------------------------
tax = pd.concat([read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F)], ignore_index=True)
phylum = tax["classification"].fillna("").str.extract(r"p__([^;]*)", expand=False).fillna("")
tax["phylum"] = phylum.replace("", "Unclassified")
g = sh.merge(tax[["user_genome", "phylum"]], left_on="representative",
             right_on="user_genome", how="left")
n_miss = g["phylum"].isna().sum()
print(f"Representatives: {len(g)} | without GTDB-Tk phylum: {n_miss} ({100 * n_miss / len(g):.1f}%)")
if n_miss / len(g) > 0.5:
    raise ValueError("Over half of the genomes have no phylum: check genome ID format.")
g["phylum"] = g["phylum"].fillna("Unclassified")

counts = g.loc[g["phylum"] != "Unclassified", "phylum"].value_counts()
big = [p for p in counts.index if counts[p] >= MIN_PHYLUM_N]
g["phylum_show"] = np.where(g["phylum"].isin(big), g["phylum"], "Other phyla")
phy_levels = big + ["Other phyla"]
print("Phyla shown:", big)

# ---- carriage per family x phylum x tier ---------------------------------------
pres = pd.read_csv(PRES_F, sep="\t", usecols=["genome", "cazy_family"], dtype=str)
pres = pres[pres["cazy_family"].isin(FAM_ORDER)].drop_duplicates()
miss = set(FAM_ORDER) - set(pres["cazy_family"])
if miss:
    print("WARNING: no carriers found for:", sorted(miss))

den = g.groupby(["phylum_show", "tier"]).size().rename("n").reset_index()

rows = []
for fam in FAM_ORDER:
    carriers = pres.loc[pres["cazy_family"] == fam, "genome"].unique()
    car = g.assign(carrier=g["representative"].isin(carriers))
    agg = car.groupby(["phylum_show", "tier"]).agg(carriers=("carrier", "sum"),
                                                    n=("carrier", "size")).reset_index()
    agg["family"] = fam
    rows.append(agg)
res = pd.concat(rows, ignore_index=True)
res["pct"] = 100 * res["carriers"] / res["n"]
res["group"] = res["family"].map(FAM_GROUP)
res.to_csv(OUT_DIR / "figC_a_phylum_tier_cazy.csv", index=False)

# ---- layout: one column per (phylum, tier), gaps between phyla ----------------
cols = []  # (phylum, tier) or None for a gap
for i, p in enumerate(phy_levels):
    if i > 0:
        cols.append(None)
    for t in TIERS:
        cols.append((p, t))
n_cols = len(cols)

row_items = []  # family name or None for a gap between groups
for i, grp in enumerate(GROUPS):
    if i > 0:
        row_items.append(None)
    row_items.extend(f for g_, f in FAMILIES if g_ == grp)
row_items.append(None)  # gap row under the last group for its label
n_rows = len(row_items)

mat = np.full((n_rows, n_cols), np.nan)
lab = np.full((n_rows, n_cols), "", dtype=object)
res_idx = res.set_index(["family", "phylum_show", "tier"])
for r, fam in enumerate(row_items):
    if fam is None:
        continue
    for c, col in enumerate(cols):
        if col is None or (fam, col[0], col[1]) not in res_idx.index:
            continue
        row = res_idx.loc[(fam, col[0], col[1])]
        if row["n"] >= MIN_CELL_N:
            mat[r, c] = row["pct"]
            lab[r, c] = f"{row['pct']:.0f}"

# x positions: real columns are 1 unit wide, the gaps between phyla are narrower,
# so the heatmap cells get more room for their numbers
GAP_W = 0.35
XC = np.zeros(n_cols)
x_cursor = 0.0
for c_idx, col in enumerate(cols):
    if col is None:
        x_cursor += GAP_W
    else:
        XC[c_idx] = x_cursor + 0.5
        x_cursor += 1.0
XW = x_cursor

# per-phylum tier counts for the pies and phylum totals
tier_n = {p: [int(den.loc[(den["phylum_show"] == p) & (den["tier"] == t), "n"].sum())
              for t in TIERS] for p in phy_levels}

# ---- figure ------------------------------------------------------------------------
# Sized for a double-column journal figure (183 mm wide)
FIG_W, FIG_H = 7.2, 7.8
fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")
LEFT, RIGHT, TOP, BOTTOM = 0.165, 0.815, 0.905, 0.14
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 4.6], hspace=0.04,
                      left=LEFT, right=RIGHT, top=TOP, bottom=BOTTOM)
ax_n = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1])

# top strip: phylum name (staggered so long names do not collide), and a
# genome-count box per tier column with the full integer, written vertically
ax_n.set_xlim(0, XW)
ax_n.set_ylim(0, 1)
ax_n.set_axis_off()

for i, p in enumerate(phy_levels):
    cs = [c for c, col in enumerate(cols) if col is not None and col[0] == p]
    # angled names, all on one baseline above the count boxes
    ax_n.text(float(np.mean(XC[cs])) - 0.2, 0.52, p.replace("_", " "), ha="left", va="bottom",
              rotation=45, rotation_mode="anchor", fontsize=FS_LABEL)
for c, col in enumerate(cols):
    if col is None:
        continue
    n_tier = tier_n[col[0]][TIERS.index(col[1])]
    ax_n.add_patch(plt.Rectangle((XC[c] - 0.46, 0.02), 0.92, 0.42, facecolor=TIER_COL[col[1]],
                                 edgecolor="white", linewidth=0.8))
    ax_n.text(XC[c], 0.23, f"{n_tier:,}", ha="center", va="center", rotation=90,
              fontsize=FS_CELL, color="white", fontweight="bold")
ax_n.text(-0.1, 0.23, "Genomes (n)", ha="right", va="center", fontsize=FS_LABEL)

# heatmap (rows = families, top to bottom)
ax.set_xlim(0, XW)
ax.set_ylim(n_rows - 0.5, -0.5)
for r in range(n_rows):
    for c in range(n_cols):
        if row_items[r] is None or cols[c] is None:
            continue
        v = mat[r, c]
        # missing (too few genomes, or no genomes in that phylum x tier) -> grey hatched "n.a." tile
        is_na = np.isnan(v)
        face = "#EFEFEF" if is_na else CMAP(Normalize(0, 100)(v))
        ax.add_patch(plt.Rectangle((XC[c] - 0.48, r - 0.48), 0.96, 0.96, facecolor=face,
                                   edgecolor="#9A9A9A" if is_na else "white",
                                   linewidth=0.5, hatch="///" if is_na else None))
        if is_na:
            ax.text(XC[c], r, "n.a.", ha="center", va="center", fontsize=FS_CELL, color="#444444")
        elif lab[r, c]:
            txt_col = "white" if v > 70 else "#1A1A1A"
            ax.text(XC[c], r, lab[r, c], ha="center", va="center", fontsize=FS_CELL, color=txt_col)

for r, fam in enumerate(row_items):
    if fam is None:
        continue
    ax.text(-0.5, r, fam, ha="right", va="center", fontsize=FS_FAM, fontweight="bold")
for grp in GROUPS:
    rs = [r for r, f in enumerate(row_items) if f is not None and FAM_GROUP[f] == grp]
    # group name sits in the gap row beneath its block, so the right margin stays free for the legend
    ax.text(XW / 2, max(rs) + 1.0, grp, ha="center", va="center", fontsize=FS_GROUP)
ax.set_xticks([])
ax.set_yticks([])
ax.spines[:].set_visible(False)

# tier key
handles = [Patch(facecolor=TIER_COL[t], edgecolor="none", label=t) for t in TIERS]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(RIGHT + 0.01, 0.80), ncol=1, frameon=False,
           fontsize=FS_LABEL, handlelength=1.2, handleheight=1.0, labelspacing=0.6,
           title="Sharing tier\n(left to right in\neach phylum)",
           title_fontsize=FS_LABEL, alignment="left")

# colour bar for carriage %
cax = fig.add_axes([0.33, 0.075, 0.34, 0.012])
sm = mpl.cm.ScalarMappable(norm=Normalize(0, 100), cmap=CMAP)
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("Genomes carrying the CAZy family (%)", fontsize=FS_LABEL, labelpad=3)
cb.set_ticks([0, 25, 50, 75, 100])
cb.ax.tick_params(labelsize=FS_LABEL - 0.5, length=2)
cb.outline.set_visible(False)

out_png = OUT_DIR / "figC_a_phylum_tier_cazy.png"
out_pdf = OUT_DIR / "figC_a_phylum_tier_cazy.pdf"
fig.savefig(out_png, dpi=600, facecolor="white")
fig.savefig(out_pdf, facecolor="white")
print("Written:", out_png)
