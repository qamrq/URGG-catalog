"""
Fig. 3c (Python version of 82_figC_phylum_tier_cazy.R)

Phylum x sharing tier x CAZy family. For each CAZy family: % of genomes carrying
the family, split by GTDB phylum and by tier (Core/broad, Shared, Host-specific).

Inputs (read-only)
  representative_host_range.csv        tier per representative genome
  54_cazy_family_presence_propagated.tsv   genome x cazy_family (after propagation)
  gtdbtk.bac120.summary.tsv / gtdbtk.ar53.summary.tsv   phylum per genome

Outputs (written next to this script)
  figC_a_phylum_tier_cazy.png  (600 dpi)
  figC_a_phylum_tier_cazy.pdf
  figC_a_phylum_tier_cazy.csv  one row per family x phylum x tier
"""
from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

BASE = Path(r"D:\QRQ_ZJU\MAG-SAG\Manuscript\Figures")
SHARING_F = BASE / r"Fig2\supp\representative_host_range.csv"
PRES_F = BASE / r"Fig3\supp\54_cazy_family_presence_propagated.tsv"
GTDB_BAC_F = BASE / r"Fig3\supp\gtdbtk.bac120.summary.tsv"
GTDB_AR_F = BASE / r"Fig3\supp\gtdbtk.ar53.summary.tsv"
OUT_DIR = BASE / r"Fig3\scrap"

MIN_PHYLUM_N = 300   # phyla with fewer genomes -> "Other phyla"
MIN_CELL_N = 10      # phylum x tier cells with fewer genomes are left blank

TIERS = ["Core/broad", "Shared", "Host-specific"]
TIER_COL = {"Core/broad": "#C0392B", "Shared": "#E8A33D", "Host-specific": "#1B7F9E"}

FAMILIES = [  # (group, family), in display order
    ("Strongest core", "GT19"), ("Strongest core", "CE11"),
    ("Strongest core", "GT3"), ("Strongest core", "GT117"),
    ("Strongest host-specific", "CE5"), ("Strongest host-specific", "GH112"),
    ("Strongest host-specific", "GH190"),
    ("Within lineages", "CBM98"), ("Within lineages", "CE2"),
    ("Within lineages", "GH123"), ("Within lineages", "GH1"),
    ("Within lineages", "GH38"), ("Within lineages", "GT30"),
    ("Loss:gain outlier", "GH189"),
]
GROUPS = list(dict.fromkeys(g for g, _ in FAMILIES))
FAM_ORDER = [f for _, f in FAMILIES]
FAM_GROUP = {fam: grp for grp, fam in FAMILIES}

CMAP = LinearSegmentedColormap.from_list("carry", ["#F5F3EE", "#8FBFC4", "#1B6B76", "#0B3D45"])

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def tier_of(cat: str):
    if cat.startswith("Core"):
        return "Core/broad"
    if cat.startswith("Shared"):
        return "Shared"
    if cat.startswith("Host"):
        return "Host-specific"
    return None


def read_gtdb(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", usecols=["user_genome", "classification"], dtype=str)
    return df


# ---- tier per representative -------------------------------------------------
sh = pd.read_csv(SHARING_F, usecols=["representative", "sharing_category"], dtype=str)
sh["tier"] = sh["sharing_category"].map(tier_of)
if sh["tier"].isna().any():
    raise ValueError("Unrecognised sharing_category values present")
print("Tier counts:\n", sh["tier"].value_counts())

# ---- phylum per genome ---------------------------------------------------------
tax = pd.concat([read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F)], ignore_index=True)
cls = tax["classification"].fillna("")
phylum = cls.str.extract(r"p__([^;]*)", expand=False).fillna("")
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

cell_of = g.set_index("representative")[["phylum_show", "tier"]]
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
cols = []  # list of (phylum, tier) or None for a gap
for i, p in enumerate(phy_levels):
    if i > 0:
        cols.append(None)
    for t in TIERS:
        cols.append((p, t))
n_cols = len(cols)

# rows: families with a gap row between groups
row_items = []
for i, grp in enumerate(GROUPS):
    if i > 0:
        row_items.append(None)
    for fam in [f for g_, f in FAMILIES if g_ == grp]:
        row_items.append(fam)
n_rows = len(row_items)

mat = np.full((n_rows, n_cols), np.nan)
lab = np.full((n_rows, n_cols), "", dtype=object)
res_idx = res.set_index(["family", "phylum_show", "tier"])
for r, fam in enumerate(row_items):
    if fam is None:
        continue
    for c, col in enumerate(cols):
        if col is None:
            continue
        p, t = col
        if (fam, p, t) not in res_idx.index:
            continue
        row = res_idx.loc[(fam, p, t)]
        if row["n"] >= MIN_CELL_N:
            mat[r, c] = row["pct"]
            lab[r, c] = f"{row['pct']:.0f}"

# totals for the top bar
top_n = np.full(n_cols, np.nan)
for c, col in enumerate(cols):
    if col is None:
        continue
    m = den[(den["phylum_show"] == col[0]) & (den["tier"] == col[1])]
    top_n[c] = m["n"].iloc[0] if len(m) else 0


def fmt_n(v):
    v = int(v)
    return f"{v / 1000:.1f}k" if v >= 1000 else str(v)


# ---- figure ----------------------------------------------------------------------
fig = plt.figure(figsize=(18, 14.5), facecolor="white")
gs = fig.add_gridspec(3, 1, height_ratios=[0.5, 0.45, 9.0], hspace=0.12)
ax_leg = fig.add_subplot(gs[0])
ax_n = fig.add_subplot(gs[1])
ax = fig.add_subplot(gs[2])

# legend (tier colours)
ax_leg.axis("off")
handles = [Patch(facecolor=TIER_COL[t], edgecolor="none", label=t) for t in TIERS]
ax_leg.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=13,
              title="Sharing tier (columns in each phylum, left to right)", title_fontsize=13,
              handlelength=1.6, handleheight=1.2, columnspacing=2.5)

# genome-count row
ax_n.set_xlim(-0.5, n_cols - 0.5)
ax_n.set_ylim(0, 1)
for c, col in enumerate(cols):
    if col is None:
        continue
    ax_n.add_patch(plt.Rectangle((c - 0.45, 0.05), 0.9, 0.9, facecolor=TIER_COL[col[1]],
                                 edgecolor="white", linewidth=1.5))
    ax_n.text(c, 0.5, fmt_n(top_n[c]), ha="center", va="center", color="white",
              fontsize=9.5, fontweight="bold")
ax_n.set_axis_off()
ax_n.text(-0.7, 0.5, "Genomes (n)", ha="right", va="center", fontsize=13)

# heatmap (rows = families, top to bottom)
masked = np.ma.masked_invalid(mat)
ax.set_facecolor("white")
for r in range(n_rows):
    for c in range(n_cols):
        if row_items[r] is None or cols[c] is None:
            continue
        v = mat[r, c]
        face = "#FFFFFF" if np.isnan(v) else CMAP(Normalize(0, 100)(v))
        ax.add_patch(plt.Rectangle((c - 0.45, r - 0.45), 0.9, 0.9, facecolor=face,
                                   edgecolor="white", linewidth=1.5))
        if lab[r, c]:
            txt_col = "white" if (not np.isnan(v) and v > 55) else "#222222"
            ax.text(c, r, lab[r, c], ha="center", va="center", fontsize=12, color=txt_col)
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(n_rows - 0.5, -0.5)

# y labels: family names (bold) and group labels on the right
for r, fam in enumerate(row_items):
    if fam is None:
        continue
    ax.text(-0.75, r, fam, ha="right", va="center", fontsize=14, fontweight="bold", color="#333333")
for grp in GROUPS:
    rs = [r for r, f in enumerate(row_items) if f is not None and FAM_GROUP[f] == grp]
    ax.text(n_cols - 0.3, np.mean(rs), grp, ha="left", va="center", fontsize=13, color="#333333")
ax.set_yticks([])

# x labels: phylum names centred over their three tier columns
ax.set_xticks([])
for p in phy_levels:
    cs = [c for c, col in enumerate(cols) if col is not None and col[0] == p]
    label = textwrap.fill(p.replace("_", " "), 20)
    ax_n.text(np.mean(cs), 1.05, label, ha="center", va="bottom", fontsize=11,
              fontweight="bold", linespacing=0.95, transform=ax_n.transData)

ax.spines[:].set_visible(False)

# colourbar for carriage %
cax = fig.add_axes([0.36, 0.025, 0.28, 0.012])
sm = mpl.cm.ScalarMappable(norm=Normalize(0, 100), cmap=CMAP)
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("% of genomes carrying the CAZy family", fontsize=13, labelpad=6)
cb.set_ticks([0, 25, 50, 75, 100])
cb.ax.tick_params(labelsize=12)
cb.outline.set_visible(False)

fig.subplots_adjust(left=0.12, right=0.84, top=0.95, bottom=0.10)
out_png = OUT_DIR / "figC_a_phylum_tier_cazy.png"
out_pdf = OUT_DIR / "figC_a_phylum_tier_cazy.pdf"
fig.savefig(out_png, dpi=600, facecolor="white")
fig.savefig(out_pdf, facecolor="white")
print("Written:", out_png)
