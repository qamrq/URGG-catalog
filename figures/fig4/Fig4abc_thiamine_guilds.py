"""
Figure A: thiamine producer and salvager guilds (Python port of 67_thiamine_figures.R, Figure A).

a  pathway schematic: de novo route (KEGG M00127, ThiH type) | membrane ring around the thiamine pool |
   salvage routes (thiamine import through a transporter + TPK; HMP/HET salvage, KEGG M00899),
   with representative genera from this catalogue
b  guild composition by host-sharing tier
c  raw (Fisher) vs within-lineage (GLMM) odds ratios, core/broad vs host-specific

Inputs (read-only; override any path with an environment variable of the same name)
  G63_F        63_genome_thiamine_guilds.tsv             one row per genome (guild, sharing_category, traits)
  GLMM_F       table_s49_thiamine_guild_presence_models.csv
  SHARING_F    representative_host_range.csv             representative, host_species_list
  GTDB_BAC_F   gtdbtk.bac120.summary.tsv                 genus for each genome
  GTDB_AR_F    gtdbtk.ar53.summary.tsv
  LINEAGE_F    lineage_full_classification_table.csv     RUC lineage membership

Outputs (OUT_DIR)
  Fig_thiamine_guilds.png (600 dpi), .pdf (vector, Times New Roman embedded), .svg
  panel_c_plot_table.csv                                 the values plotted in panel c, for audit
  figA_genus_thiamine_summary.csv                        genus-level counts used to pick representative genera
"""
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "_setup"))
from urgg_paths import urgg_in, urgg_out  # inputs: <URGG_PROJECT>/figure_inputs/..., outputs: <URGG_FIGURES>/figN

import os
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Patch

# =============================================================================
# PART 1  Paths, palettes, helpers
# =============================================================================
def P(var, default):
    return os.environ.get(var, default)


ROOT = urgg_in("figure_inputs")
G63_F = Path(P("G63_F", str(ROOT / "Fig4/supp/63_genome_thiamine_guilds.tsv")))
GLMM_F = Path(P("GLMM_F", str(ROOT / "Fig4/supp/table_s49_thiamine_guild_presence_models.csv")))
SHARING_F = Path(P("SHARING_F", str(ROOT / "Fig2/supp/representative_host_range.csv")))
GTDB_BAC_F = Path(P("GTDB_BAC_F", str(ROOT / "Fig3/supp/gtdbtk.bac120.summary.tsv")))
GTDB_AR_F = Path(P("GTDB_AR_F", str(ROOT / "Fig3/supp/gtdbtk.ar53.summary.tsv")))
LINEAGE_F = Path(P("LINEAGE_F", str(ROOT / "Fig5/supp/lineage_full_classification_table.csv")))
OUT_DIR = Path(P("OUT_DIR", str(urgg_out("fig4"))))
EXCLUDE = [g.strip() for g in P("EXCLUDE_GENERA", "").split(",") if g.strip()]   # optional manual override
OUT_DIR.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
PT_PER_MM = 72.27 / 25.4            # sizes in mm -> points

TIER_HOST = "Host-specific (1 host species)"
TIER_SHARED = "Shared (2-3 host species)"
TIER_CORE = "Core / broad (>3 host species)"
TIER_LAB = {TIER_CORE: "Core/broad", TIER_SHARED: "Shared", TIER_HOST: "Host-specific"}

# palettes: refined, same hue families across all panels
PRODUCER = dict(dark="#A9423F", mid="#CF6A5F", light="#EDA99E", wash="#FCF1EE", line="#E2A097")
SALVAGER = dict(dark="#1B6B70", mid="#3A9CA0", light="#93CDC9", wash="#EBF6F5", line="#86C4C0")
NEUTRAL = dict(ink="#222222", grey="#7F7D76", light="#BDBAB1", wash="#F4F3EF", pool="#3D4652")
GUILD_COL = {"1_de_novo": PRODUCER["dark"], "2_partial_de_novo": PRODUCER["light"],
             "3_precursor_salvage": SALVAGER["light"], "4_thiamine_salvage": SALVAGER["dark"],
             "5_none_detected": NEUTRAL["light"]}
GUILD_LAB = {"1_de_novo": "Complete de novo", "2_partial_de_novo": "Partial de novo (thiC+)",
             "3_precursor_salvage": "Precursor salvage", "4_thiamine_salvage": "Thiamine salvage",
             "5_none_detected": "No route detected"}
GUILDS = list(GUILD_COL)


def pick_font():
    for fam in ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"]:
        try:
            font_manager.findfont(fam, fallback_to_default=False)
            return fam
        except Exception:
            pass
    return "DejaVu Serif"


FONT = pick_font()
mpl.rcParams.update({"font.family": FONT, "mathtext.fontset": "stix", "pdf.fonttype": 42, "svg.fonttype": "none",
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2, "ytick.major.size": 2})


def host_count(s):
    vals = s.dropna().astype(str)
    if vals.empty:
        return np.nan
    names = set()
    for v in vals:
        names.update(x.strip() for x in re.split(r"[;,|]", v) if x.strip())
    return len(names)


def fmt_int(x):
    return f"{int(round(x)):,}"


# =============================================================================
# PART 2  Data (same selections as the R script)
# =============================================================================
if not G63_F.exists():
    raise SystemExit(f"Missing genome guild table:\n  {G63_F}\nCopy 63_genome_thiamine_guilds.tsv there or set G63_F.")
for f in (GLMM_F, SHARING_F):
    if not f.exists():
        raise SystemExit(f"Missing: {f}")

g = pd.read_csv(G63_F, sep="\t")
sh = pd.read_csv(SHARING_F, usecols=lambda c: c in ("representative", "host_species_list"))
sh = sh.rename(columns={"representative": "genome"})
if "host_species_list" in sh.columns:
    g = g.merge(sh[["genome", "host_species_list"]], on="genome", how="left")
else:
    g["host_species_list"] = np.nan

# genus and domain from GTDB-Tk
tax_parts = []
for f in (GTDB_BAC_F, GTDB_AR_F):
    if f.exists():
        tax_parts.append(pd.read_csv(f, sep="\t", quoting=3, usecols=["user_genome", "classification"],
                                     dtype=str).rename(columns={"user_genome": "genome"}))
if tax_parts:
    tax = pd.concat(tax_parts, ignore_index=True)
    tax["genus"] = tax.classification.str.extract(r"g__([^;]*)", expand=False).fillna("")
    tax["genus"] = tax.genus.where(tax.genus.str.len() > 0, "(no genus)")
    tax["domain"] = np.where(tax.classification.str.contains("d__Archaea", na=False), "Archaea", "Bacteria")
    g = g.merge(tax[["genome", "genus", "domain"]], on="genome", how="left")
g["genus"] = g.get("genus", pd.Series("(no genus)", index=g.index)).fillna("(no genus)")
g["domain"] = g.get("domain", pd.Series("Bacteria", index=g.index)).fillna("Bacteria")

N_ALL = len(g)
print(f"Genomes: {fmt_int(N_ALL)}")

gpct = g.guild.value_counts()


def gp(x):
    return 100 * gpct.get(x, 0) / N_ALL


pct_thiC = 100 * g.hmp_denovo.mean()
tr_by_guild = g.groupby("guild").transport.mean() * 100


def tp(x):
    return tr_by_guild.get(x, np.nan)


# representative genera per guild
g["is_complete"] = g.guild.eq("1_de_novo")
g["is_salv"] = g.guild.eq("4_thiamine_salvage")
g["is_prec"] = g.guild.eq("3_precursor_salvage")
gg = g[g.genus != "(no genus)"]
gen = gg.groupby(["genus", "domain"]).agg(
    n=("genome", "size"), n_thiC=("hmp_denovo", "sum"), n_complete=("is_complete", "sum"),
    n_salv=("is_salv", "sum"), n_prec=("is_prec", "sum"), n_tr=("transport", "sum"),
    n_hosts=("host_species_list", host_count)).reset_index()
for k, c in [("f_thiC", "n_thiC"), ("f_complete", "n_complete"), ("f_salv", "n_salv"), ("f_tr", "n_tr")]:
    gen[k] = gen[c] / gen["n"]
gen["named"] = gen.genus.str.match(r"^[A-Z][a-z]+(_[A-Z])?$")
gen["pct_samples"] = np.nan


def pick_distinct(d, k):
    """skip sister placeholders (Phocaeicola_A after Phocaeicola)"""
    d = d[~d.genus.isin(EXCLUDE)].copy()
    d["base"] = d.genus.str.replace(r"_[A-Z]+$", "", regex=True)
    return d.drop_duplicates("base").head(k)


rep_prod = pick_distinct(gen[(gen.n >= 20) & (gen.f_thiC >= 0.5) & (gen.domain == "Bacteria")]
                         .sort_values(["n_complete", "n_thiC"], ascending=False), 4)
rep_salv = pick_distinct(gen[(gen.n >= 20) & (gen.f_salv >= 0.6) & (gen.f_thiC < 0.05) & (gen.f_tr >= 0.5)]
                         .sort_values("n_salv", ascending=False), 3)
gen.sort_values("n_thiC", ascending=False).to_csv(OUT_DIR / "figA_genus_thiamine_summary.csv", index=False)
print("Representative producer genera:", list(rep_prod.genus))
print("Representative salvager genera:", list(rep_salv.genus))

# RUC lineage
ruc = None
if LINEAGE_F.exists():
    lin = set(pd.read_csv(LINEAGE_F, usecols=["genome"]).genome)
    r = g[g.genome.isin(lin)]
    if len(r):
        ruc = dict(n=len(r), tr=100 * r.transport.mean(), tpk=100 * r.tpk_thiN.mean(),
                   thiL=100 * r.tmp_kinase_thiL.mean(), salv=100 * r.is_salv.mean(),
                   n_hosts=host_count(r.host_species_list))

# GLMM (panel c)
trait_lab = {"hmp_denovo": "thiC (de novo, committed step)", "tpk_thiN": "TPK (activates imported thiamine)",
             "hmp_kinase": "HMP kinase (ThiD/DE/DN)", "denovo_module_ge_thresh": "De novo module (at least 50% complete)",
             "tmp_kinase_thiL": "ThiL (TMP kinase)", "transport": "Thiamine transporter",
             "guild:1_de_novo": "Complete de novo guild", "guild:2_partial_de_novo": "Partial de novo guild (thiC+)",
             "guild:4_thiamine_salvage": "Thiamine-salvage guild", "guild:5_none_detected": "No route detected"}
gm = pd.read_csv(GLMM_F)
full = gm[(gm.trait.isin(trait_lab)) & (gm.adjustment == "adjusted for log genome size + MAG completeness")].copy()
full = full.set_index("trait").loc[list(trait_lab)].reset_index()
full["label"] = full.trait.map(trait_lab)
# fill rule: BH-adjusted LRT p < 0.05 and the within-lineage CI excludes 1
# non-converged fits (status check_convergence) have LRT p = NA: drawn open with a dagger, never filled
full["nonconv"] = full.status.str.contains("check_convergence", na=False)
full["sig"] = (~full.nonconv) & (full.LRT_p_BH < 0.05) & ((full.OR_glmm_low > 1) | (full.OR_glmm_high < 1))
full["plot_label"] = full.label + np.where(full.nonconv, " \u2020", "")
full[["trait", "label", "OR_fisher", "OR_glmm", "OR_glmm_low", "OR_glmm_high", "LRT_p_BH", "status", "nonconv", "sig"]] \
    .rename(columns={"OR_glmm_low": "CI_low", "OR_glmm_high": "CI_high", "sig": "filled"}) \
    .to_csv(OUT_DIR / "panel_c_plot_table.csv", index=False)

# =============================================================================
# PART 3  Figure
# =============================================================================
FIG_W, FIG_H = 183, 170                       # mm
fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))


def ax_at(x0, y0, w, h):
    """axes placed in millimetres from the bottom-left corner of the figure"""
    return fig.add_axes([x0 / FIG_W, y0 / FIG_H, w / FIG_W, h / FIG_H])


# schematic: 158 x 90 mm, so 200 x 114 data units give equal x and y scales (round pool sphere)
ax_s = ax_at(10, 79, 158, 90)
ax_b = ax_at(32, 34, 60, 40)                  # guild composition; room for tier labels on the left
ax_c = ax_at(139, 34, 39, 40)                 # OR panel; room for trait labels on the left

# ---- a: schematic (illustrated style) ----------------------------------------
ax_s.set_xlim(0, 200)
ax_s.set_ylim(-9, 111)
ax_s.set_aspect("equal", adjustable="box")
ax_s.axis("off")

PD, SD = PRODUCER["dark"], SALVAGER["dark"]
SHADOW = "#3A3A3A"


def fs(mm):
    return mm * PT_PER_MM


def rounded(x, y, w, h, r, fc, ec, lw, z, alpha=1.0):
    ax_s.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec,
                                  lw=lw, zorder=z, alpha=alpha))


def node(x, y, label, side, w=15, h=6.4, size=2.5):
    fc = PRODUCER["wash"] if side == "p" else SALVAGER["wash"]
    ec = PRODUCER["line"] if side == "p" else SALVAGER["line"]
    rounded(x - w / 2 + 0.7, y - h / 2 - 0.8, w, h, 1.4, SHADOW, "none", 0, 2, alpha=0.13)    # drop shadow
    rounded(x - w / 2, y - h / 2, w, h, 1.4, fc, ec, 1.0, 3)
    ax_s.text(x, y, label, fontsize=fs(size), color=NEUTRAL["ink"], ha="center", va="center", zorder=4,
              linespacing=0.95)


def arrow(x1, y1, x2, y2, col, lw=1.3):
    ax_s.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>,head_length=0.5,head_width=0.24",
                                   mutation_scale=7.5, color=col, lw=lw, shrinkA=0, shrinkB=0, zorder=2.5))


def label_h(x, y, enzyme, ko, col):
    """enzyme above the arrow, KO below it"""
    ax_s.text(x, y + 2.2, enzyme, fontsize=fs(2.25), color=col, style="italic", weight="bold", ha="center",
              va="center")
    if ko:
        ax_s.text(x, y - 2.3, ko, fontsize=fs(1.85), color=NEUTRAL["grey"], ha="center", va="center")


def label_v(x, y, enzyme, ko, col, ha="left"):
    ax_s.text(x, y, enzyme, fontsize=fs(2.25), color=col, style="italic", weight="bold", ha=ha, va="center",
              linespacing=0.9)
    if ko:
        ax_s.text(x, y - 3.6 - 1.9 * ko.count("\n"), ko, fontsize=fs(1.85), color=NEUTRAL["grey"], ha=ha,
                  va="center", linespacing=0.9)


# guild backgrounds (soft panels)
rounded(2, -7.5, 88, 117, 3.0, PRODUCER["wash"], PRODUCER["line"], 0.7, 0.5, alpha=0.8)
rounded(112, -7.5, 86, 117, 3.0, SALVAGER["wash"], SALVAGER["line"], 0.7, 0.5, alpha=0.8)

# producer panel
ax_s.text(5, 104.6, "De novo producer guild", fontsize=fs(3.3), color=PD, weight="bold", va="center")
ax_s.text(5, 100.4, f"Complete route {gp('1_de_novo'):.1f}%  |  thiC+ (complete or partial) "
                    f"{pct_thiC:.1f}% of {fmt_int(N_ALL)} genomes", fontsize=fs(1.9), color=NEUTRAL["ink"],
          va="center")
node(11, 93, "AIR", "p", w=13)
node(44, 93, "HMP-P", "p")
node(70, 93, "HMP-PP", "p")
node(11, 77, "DXP +\ntyrosine", "p", w=14, h=9.6, size=2.3)
node(44, 77, "HET-P", "p")
node(70, 69, "TMP", "p")
node(70, 54, "TPP", "p")
arrow(17.5, 93, 36.5, 93, PD); label_h(27, 93, "ThiC", "K03147", PD)
arrow(51.5, 93, 62.5, 93, PD); label_h(57, 93, "ThiD", "K00941", PD)
arrow(17.5, 77, 36.5, 77, PD); label_h(27, 77, "ThiH + ThiG", "K03150, K03149", PD)
arrow(70, 89.9, 70, 72.1, PD)
label_v(72.3, 84.6, "ThiE /\nThiDE /\nThiDN", "K00788, K14153,\nK21219", PD)
arrow(51.5, 77, 62.5, 70.6, PD)                                   # HET-P joins TMP (condensation with HMP-PP)
arrow(70, 65.9, 70, 57.1, PD)
label_v(72.3, 64.4, "ThiL", "K00946", PD)                       # raised so K00946 clears the TPP box
ax_s.text(30, 85.2, "de novo synthesis (KEGG M00127, ThiH type)", fontsize=fs(1.95), color=PRODUCER["mid"],
          style="italic", ha="center", va="center")

# salvager panel
ax_s.text(115, 104.6, "Salvager guilds", fontsize=fs(3.3), color=SD, weight="bold", va="center")
ax_s.text(115, 100.4, f"Thiamine salvage {gp('4_thiamine_salvage'):.1f}%  |  precursor salvage "
                      f"{gp('3_precursor_salvage'):.1f}%  |  no route {gp('5_none_detected'):.1f}%",
          fontsize=fs(1.8), color=NEUTRAL["ink"], va="center")
node(122, 93, "Thiamine\n(external)", "s", w=17.5, h=8.4)
node(160, 93, "Thiamine", "s")
node(190, 93, "TPP", "s", w=12)
node(124, 77, "HMP", "s", w=13)
node(124, 63, "HET", "s", w=13)
node(159, 77, "HMP-PP", "s")
node(159, 63, "HET-P", "s")
node(189, 70, "TMP", "s", w=13)

# thiamine transporter (ThiT / ThiBPQ) on the membrane, between external thiamine and the cell
ax_s.add_patch(Ellipse((141.3, 93), 7.2, 4.4, fc=SALVAGER["light"], ec=SALVAGER["dark"], lw=0.9, zorder=3.2))
ax_s.add_patch(Ellipse((141.3, 93), 3.4, 1.9, fc=SALVAGER["dark"], ec="none", alpha=0.45, zorder=3.3))
arrow(130.75, 93, 151.5, 93, SD)
ax_s.text(141.3, 97.6, "ThiT / ThiBPQ", fontsize=fs(2.25), color=SD, style="italic", weight="bold", ha="center",
          va="center")
ax_s.text(141.3, 89.4, "K16789; K02062-4", fontsize=fs(1.85), color=NEUTRAL["grey"], ha="center", va="center")
arrow(167.5, 93, 182, 93, SD); label_h(174.8, 93, "TPK", "K00949", SD)
arrow(130.5, 77, 151.5, 77, SD); label_h(141, 77, "ThiD", "K00941", SD)
arrow(130.5, 63, 151.5, 63, SD); label_h(141, 63, "ThiM", "K00878", SD)
arrow(166.5, 77, 182.5, 71.6, SD)
arrow(166.5, 63, 182.5, 68.4, SD)
ax_s.text(172.5, 70.6, "ThiE", fontsize=fs(2.25), color=SD, style="italic", weight="bold", ha="center", va="center")
ax_s.text(172.5, 67.6, "K00788", fontsize=fs(1.85), color=NEUTRAL["grey"], ha="center", va="center")   # between the two arrows into TMP
arrow(189, 73.1, 189, 89.6, SD)                                  # TMP -> TPP (ThiL)
label_v(187.2, 84.5, "ThiL", "K00946", SD, ha="right")
ax_s.text(141, 84.4, "thiamine salvage", fontsize=fs(1.95), color=SALVAGER["mid"], style="italic",
          ha="center", va="center")
ax_s.text(141, 56.4, "precursor salvage (KEGG M00899)", fontsize=fs(1.95), color=SALVAGER["mid"], style="italic",
          ha="center", va="center")

# thiamine pool with a membrane ring: producer side (left) and salvage side (right) coloured separately
ax_s.add_patch(Circle((102.6, 77.6), 9.2, fc=SHADOW, ec="none", alpha=0.13, zorder=2))
ax_s.add_patch(Circle((102, 79), 9.2, fc=NEUTRAL["pool"], ec="none", zorder=3))
ax_s.add_patch(Arc((102, 79), 21.4, 21.4, theta1=90, theta2=270, color=PRODUCER["mid"], lw=2.4, zorder=3.6))
ax_s.add_patch(Arc((102, 79), 21.4, 21.4, theta1=-90, theta2=90, color=SALVAGER["mid"], lw=2.4, zorder=3.6))
ax_s.add_patch(Circle((102, 79), 9.9, fc="none", ec="white", lw=0.5, alpha=0.6, zorder=3.7))
ax_s.add_patch(Circle((99.6, 81.6), 5.2, fc="white", ec="none", alpha=0.07, zorder=3.1))
ax_s.add_patch(Circle((98.2, 83.2), 2.4, fc="white", ec="none", alpha=0.10, zorder=3.1))
ax_s.text(102, 79, "Thiamine\npool", fontsize=fs(2.4), color="white", weight="bold", ha="center", va="center",
          linespacing=0.85, zorder=4)
ax_s.add_patch(Ellipse((107.0, 88.4), 4.0, 2.4, angle=-28, fc=SALVAGER["light"], ec=SALVAGER["dark"], lw=0.9,
                       zorder=3.8))                                 # transporter on the ring
arrow(77.5, 54, 94.5, 72.5, NEUTRAL["pool"])                     # producers release thiamine
ax_s.text(91.5, 60, "release", fontsize=fs(1.9), color=NEUTRAL["grey"], style="italic", va="center")
arrow(108.0, 86.0, 113.5, 89.6, NEUTRAL["pool"])                 # pool -> external thiamine (via transporter)
arrow(111.2, 77.6, 117, 77, NEUTRAL["pool"])                     # pool -> HMP
arrow(109, 73.2, 117, 65, NEUTRAL["pool"])                       # pool -> HET


# representative genera
def genus_line(r, kind):
    if kind == "p":
        stat = f"complete route {100 * r.f_complete:.0f}%, thiC+ {100 * r.f_thiC:.0f}% of {fmt_int(r.n)} genomes"
    else:
        stat = f"salvage {100 * r.f_salv:.0f}%, transporter {100 * r.f_tr:.0f}% of {fmt_int(r.n)} genomes"
    if pd.notna(r.pct_samples):
        reach = f"{r.pct_samples:.1f}% of samples, {int(r.n_hosts)} host species"
    elif pd.notna(r.n_hosts):
        reach = f"{int(r.n_hosts)} host species"
    else:
        reach = ""
    return "; ".join(x for x in (stat, reach) if x)


def draw_genera(tbl, kind, x0, ytop):
    col = PD if kind == "p" else SD
    y = ytop
    for _, r in tbl.iterrows():
        ax_s.scatter([x0], [y], s=13, marker="o", facecolor=col if r.named else "white", edgecolor=col,
                     linewidth=0.8, zorder=4)
        ax_s.text(x0 + 2.4, y, r.genus, fontsize=fs(2.3), color=NEUTRAL["ink"], style="italic", weight="bold",
                  va="center")
        ax_s.text(x0 + 2.4, y - 3.0, genus_line(r, kind), fontsize=fs(1.85), color=NEUTRAL["grey"], va="center")
        y -= 7.0
    return y


LABEL_MM, PLOT_MM = 14, 34                    # genus-name column and plot width inside each guild panel (mm)
S_X, S_Y = 158 / 200, 90 / 114                # mm per data unit on the schematic axes


def genus_plot(tbl, kind, x_data, ytop, step=7.0):
    """representative genera as a small dot plot: circle = first share, diamond = second share, counts at right"""
    col = PD if kind == "p" else SD
    light = PRODUCER["light"] if kind == "p" else SALVAGER["light"]
    a_key, b_key = ("f_complete", "f_thiC") if kind == "p" else ("f_salv", "f_tr")
    a_name, b_name = ("complete route", "thiC+") if kind == "p" else ("salvage", "transporter")
    n = len(tbl)
    ylo, yhi = ytop - step * (n - 1) - 3.5, ytop + 3.5
    x_mm = 10 + x_data * S_X
    y_mm = 79 + (ylo + 3) * S_Y
    h_mm = (yhi - ylo) * S_Y
    ax = fig.add_axes([(x_mm + LABEL_MM) / FIG_W, y_mm / FIG_H, PLOT_MM / FIG_W, h_mm / FIG_H])
    ax.set_facecolor("#F6DFD9" if kind == "p" else "#D9EEEC")       # tint of each guild's colour, not white
    ys = [ytop - step * i for i in range(n)]
    for y, (_, r) in zip(ys, tbl.iterrows()):
        a, b = 100 * r[a_key], 100 * r[b_key]
        ax.plot([b, a], [y, y], color=NEUTRAL["light"], lw=0.8, zorder=1)
        ax.scatter([a], [y], s=14, facecolor=col if r.named else "white", edgecolor=col, linewidth=0.7, zorder=3)
        ax.scatter([b], [y], marker="D", s=9, facecolor="white", edgecolor=col, linewidth=0.7, zorder=3)
        ax.text(1.04, y, f"n = {fmt_int(r.n)}; {int(r.n_hosts)} hosts", transform=ax.get_yaxis_transform(),
                fontsize=5, va="center", color=NEUTRAL["grey"])
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(["0%", "50%", "100%"], fontsize=5)
    ax.set_yticks(ys)
    ax.set_yticklabels(list(tbl.genus), fontsize=5.8, style="italic", weight="bold", color=NEUTRAL["ink"])
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", length=1.5, pad=1)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x", color="#E6E4DD", linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"% of genomes (circle: {a_name}; diamond: {b_name})", fontsize=5, labelpad=1.5)


ax_s.text(5, 46.5, "Representative producer genera (this catalogue)", fontsize=fs(2.35), color=PD, weight="bold",
          va="center")
if len(rep_prod):
    genus_plot(rep_prod, "p", 7, 41.6)
ax_s.text(115, 46.5, "Representative salvager genera (this catalogue)", fontsize=fs(2.35), color=SD, weight="bold",
          va="center")
nr = min(len(rep_salv), 4 if ruc is None else 3)
salv_tbl = rep_salv.head(nr)
if ruc is not None:                           # RUC lineage shown as a fourth row, same two shares as the genera
    salv_tbl = pd.concat([salv_tbl, pd.DataFrame([dict(genus="RUC lineage", n=ruc["n"], n_hosts=ruc["n_hosts"],
                                                       f_salv=ruc["salv"] / 100, f_tr=ruc["tr"] / 100, named=False)])],
                         ignore_index=True)
if len(salv_tbl):
    genus_plot(salv_tbl, "s", 117, 41.6)

# footnotes

# ---- b: guild composition by tier -------------------------------------------
tiers = [TIER_CORE, TIER_SHARED, TIER_HOST]
ypos = {t: 2 - i for i, t in enumerate(tiers)}
counts = pd.crosstab(g.sharing_category, g.guild).reindex(index=tiers, columns=GUILDS, fill_value=0)
pct = counts.div(counts.sum(axis=1), axis=0) * 100
ntier = counts.sum(axis=1)
left = {t: 0.0 for t in tiers}
for guild in GUILDS:
    for t in tiers:
        w = pct.loc[t, guild]
        if w <= 0:
            continue
        ax_b.barh(ypos[t], w, left=left[t], height=0.66, color=GUILD_COL[guild], edgecolor="white", linewidth=0.4)
        txt_col = "white" if guild in ("1_de_novo", "4_thiamine_salvage") else NEUTRAL["ink"]
        if w >= 6:
            ax_b.text(left[t] + w / 2, ypos[t], f"{w:.0f}", ha="center", va="center", fontsize=6.6, color=txt_col)
        else:
            # segments too narrow for an inside label: small label above the bar
            ax_b.text(left[t] + w / 2, ypos[t] + 0.36, f"{w:.1f}", ha="center", va="bottom", fontsize=5.2,
                      color=NEUTRAL["ink"])
        left[t] += w
ax_b.set_yticks([ypos[t] for t in tiers])
ax_b.set_yticklabels([f"{TIER_LAB[t]}\n(n = {fmt_int(ntier[t])})" for t in tiers], fontsize=6.8)
ax_b.set_ylim(-0.6, 2.7)
ax_b.set_xlim(0, 100)
ax_b.set_xticks([0, 25, 50, 75, 100])
ax_b.set_xticklabels([f"{x}%" for x in [0, 25, 50, 75, 100]], fontsize=6.8)
ax_b.set_xlabel("Genomes in each guild (%)", fontsize=7.4)
for side in ("top", "right", "left"):
    ax_b.spines[side].set_visible(False)
ax_b.tick_params(axis="y", length=0)
ax_b.grid(axis="x", color="#E6E4DD", linewidth=0.4, zorder=0)
ax_b.set_axisbelow(True)
handles_b = [Patch(facecolor=GUILD_COL[k], edgecolor="none", label=GUILD_LAB[k]) for k in GUILDS]
ax_b.legend(handles=handles_b, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2, frameon=False,
            fontsize=6.6, handlelength=1.2, columnspacing=1.2, handletextpad=0.5)

# ---- c: raw vs within-lineage OR ---------------------------------------------
ylab = list(full.plot_label)
yy = {lab: len(ylab) - 1 - i for i, lab in enumerate(ylab)}
for _, r in full.iterrows():
    y = yy[r.plot_label]
    ax_c.errorbar(r.OR_glmm, y, xerr=[[r.OR_glmm - r.OR_glmm_low], [r.OR_glmm_high - r.OR_glmm]], fmt="none",
                  ecolor=NEUTRAL["ink"], elinewidth=0.6, capsize=1.6, capthick=0.6, zorder=2)
    ax_c.scatter([r.OR_fisher], [y], marker="D", s=12, facecolor="none", edgecolor=NEUTRAL["ink"], linewidth=0.6,
                 zorder=5)                                              # on top, so it stays visible over a near-equal circle
    if r.nonconv:
        face = "white"
    else:
        face = PRODUCER["dark"] if r.sig else NEUTRAL["light"]
    ax_c.scatter([r.OR_glmm], [y], marker="o", s=16, linewidth=0.5 if not r.nonconv else 0.8,
                 edgecolor=NEUTRAL["ink"], facecolor=face, zorder=4)
ax_c.axvline(1, color=NEUTRAL["grey"], linewidth=0.5, linestyle=(0, (2, 2)), zorder=1)
ax_c.set_xscale("log")
ax_c.set_xlim(0.2, 7)
ax_c.set_xticks([0.25, 0.5, 1, 2, 4])
ax_c.set_xticklabels(["0.25", "0.5", "1", "2", "4"], fontsize=6.8)
ax_c.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax_c.set_yticks(range(len(ylab)))
ax_c.set_yticklabels(list(reversed(ylab)), fontsize=6.8)
ax_c.set_ylim(-0.6, len(ylab) - 0.4)
ax_c.set_xlabel("Odds ratio, core/broad vs host-specific\n(log scale)", fontsize=7.0)
for side in ("top", "right", "left"):
    ax_c.spines[side].set_visible(False)
ax_c.tick_params(axis="y", length=0)
ax_c.grid(axis="x", color="#E6E4DD", linewidth=0.4, zorder=0)
ax_c.set_axisbelow(True)
handles_c = [Line2D([], [], marker="D", ls="none", markersize=3.6, markerfacecolor="white",
                    markeredgecolor=NEUTRAL["ink"], markeredgewidth=0.6, label="Unadjusted (Fisher)"),
             Line2D([], [], marker="o", ls="none", markersize=3.8, markerfacecolor=NEUTRAL["light"],
                    markeredgecolor=NEUTRAL["ink"], markeredgewidth=0.5, label="Within lineages (GLMM)")]
ax_c.legend(handles=handles_c, loc="upper left", bbox_to_anchor=(-0.05, -0.34), ncol=1, frameon=False,
            fontsize=6.6, handletextpad=0.3, labelspacing=0.5)

# panel tags
for tag_txt, x_mm, y_mm in (("a", 1, 169), ("b", 1, 76), ("c", 97, 76)):
    fig.text(x_mm / FIG_W, y_mm / FIG_H, tag_txt, fontsize=11, weight="bold", va="top", ha="left")

# ---- export ----
base = OUT_DIR / "Fig_thiamine_guilds"
fig.savefig(f"{base}.png", dpi=600, facecolor="white")
fig.savefig(f"{base}.pdf")
fig.savefig(f"{base}.svg")
plt.close(fig)
print("Font:", FONT)
print("Written:")
for f in sorted(OUT_DIR.glob("Fig_thiamine_guilds.*")) + [OUT_DIR / "panel_c_plot_table.csv",
                                                         OUT_DIR / "figA_genus_thiamine_summary.csv"]:
    print(" ", f)
