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
#  SCRIPT  : make_workflow_infographic.py  -  Catalog construction workflow infographic
#  AUTHOR  : QRQ
#  INPUT   : Host icons (ICON_DIR), two-stage dRep method image, RUC lineage image
#  OUTPUT  : Microbial_Genome_Catalog_Workflow_v7.png at OUT (the figure in the README)
# ==========================================================

"""Microbial genome catalog construction workflow (URGG catalog), v6.

Layout: summary on top; host panel on the left (host species + GI sampling sites);
Dataset -> QC in the top row; two-stage dereplication under both; Functional and
Taxonomic side by side; Cooperative-trait genes and the RUC lineage figure side by
side at the bottom.

Numbers follow Table 1 of the manuscript and the virus section draft (v15).
Host silhouettes come from C:\\qrq_metagenome_pipeline\\icons (see test\\make_icons.py).
Layout is in pixel-like units (768 x H, origin top-left).
"""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, Polygon, Wedge
from PIL import Image, ImageOps

OUT = r"C:\qrq_metagenome_pipeline\Plots and Figures\Microbial_Genome_Catalog_Workflow_v7.png"
ICON_DIR = r"C:\qrq_metagenome_pipeline\icons"
DREP_IMG = r"C:\qrq_metagenome_pipeline\manuscript\Figures\Fig1\supp\2stage drep method.png"
RUC_IMG = r"C:\qrq_metagenome_pipeline\images\4.png"

# ---- palette -------------------------------------------------------------
TITLE_BG = "#E4EAF3"
EDGE = "#A7B1C5"
DARK = "#3B4252"
TEXT = "#222B3A"
ARROW = "#5C64A8"
HDR = "#EDF1F6"

W, H = 768, 1170
fig = plt.figure(figsize=(W / 100, H / 100), dpi=150)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)  # y grows downward, like the reference image
ax.axis("off")


# ---- helpers -------------------------------------------------------------
def rbox(x, y, w, h, fc="white", ec=EDGE, lw=1.2, r=10, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=z))


def txt(x, y, s, size=9, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    kw.setdefault("color", TEXT)
    kw.setdefault("zorder", 5)
    ax.text(x, y, s, fontsize=size, **kw)


def text_width(s, size):
    """Width of a string in data units, measured with the real font."""
    t = ax.text(0, 0, s, fontsize=size)
    w = t.get_window_extent(fig.canvas.get_renderer()).width * (W / fig.bbox.width)
    t.remove()
    return w


def justified_paragraph(s, x0, x1, y0, step, size):
    """Wrap s into lines of width x1 - x0; justify every line except the last."""
    words = s.split()
    space = text_width("a a", size) - text_width("aa", size)
    lines, current, used = [], [], 0.0
    for word in words:
        wd = text_width(word, size)
        extra = wd if not current else space + wd
        if current and used + extra > x1 - x0:
            lines.append((current, used))
            current, used = [word], wd
        else:
            current.append(word)
            used += extra
    lines.append((current, used))
    for i, (line_words, _) in enumerate(lines):
        y = y0 + i * step
        if i == len(lines) - 1 or len(line_words) == 1:
            txt(x0, y, " ".join(line_words), size=size, ha="left")
            continue
        widths = [text_width(w_, size) for w_ in line_words]
        gap = (x1 - x0 - sum(widths)) / (len(line_words) - 1)
        x = x0
        for w_, wd in zip(line_words, widths):
            txt(x, y, w_, size=size, ha="left")
            x += wd + gap


def header(x, y, w, h, s, fc, size=11):
    rbox(x, y, w, h, fc=fc, ec="none", r=8, z=2)
    txt(x + w / 2, y + h / 2, s, size=size, weight="bold")


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=2,
                                mutation_scale=14, shrinkA=0, shrinkB=0), zorder=4)


def link(xa, ya, xb, yb, lw=1.2, color=DARK):
    ax.plot([xa, xb], [ya, yb], color=color, lw=lw, zorder=4, solid_capstyle="round")


def helix(cx, cy, h=30, w=8, c1="#7C6FC4", c2="#4C8FC9"):
    t = np.linspace(0, 3 * np.pi, 120)
    y = cy - h / 2 + h * t / (3 * np.pi)
    x1 = cx + w * np.sin(t)
    x2 = cx - w * np.sin(t)
    ax.plot(x1, y, color=c1, lw=2, zorder=4)
    ax.plot(x2, y, color=c2, lw=2, zorder=4)
    for i in range(0, 120, 14):
        ax.plot([x1[i], x2[i]], [y[i], y[i]], color="#9AA3B5", lw=1, zorder=4)


def icon_img(name, cx, cy, box_w=40, box_h=24):
    """Place a silhouette PNG centred on (cx, cy), fitted inside box_w x box_h."""
    im = Image.open(os.path.join(ICON_DIR, name + ".png")).convert("RGBA")
    scale = min(box_w / im.width, box_h / im.height)
    w, h = im.width * scale, im.height * scale
    ax.imshow(np.asarray(im), extent=(cx - w / 2, cx + w / 2, cy + h / 2, cy - h / 2),
              zorder=5, aspect="auto", interpolation="antialiased")


# ---- small icons ---------------------------------------------------------
def icon_db(cx, cy):
    ax.add_patch(Ellipse((cx, cy - 9), 22, 7, fc=DARK, zorder=4))
    ax.add_patch(Rectangle((cx - 11, cy - 9), 22, 18, fc=DARK, zorder=4))
    ax.add_patch(Ellipse((cx, cy + 9), 22, 7, fc=DARK, zorder=4))
    for dy in (-2, 5):
        ax.plot([cx - 11, cx + 11], [cy + dy, cy + dy], color="white", lw=1, zorder=5)


def icon_doc(cx, cy):
    ax.add_patch(Rectangle((cx - 11, cy - 13), 22, 26, fc=DARK, zorder=4))
    for dy in (-6, -1, 4, 9):
        ax.plot([cx - 6, cx + 6], [cy + dy, cy + dy], color="white", lw=1.2, zorder=5)


def icon_people(cx, cy):
    for dx in (-8, 0, 8):
        ax.add_patch(Circle((cx + dx, cy - 6), 3.5, fc=DARK, zorder=4))
        ax.add_patch(Ellipse((cx + dx, cy + 8), 11, 12, fc=DARK, zorder=4))


def icon_funnel(cx, cy):
    ax.add_patch(Polygon([(cx - 13, cy - 12), (cx + 13, cy - 12), (cx + 3, cy + 1),
                          (cx + 3, cy + 12), (cx - 3, cy + 12), (cx - 3, cy + 1)],
                         closed=True, fc="#4A5B9A", zorder=4))


def icon_qc(cx, cy):
    ax.add_patch(Circle((cx, cy), 9, fc="white", ec=DARK, lw=1.8, zorder=4))
    txt(cx, cy, "QC", size=6.5, weight="bold", color=DARK)


# ---- title ---------------------------------------------------------------
rbox(10, 8, 748, 56, fc=TITLE_BG, ec="none", r=8)
txt(384, 24, "Unified Ruminant Gastrointestinal Genome (URGG)", size=15, weight="bold", color=DARK)
txt(384, 47, "Catalog Construction", size=15, weight="bold", color=DARK)

# ---- summary --------------------------------------------------------------
rbox(12, 68, 744, 232)
header(12, 68, 744, 30, "Source Data & Project Summary", HDR, size=12)

cards = [
    (20, "#D9E8F6", "Raw MAG bins scanned:", "207,159", icon_db),
    (264, "#F3D9EC", "Studies (BioProjects):", "39 (54)", icon_doc),
    (508, "#E2E8F2", "Quality-filtered metagenomes:", "1,447", icon_people),
]
for x, fc, label, value, icon in cards:
    rbox(x, 104, 236, 88, fc=fc, ec="none", r=8)
    cx = x + 118
    icon(cx, 124)
    txt(cx, 147, label, size=8.6)
    txt(cx, 165, value, size=10.5, weight="bold")

summary_text = (
    "25 studies were reprocessed from raw reads; public MAGs came from 14 further studies (29 BioProjects). "
    "1,447 quality-filtered metagenomes: 800 were assembled and binned; 1,369 were used for abundance profiling. "
    "Of 207,159 bins, 117,145 passed QC (\u226570% completeness, \u22645% contamination); 38,408 are high-quality. "
    "These yield 25,451 species-level genomes (95% ANI) across 17 host species: 15 ruminant, 2 pseudo-ruminant. "
    "Viral catalog: 15,301 viral sequences (before clustering) grouped into 8,980 vOTUs (non-redundant). "
    "Bacterial hosts were predicted for 6,619 vOTUs (73.7%): BLASTN 3,725, CRISPR-spacer 2,894; "
    "308 (3.4%) carry depolymerase candidates."
)
justified_paragraph(summary_text, x0=24, x1=744, y0=208, step=13, size=7.8)

# ---- left panel: host species and GI sampling sites ------------------------
rbox(12, 330, 236, 825)
header(12, 330, 236, 28, "Host Diversity and GI Sites", HDR, size=10)

HOST_GROUPS = [
    ("Bovidae", "#F6D7D5", "#8E3B36", [
        ("Ovis_aries", "Ovis aries", "Sheep"),
        ("Pantholops_hodgsonii", "Pantholops hodgsonii", "Tibetan antelope"),
        ("Bos_taurus", "Bos taurus", "Cattle"),
        ("Bos_grunniens", "Bos grunniens", "Yak"),
        ("Capra_hircus", "Capra hircus", "Goat"),
        ("Bubalus_bubalis", "Bubalus bubalis", "Buffalo"),
        ("Bos_indicus", "Bos indicus", "Zebu cattle"),
        ("Procapra_picticaudata", "Procapra picticaudata", "Tibetan gazelle"),
        ("Bos_frontalis", "Bos frontalis", "Gaur"),
    ]),
    ("Cervidae", "#DDEFCB", "#3F6B2A", [
        ("Hydropotes_inermis", "Hydropotes inermis", "Water deer"),
        ("Alces_alces", "Alces alces", "Moose"),
        ("Cervus_elaphus", "Cervus elaphus", "Red deer"),
        ("Rangifer_tarandus", "Rangifer tarandus", "Reindeer"),
        ("Alces_americanus", "Alces americanus", "North American moose"),
    ]),
    ("Giraffidae", "#F3E6C8", "#7A5A1E", [
        ("Giraffa_camelopardalis", "Giraffa camelopardalis", "Giraffe"),
    ]),
    ("Camelidae", "#E4DCF2", "#5B4A8A", [
        ("Camelus_bactrianus", "Camelus bactrianus", "Bactrian camel"),
        ("Camelus_dromedarius", "Camelus dromedarius", "Dromedary"),
    ]),
]

ROW = 30
y = 366
for family, bg, fg, species in HOST_GROUPS:
    group_h = 16 + ROW * len(species) + 4
    rbox(14, y, 198, group_h, fc=bg, ec="none", r=8, z=2)
    txt(22, y + 8, f"{family}  \u00b7  {len(species)} species", size=8, weight="bold",
        color=fg, ha="left", zorder=6)
    for i, (icon_name, sci, common) in enumerate(species):
        cy = y + 16 + ROW * i + ROW / 2
        icon_img(icon_name, 40, cy, box_w=40, box_h=24)
        txt(66, cy - 5, sci, size=8, style="italic", ha="left", zorder=6)
        txt(66, cy + 7, common, size=7.4, color="#4B5567", ha="left", zorder=6)
    y += group_h + 6

# GI sampling sites (names as shown in the foregut/hindgut schematic)
GI_SITES = [
    "Rumen contents", "Omasum", "Reticulum", "Abomasum", "Duodenum", "Jejunum",
    "Ileum", "Cecum", "Colon", "Rectum", "Faeces",
]
gi_top = y + 8
txt(22, gi_top, "GI sampling sites", size=8.4, weight="bold", ha="left", color=DARK)
for i, site in enumerate(GI_SITES):
    col, row = (0, i) if i < 6 else (1, i - 6)
    txt(26 + col * 108, gi_top + 22 + row * 17, site, size=8, ha="left", color=TEXT)

# ---- top row: dataset assembly -> QC -------------------------------------
rbox(268, 330, 234, 160, fc="white", ec=EDGE)
txt(385, 350, "Dataset assembly and", size=10, weight="bold")
txt(385, 366, "metadata curation", size=10, weight="bold")
icon_db(310, 404)
icon_doc(385, 404)
icon_people(460, 404)
for x, lines in [(310, ("207,159", "bins")), (385, ("54", "BioProjects")),
                 (460, ("1,447", "metagenomes"))]:
    txt(x, 437, lines[0], size=8, weight="bold")
    txt(x, 451, lines[1], size=7)

rbox(520, 330, 236, 160, fc="white", ec=EDGE)
txt(638, 350, "Quality control, assembly,", size=10, weight="bold")
txt(638, 366, "and binning", size=10, weight="bold")
icon_funnel(578, 404)
icon_qc(600, 402)
helix(690, 404, h=34, w=9)
rbox(538, 436, 200, 18, fc="#F4F7FC", ec=EDGE, lw=0.9, r=5, z=3)
txt(638, 445, "800 samples assembled and binned", size=7.2)

arrow(250, 410, 266, 410)
arrow(504, 410, 518, 410)
arrow(638, 490, 638, 508)

# ---- dereplication: under both boxes ---------------------------------------
rbox(268, 510, 488, 260, fc="white", ec=EDGE)
header(268, 510, 488, 26, "Two-stage genome dereplication", HDR, size=10.5)
drep = Image.open(DREP_IMG).convert("RGB")
drep = drep.crop((0, 185, drep.width, drep.height))  # drop the figure's own title
dw = 472
dh = dw * drep.height / drep.width
ax.imshow(np.asarray(drep), extent=(276, 276 + dw, 540 + dh, 540), zorder=5,
          aspect="auto", interpolation="antialiased")

arrow(386, 770, 386, 783)

# ---- row 1: functional annotation | taxonomic classification ---------------
R1_TOP, R1_H = 785, 132
# Functional annotation
rbox(268, R1_TOP, 236, R1_H, fc="white", ec=EDGE)
header(268, R1_TOP, 236, 22, "Functional annotation", "#CDEAE6", size=10)
cols = ["#4C8FC9", "#E08F5B", "#8FC9C2"]
for cx, cy, labels in [
    (292, 840, ["Cohorts", "Functional", "Categories"]),
    (292, 892, ["Basic", "Functional", "Categories"]),
]:
    for (a1, a2), c in zip([(0, 150), (150, 230), (230, 360)], cols):
        ax.add_patch(Wedge((cx, cy), 16, a1, a2, width=7, fc=c, ec="white", lw=0.6, zorder=4))
    for i, (lab, c) in enumerate(zip(labels, cols)):
        ly = cy - 12 + i * 12
        ax.add_patch(Rectangle((318, ly - 3), 6, 6, fc=c, zorder=4))
        txt(328, ly, lab, size=7.8, ha="left")
link(424, 808, 424, 906, lw=1, color=DARK)
link(424, 906, 470, 906, lw=1, color=DARK)
for i, (h, c) in enumerate(zip([44, 34, 22, 12], cols + ["#C9A7E0"])):
    ax.add_patch(Rectangle((432 + i * 9, 906 - h), 7, h, fc=c, zorder=4))
txt(412, 858, "Genes", size=6.8, rotation=90)

# Taxonomic classification
rbox(520, R1_TOP, 236, R1_H, fc="white", ec=EDGE)
header(520, R1_TOP, 236, 22, "Taxonomic classification", "#D3EBC4", size=10)


def cherry(xc, y1, y2, xp):
    link(xc, y1, xp, y1)
    link(xc, y2, xp, y2)
    link(xp, y1, xp, y2)
    return xp, (y1 + y2) / 2


def merge(a, b, xp):
    (xa, ya), (xb, yb) = a, b
    link(xa, ya, xp, ya)
    link(xb, yb, xp, yb)
    link(xp, ya, xp, yb)
    return xp, (ya + yb) / 2


n1 = cherry(596, 830, 842, 576)
n2 = cherry(596, 858, 870, 576)
n3 = cherry(596, 886, 898, 576)
n4 = merge(n1, n2, 558)
n5 = merge(n4, n3, 540)
link(540, n5[1], 528, n5[1])

for base in (838, 906):
    heights = [24, 18, 12, 7, 4]
    link(610, base - 28, 610, base, lw=1, color=DARK)
    link(610, base, 664, base, lw=1, color=DARK)
    for i, (h, c) in enumerate(zip(heights, ["#4C8FC9", "#E08F5B", "#8FC9C2", "#A7D68E", "#C9A7E0"])):
        ax.add_patch(Rectangle((616 + i * 9, base - h), 6, h, fc=c, zorder=4))
txt(700, 872, "Novel", size=7.4)
txt(700, 884, "phyla", size=7.4)

# ---- row 2: cooperative-trait genes | RUC lineage (figure only) ------------
R2_TOP, R2_H = R1_TOP + R1_H + 12, 226
rbox(268, R2_TOP, 236, R2_H, fc="white", ec=EDGE)
header(268, R2_TOP, 236, 34, "Cooperative-trait gene\nclassification", "#E4CFEB", size=9.5)
cx0, cy0 = 386, R2_TOP + 80
nodes = [(cx0 - 36, cy0 - 6), (cx0 - 16, cy0 - 20), (cx0 + 4, cy0 - 6),
         (cx0 + 22, cy0 - 20), (cx0 - 14, cy0 + 12), (cx0 + 6, cy0 + 18)]
for a, b in [(0, 1), (1, 2), (2, 3), (0, 4), (4, 5), (5, 2)]:
    link(*nodes[a], *nodes[b], lw=1.3, color="#6D7590")
for (nx, ny), c in zip(nodes, ["#7C6FC4", "#4C8FC9", "#7C6FC4", "#E08F5B", "#4C8FC9", "#E08F5B"]):
    ax.add_patch(Circle((nx, ny), 4.5, fc=c, zorder=5))
txt(cx0, cy0 + 36, "Metabolic pathways", size=7.6)
for i, line in enumerate(["46,918,358 predicted genes", "Median 1,784 per genome",
                          "37,129,697 non-redundant genes", "(95% identity, 90% coverage)"]):
    txt(cx0, cy0 + 60 + i * 14, line, size=7.6)

rbox(520, R2_TOP, 236, R2_H, fc="white", ec=EDGE)
header(520, R2_TOP, 236, 26, "RUC lineage", "#F1E2C6", size=10.5)
ruc = Image.open(RUC_IMG).convert("RGB")
bbox = ImageOps.invert(ruc.convert("L")).point(lambda v: 255 if v > 12 else 0).getbbox()
if bbox:
    ruc = ruc.crop(bbox)
max_w, max_h = 220, R2_H - 26 - 10
scale = min(max_w / ruc.width, max_h / ruc.height)
rw, rh = ruc.width * scale, ruc.height * scale
rx = 520 + (236 - rw) / 2
ry = R2_TOP + 26 + (max_h - rh) / 2 + 4
ax.imshow(np.asarray(ruc), extent=(rx, rx + rw, ry + rh, ry), zorder=5,
          aspect="auto", interpolation="antialiased")

ax.set_xlim(0, W)
ax.set_ylim(H, 0)
fig.savefig(OUT, dpi=150, facecolor="white")
print("saved", OUT)

