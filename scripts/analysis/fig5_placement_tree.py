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
#  SCRIPT  : fig5_placement_tree.py  -  Target placement summary, tree layout (Fig. 5 supplement)
#  AUTHOR  : QRQ
#  INPUT   : placement_support_table.csv
#  OUTPUT  : fig_target_placement_summary.png at OUT_PATH (run one fig5 script at a time)
# ==========================================================

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

CSV_PATH = Path(r"C:\qrq_metagenome_pipeline\manuscript\Figures\Fig5\supp\placement_support_table.csv")
OUT_PATH = Path(r"C:\qrq_metagenome_pipeline\manuscript\Figures\Fig5\supp\fig_target_placement_summary.png")

GROUP_ORDER = ["RUC + GTDB UMGS1840", "GTDB SK-Y3", "GTDB RGIG8785", "Outgroup"]
GROUP_COLOR = {
    "RUC + GTDB UMGS1840": "#1F78B4",
    "GTDB SK-Y3": "#E31A1C",
    "GTDB RGIG8785": "#33A02C",
    "Outgroup": "#999999",
}

# Fallback copy of placement_support_table.csv as read on 2026-10-15 (file was missing at re-run)
FALLBACK_CSV = """run,group,n_tips,group_clade_UFBoot
Run A: Bacilli outgroup,RUC + GTDB UMGS1840,101,100
Run A: Bacilli outgroup,GTDB SK-Y3,11,82
Run A: Bacilli outgroup,GTDB RGIG8785,2,100
Run A: Bacilli outgroup,Outgroup,3,NA
Run B: Oscillospirales outgroup,RUC + GTDB UMGS1840,101,100
Run B: Oscillospirales outgroup,GTDB SK-Y3,11,100
Run B: Oscillospirales outgroup,GTDB RGIG8785,2,100
Run B: Oscillospirales outgroup,Outgroup,4,100"""

if CSV_PATH.exists():
    source = CSV_PATH.open(encoding="utf-8")
    print("reading", CSV_PATH)
else:
    import io
    source = io.StringIO(FALLBACK_CSV)
    print("CSV not found; using fallback copy")

rows = {}
with source:
    for r in csv.DictReader(source):
        rows.setdefault(r["run"], {})[r["group"]] = (int(r["n_tips"]), r["group_clade_UFBoot"])

runs = list(rows.keys())
fig, axes = plt.subplots(1, len(runs), figsize=(11, 5.2), sharey=True)
plt.rcParams["font.family"] = "DejaVu Sans"

for ax, run in zip(axes, runs):
    groups = rows[run]
    ys = {g: len(GROUP_ORDER) - 1 - i for i, g in enumerate(GROUP_ORDER) if g in groups}
    root_x, node_x, tip_x = 0.0, 1.0, 2.2
    y_mid = sum(ys[g] for g in ys if g != "Outgroup") / max(1, len([g for g in ys if g != "Outgroup"]))

    # Root and ingroup/outgroup split
    ax.plot([root_x, node_x], [y_mid, y_mid], color="black", lw=1.4)
    ing = [ys[g] for g in ys if g != "Outgroup"]
    ax.plot([node_x, node_x], [min(ing), max(ing)], color="black", lw=1.4)

    for g, y in ys.items():
        n, sup = groups[g]
        color = GROUP_COLOR[g]
        if g == "Outgroup":
            ax.plot([root_x, root_x], [y_mid, y], color="black", lw=1.4)
            ax.plot([root_x, tip_x], [y, y], color="black", lw=1.4)
        else:
            ax.plot([node_x, tip_x], [y, y], color="black", lw=1.4)

        # Collapsed-clade triangle, width scaled to tip count (log-compressed)
        height = 0.32
        width = 0.35 + 0.22 * (n ** 0.5) / (101 ** 0.5) * 3
        tri = Polygon([[tip_x, y - height / 2], [tip_x, y + height / 2], [tip_x + width, y]],
                      closed=True, facecolor=color, edgecolor="black", lw=0.6, alpha=0.9)
        ax.add_patch(tri)
        # Support value goes into the group label, clear of all branch lines
        if sup == "NA":
            sup_text, lcolor = "UFBoot NA", "#666666"
        else:
            v = int(float(sup))
            sup_text = f"UFBoot {v}"
            lcolor = "black" if v >= 95 else "#D95F02"
        ax.text(tip_x + width + 0.06, y, f"{g}\n(n = {n} tips)", va="center", fontsize=9)
        ax.text(tip_x + width + 0.06, y - 0.27, sup_text, va="center", fontsize=8,
                color=lcolor, fontweight="bold" if lcolor != "#666666" else "normal")

    ax.set_title(run, fontsize=11, loc="left")
    ax.set_xlim(-0.2, 5.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

axes[0].set_ylim(-0.6, len(GROUP_ORDER) - 0.4)
fig.suptitle("Placement of the target genome and reference groups across two outgroup choices",
             fontsize=12, x=0.02, ha="left")
fig.text(0.02, 0.02,
         "Collapsed clades; triangle width scales with tip count. UFBoot support under each group; "
         "values below 95 in orange. Relationships among groups are not shown in the source table.",
         fontsize=8, color="#444444")
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PATH, dpi=300)
print("saved", OUT_PATH)
