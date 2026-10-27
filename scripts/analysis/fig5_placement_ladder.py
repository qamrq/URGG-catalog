import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

CSV_PATH = Path(r"D:\QRQ_ZJU\MAG-SAG\Manuscript\Figures\Fig5\supp\placement_support_table.csv")
OUT_PATH = Path(r"D:\QRQ_ZJU\MAG-SAG\Manuscript\Figures\Fig5\supp\fig_target_placement_summary.png")

# Group colours follow the existing Fig. 5 palette in the reference layout
COLOR = {
    "Outgroup": "#595959",
    "GTDB SK-Y3": "#E8A02A",
    "RUC + GTDB UMGS1840": "#8E4EA8",
    "GTDB RGIG8785": "#4A9ACB",
}
TARGET_COLOR = "#E31A1C"

# Internal node UFBoot values, read from the reference layout (not in the CSV).
# Tip order top to bottom for each run, and the node values A, B, C, D along the
# ladder: A = outgroup node, B = SK-Y3 node, C = RUC node, D = RGIG8785 + sheep node.
LADDER = {
    "Run A: Bacilli outgroup": {
        "order": ["Outgroup", "GTDB SK-Y3", "RUC + GTDB UMGS1840", "GTDB RGIG8785"],
        "nodes": {"A": 100, "B": 100, "C": 47},
    },
    "Run B: Oscillospirales outgroup": {
        "order": ["Outgroup", "GTDB RGIG8785", "RUC + GTDB UMGS1840", "GTDB SK-Y3"],
        "nodes": {"A": 100, "B": 79, "C": 58},
    },
}

# Read group sizes and support values from the placement table
rows = {}
with CSV_PATH.open(encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        rows.setdefault(r["run"], {})[r["group"]] = (int(r["n_tips"]), r["group_clade_UFBoot"])

# Ladder geometry in image pixel-like units (y grows downward)
NODE_X = {"A": 70, "B": 157, "C": 243, "D": 330}
NODE_Y = {"A": 345, "B": 393, "C": 441, "D": 490}
TIP_X = 414
TIP_Y = [150, 248, 345, 442, 540]  # 5 tips: four groups + unclassified sheep genome
TRI_W, TRI_H = 69, 50
LABEL_X = 496

fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
fig.patch.set_facecolor("white")

for ax, (run, spec) in zip(axes, LADDER.items()):
    order = spec["order"]
    nodes = spec["nodes"]
    groups = rows[run]

    # Spine: A -> B -> C -> D
    for a, b in [("A", "B"), ("B", "C"), ("C", "D")]:
        ax.plot([NODE_X[a], NODE_X[b]], [NODE_Y[a], NODE_Y[b]], color="black", lw=2.4)

    # Tip branches: A to outgroup, B to SK-Y3 slot, C to RUC slot, D to RGIG8785 and sheep
    attach = {"A": 0, "B": 1, "C": 2, "D": 3}
    tip_nodes = [("A", 0), ("B", 1), ("C", 2), ("D", 3)]
    for node, tip_idx in tip_nodes:
        ax.plot([NODE_X[node], TIP_X], [NODE_Y[node], TIP_Y[tip_idx]], color="black", lw=2.4)
    ax.plot([NODE_X["D"], TIP_X], [NODE_Y["D"], TIP_Y[4]], color="black", lw=2.4)

    # Collapsed group triangles and labels
    for idx, g in enumerate(order):
        y = TIP_Y[idx]
        n, sup = groups[g]
        tri = Polygon([[TIP_X, y], [TIP_X + TRI_W, y - TRI_H / 2], [TIP_X + TRI_W, y + TRI_H / 2]],
                      closed=True, facecolor=COLOR[g], edgecolor=COLOR[g], lw=0.8)
        ax.add_patch(tri)
        if sup == "NA":
            text = f"{g}\nn = {n}"
        else:
            text = f"{g}\nn = {n}; UFBoot {int(float(sup))}"
        ax.text(LABEL_X, y, text, va="center", fontsize=10.5, linespacing=1.25)

    # Unclassified sheep genome tip (star + bold red label)
    y_s = TIP_Y[4]
    ax.plot([TIP_X + 4], [y_s], marker=(6, 2, 0), markersize=13, color=TARGET_COLOR,
            markeredgewidth=2.2)
    ax.text(TIP_X + 30, y_s, "Unclassified sheep genome", va="center", fontsize=11,
            color=TARGET_COLOR, fontweight="bold")

    # Node support values: black when >= 95, red when < 95
    for key, val in nodes.items():
        color = "black" if val >= 95 else "#D0312D"
        # Place on the lower-left side of the spine so the text does not sit on a branch
        ax.text(NODE_X[key] - 10, NODE_Y[key] + 20, str(val), fontsize=11, fontweight="bold",
                color=color, ha="right", va="center")

    ax.set_title(run, fontsize=13, fontweight="bold", loc="left", x=0.0)
    ax.set_xlim(40, 1000)
    ax.set_ylim(600, 100)
    ax.axis("off")

fig.text(0.02, 0.02,
         "Red star: unclassified sheep genome (target).\n"
         "Node values: UFBoot support; black = 95 or higher, red = below 95 (unresolved).\n"
         "Triangles: collapsed groups.",
         fontsize=10)
fig.tight_layout(rect=[0, 0.05, 1, 1])
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PATH, dpi=300, facecolor="white")
print("saved", OUT_PATH)
