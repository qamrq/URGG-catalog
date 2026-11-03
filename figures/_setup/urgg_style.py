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
#  SCRIPT  : urgg_style.py  -  Shared colours, fonts and save() for Python figure scripts
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : style constants and helpers
# ==========================================================

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
