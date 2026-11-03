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
#  SCRIPT  : 3f_fig1f_phylum_by_host.py  -  Fig. 1f, phylum composition per host species
#  AUTHOR  : QRQ
#  INPUT   : catalog_tables/web_genomes.tsv
#  OUTPUT  : fig1/Fig1f_phylum_by_host.{png,pdf,svg}, Fig1f_phylum_by_host.csv
# ==========================================================

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "_setup"))
from urgg_paths import urgg_in, urgg_out  # noqa: E402
from urgg_style import PALETTE8, save, set_font  # noqa: E402

g = pd.read_csv(urgg_in("catalog_tables/web_genomes.tsv"), sep="\t",
                usecols=["urgg_id", "rep_host_species", "phylum"])
g["phylum"] = g.phylum.fillna("NA")
top = g.phylum.value_counts().index[:7].tolist()
cats = top + ["Other"]
g["ph"] = np.where(g.phylum.isin(top), g.phylum, "Other")
vc = g.rep_host_species.value_counts()
hosts = sorted(vc.index, key=lambda h: (-vc[h], h))  # ties broken by name
tab = pd.crosstab(g.rep_host_species, g.ph).reindex(index=hosts, columns=cats, fill_value=0)
pct = 100 * tab.div(tab.sum(axis=1), axis=0)
out = urgg_out("fig1")
tab.assign(total=tab.sum(axis=1)).join(pct.add_suffix(" (%)").round(2)).to_csv(out / "Fig1f_phylum_by_host.csv")

set_font(size=8)
fig, ax = plt.subplots(figsize=(7.1, 5.2))
fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.36)
x = np.arange(len(hosts)); bottom = np.zeros(len(hosts))
for c, col in zip(cats, PALETTE8):
    ax.bar(x, pct[c].values, bottom=bottom, color=col, edgecolor="white", linewidth=0.4, width=0.7, label=c)
    bottom += pct[c].values
ax.set_xticks(x)
ax.set_xticklabels(hosts, rotation=45, ha="right", rotation_mode="anchor", style="italic", fontsize=8.5)
ax.set_ylim(0, 100); ax.set_ylabel("Representatives by phylum (%)", fontsize=10.5)
ax.set_xlabel("Host species", fontsize=10.5, labelpad=2)
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="x", length=0)
fig.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.54, 0.0), frameon=False, title="Phylum")
save(fig, out, "Fig1f_phylum_by_host", tight=False)
