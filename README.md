# URGG catalog: construction and analysis scripts

Scripts behind the ruminant and pseudo-ruminant genome catalog (URGG catalog).
Website: The URGG genome, gene and viral catalogs are available at https://alphaindex.zju.edu.cn/urgg

## Catalog at a glance



## Repository layout

```
scripts/
  analysis/   figure and analysis scripts (thiamine guilds, CAZy, AMR, placement)
  figures/    workflow infographic and host-icon generation
  tests/      layout-test stand-in data (not real data)
docs/         GitHub Pages site and workflow image
```

## Pipeline

1. **Reprocessing**: raw reads from 25 studies were run through one uniform pipeline.
2. **Public MAGs**: collected from 14 studies (29 BioProjects).
3. **QC, assembly and binning**: 1,447 quality-filtered metagenomes, 800 assembled and binned.
4. **Dereplication**: two stages, species (95% ANI) then strain (99% ANI).
5. **Annotation and analysis**: functional, CAZy, thiamine, AMR and placement analyses (`scripts/analysis/`).

> The construction-stage scripts (steps 1–4) are not yet in this repository. They will be added here.

## Script index

| Script | Produces | Inputs |
|---|---|---|
| `analysis/figA_thiamine_guilds.py` | Thiamine producer and salvager guilds (Fig. A) | `63_genome_thiamine_guilds.tsv`, GLMM table |
| `analysis/figC_a_phylum_tier_cazy.py` | Phylum × tier × CAZy family (Fig. 3c, working version) | tier table, CAZy presence, GTDB-Tk summaries |
| `analysis/figC_a_phylum_tier_cazy_pub.py` | Same as above, publication layout with phylum pies | as above |
| `analysis/figC_b_subpathway_phylum.py` | Sub-pathway phylum composition and tier prevalence (Fig. 3b) | sub-pathway presence, host-range table |
| `analysis/figure6b_amr_heatmap.py` | AMR gene families across host species (Fig. 6b) | AMRFinderPlus output, host-range table |
| `analysis/fig5_placement_ladder.py` | Target placement summary, ladder layout (Fig. 5 supplement) | placement support table |
| `analysis/fig5_placement_tree.py` | Target placement summary, tree layout (Fig. 5 supplement) | placement support table |
| `figures/make_workflow_infographic.py` | Workflow infographic (`docs/assets/workflow.png`) | none |
| `figures/make_host_icons.py` | Transparent host icons from a silhouette sheet | silhouette sheet (not included) |
| `tests/make_mock_63.py` | Mock thiamine table for layout tests only | host-range table |

## Requirements

- Python 3.10 or newer
- `numpy`, `pandas`, `matplotlib`, `scipy`, `Pillow`

```powershell
pip install numpy pandas matplotlib scipy pillow
```

## Running

Scripts read their inputs from paths set near the top of each file. Most of these
currently point to the author's local project folder, so edit those paths before running.
Data files are not in this repository: see *Data availability*.

```powershell
python scripts/analysis/figure6b_amr_heatmap.py <out_dir> 600
```

## Data availability

Input tables (genome-level tables, CAZy and AMR calls, GTDB-Tk summaries) are not stored here.
Add accession numbers or a download link once the data are deposited.

## Notes

- `tests/make_mock_63.py` generates random data for layout checks. Never use its output for results.
- Host silhouettes in the infographic come from an external image sheet. Their license has not been checked, so the icon files are not included in this repository.
- `fig5_placement_ladder.py` and `fig5_placement_tree.py` write to the same output path. Run one at a time, or change `OUT_PATH` in one of them.

## License

TBD. Add a `LICENSE` file before making the repository public.
