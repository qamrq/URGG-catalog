# URGG catalog: construction and analysis scripts

The Unified Ruminant Gastrointestinal Genome (URGG) catalog comprises 25,451 species-level genomes and 37.1 million non-redundant genes.
Website: The URGG genome, gene and viral catalogs are available at https://alphaindex.zju.edu.cn/urgg

## Catalog at a glance

<img width="1152" height="1755" alt="Microbial_Genome_Catalog_Workflow_v7" src="https://github.com/user-attachments/assets/7cd9c108-56bf-4612-b95a-d4171932372f" />


## Repository layout

```
scripts/
  analysis/   figure and analysis scripts (thiamine guilds, CAZy, AMR, placement)
  figures/    workflow infographic and host-icon generation
```

## Pipeline

1. **Reprocessing**: raw reads from 25 studies were run through one uniform pipeline.
2. **Public MAGs**: collected from 14 studies (29 BioProjects).
3. **QC, assembly and binning**: 1,447 quality-filtered metagenomes, 800 assembled and binned.
4. **Dereplication**: two stages, species (95% ANI) then strain (99% ANI).
5. **Annotation and analysis**: functional, CAZy, thiamine, AMR and cooperative gene analyses (`scripts/analysis/`).

## Script index

| Script | Produces |
|---|---|
| `analysis/figA_thiamine_guilds.py` | Thiamine producer and salvager guilds (Fig. A)
| `analysis/figC_a_phylum_tier_cazy.py` | Phylum × tier × CAZy family (Fig. 3c, working version)
| `analysis/figC_a_phylum_tier_cazy_pub.py` | Same as above, publication layout with phylum pies
| `analysis/figC_b_subpathway_phylum.py` | Sub-pathway phylum composition and tier prevalence (Fig. 3b)
| `analysis/figure6b_amr_heatmap.py` | AMR gene families across host species (Fig. 6b)
| `analysis/fig5_placement_ladder.py` | Target placement summary, ladder layout (Fig. 5 supplement)
| `analysis/fig5_placement_tree.py` | Target placement summary, tree layout (Fig. 5 supplement)
| `figures/make_workflow_infographic.py` | Workflow infographic (`docs/assets/workflow.png`)
| `figures/make_host_icons.py` | Transparent host icons from a silhouette sheet | silhouette sheet (not included) |
| `tests/make_mock_63.py` | Mock thiamine table for layout tests only | host-range table |

## Requirements

- Python 3.10 or newer and R version v4.5.3
- `numpy`, `pandas`, `matplotlib`, `scipy`, `Pillow`

```powershell
pip install numpy pandas matplotlib scipy pillow
```

## Data availability

Source data can be requested from the Author to run or test the scripts. 

## Notes

- Further information on research design, metagenome analysis, and results can be inquired from QRQ (qamar@zju.edu.cn)

## License

TBD. Add a `LICENSE` file before making the repository public.
