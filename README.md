# URGG catalog: construction and analysis scripts

The Unified Ruminant Gastrointestinal Genome (URGG) catalog comprises 25,451 species-level genomes and 37.1 million non-redundant genes.
Website: The URGG genome, gene and viral catalogs are available at https://alphaindex.zju.edu.cn/urgg

## Catalog at a glance

<img width="1152" height="1755" alt="Microbial_Genome_Catalog_Workflow_v7" src="https://github.com/user-attachments/assets/7cd9c108-56bf-4612-b95a-d4171932372f" />


## Repository layout

```
scripts/
  pipeline/   construction steps 01-15 (QC to gene annotation), R; run_pipeline.R runs them in order
  stats/      statistical tests and mixed models (GLMM, MCMCglmm, Fisher, PERMANOVA), R
  analysis/   figure and analysis scripts (thiamine guilds, CAZy, AMR, placement), Python
  figures/    workflow infographic and host-icon generation (Python); R figure scripts in figures/r/
  tests/      layout-test stand-in data (not real data), R
```

## Pipeline

1. **Reprocessing**: raw reads from 25 studies were run through one uniform pipeline.
2. **Public MAGs**: collected from 14 studies (29 BioProjects).
3. **QC, assembly and binning**: 1,447 quality-filtered metagenomes, 800 assembled and binned.
4. **Dereplication**: two stages, species (95% ANI) then strain (99% ANI).
5. **Annotation and analysis**: functional, CAZy, thiamine, AMR and cooperative gene analyses (`scripts/analysis/`).

The construction steps are in `scripts/pipeline/` and run from the repository root:

```r
Rscript scripts/pipeline/run_pipeline.R        # all steps, 01-15
Rscript scripts/pipeline/run_pipeline.R 06 09  # a range of steps
```

Paths and thread counts are set once in `scripts/pipeline/config.R`. External programs (fastp, HoCoRT, MEGAHIT, BASALT, CheckM2, dRep, GTDB-Tk, Prodigal, MMseqs2, eggNOG-mapper, run_dbCAN, DIAMOND) are called from R and must be installed or available in the conda environments named in `config.R`.

## Script index

| Script | Produces |
|---|---|
| `pipeline/01_fastp_qc.R` | Quality-filtered reads (fastp) |
| `pipeline/02_hocort_host_removal.R` | Host and diet read removal (HoCoRT) |
| `pipeline/03_megahit_assembly.R` | Per-sample assemblies (MEGAHIT) |
| `pipeline/04_basalt_binning.R` | Bins and best bin set (BASALT) |
| `pipeline/05_checkm2_quality.R` | Completeness and contamination (CheckM2) |
| `pipeline/06_prepare_drep_inputs.R` | Quality-gated genome pool and dRep inputs |
| `pipeline/07_drep_stage1_per_host.R` | Per-host 95% ANI dereplication (dRep) |
| `pipeline/08_drep_stage2_species.R` | Species catalog at 95% ANI (dRep) |
| `pipeline/09_drep_strain_99.R` | Strain set at 99% ANI (dRep) |
| `pipeline/10_gtdbtk_classify.R` | Taxonomy (GTDB-Tk) |
| `pipeline/11_gene_prediction_prodigal.R` | Genes per genome (Prodigal) |
| `pipeline/12_gene_catalog_mmseqs2.R` | Non-redundant gene catalog (MMseqs2) |
| `pipeline/13_eggnog_annotation.R` | Functional annotation (eggNOG-mapper) |
| `pipeline/14_dbcan_cazyme.R` | CAZyme annotation (run_dbCAN) |
| `pipeline/15_diamond_nr_unassigned.R` | NR search for genes without an eggNOG assignment (DIAMOND) |
| `stats/35a_lineage_cazyme_fisher.R` | Fisher's exact tests on CAZy families, per lineage |
| `stats/35c_diversity_analysis.R` | Shannon diversity, Bray–Curtis PCoA and PERMANOVA |
| `stats/47_cooperative_taxonomy_glmm.R` | Cooperative-gene GLMM with nested taxonomic random effects |
| `stats/48_cooperative_mcmcglmm_phylo.R` | Phylogenetic MCMCglmm for cooperative genes |
| `stats/49_trait_presence_taxonomy_glmm.R` | Trait-presence GLMM for any trait table (run after 47) |
| `analysis/figA_thiamine_guilds.py` | Thiamine producer and salvager guilds (Fig. A) |
| `analysis/figC_a_phylum_tier_cazy.py` | Phylum × tier × CAZy family (Fig. 3c, working version) |
| `analysis/figC_a_phylum_tier_cazy_pub.py` | Same as above, publication layout with phylum pies |
| `analysis/figC_b_subpathway_phylum.py` | Sub-pathway phylum composition and tier prevalence (Fig. 3b) |
| `analysis/figure6b_amr_heatmap.py` | AMR gene families across host species (Fig. 6b) |
| `analysis/fig5_placement_ladder.py` | Target placement summary, ladder layout (Fig. 5 supplement) |
| `analysis/fig5_placement_tree.py` | Target placement summary, tree layout (Fig. 5 supplement) |
| `figures/make_workflow_infographic.py` | Workflow infographic (the figure above) |
| `figures/make_host_icons.py` | Transparent host icons from a silhouette sheet (sheet not included) |
| `figures/r/31g_sunburst_catalogue.R` | Family to species sunburst (Fig. 1a) |
| `figures/r/31j_figure1_composite.R` | Figure 1 five-panel composite |
| `figures/r/36_fig5_diversity_composite.R` | Fig. 5 diversity composite and PCoA |
| `figures/r/37_cooperative_ko__figures.R` | Cooperative-trait figure (Fig. 7) |
| `figures/r/44_cazyme_consolidated_figure.R` | Consolidated CAZyme figure |
| `tests/make_mock_63.R` | Mock thiamine table for layout tests only (not real data) |

The other R figure scripts in `figures/r/` (`00`–`05` helpers and the remaining `31*`, `42`, `43` scripts) are listed in their own file headers.

## Requirements

- R 4.5 or newer (the analyses were run in R 4.5.3), with the packages used by the stats and figure scripts (installed by `figures/r/00_install_packages.R`)
- Python 3.10 or newer for the scripts that are still in Python (`analysis/`, `figures/*.py`):

```powershell
pip install numpy pandas matplotlib scipy pillow
```

- External programs for the pipeline (see the Pipeline section)

## Data availability

Source data can be requested from the Author to run or test the scripts. 

## Notes

- Scripts read their inputs from paths set near the top of each file (or in `config.R` for the pipeline). Most point to the author's local project folder, so edit those paths before running.
- `tests/make_mock_63.R` generates random data for layout checks. Never use its output for results.
- Host silhouettes in the infographic come from an external image sheet. Their license has not been checked, so the icon files are not included in this repository.
- `fig5_placement_ladder.py` and `fig5_placement_tree.py` write to the same output path. Run one at a time, or change `OUT_PATH` in one of them.

- Further information on research design, metagenome analysis, and results can be inquired from QRQ (qamar@zju.edu.cn)

## License

TBD. Add a `LICENSE` file before making the repository public.
