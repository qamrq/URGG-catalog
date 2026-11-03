# URGG: Unified Ruminant Gut Genome catalog — code

Code for *Ruminant gut genome catalog links microbial novelty to narrow host range and function to lineage* (Raza et al., submitted to *Nature Communications*).

The catalog holds 25,451 species-level genomes (URGG00001–URGG25451) and 23,999 strain-level genomes (URGGS00001–URGGS23999) from 117,145 quality-filtered MAGs of 17 ruminant and camelid host species, a 37,129,697-gene non-redundant catalog and 8,980 vOTUs.

- **Data:** NGDC OMIX [OMIX-ACCESSION]; raw reads reprocessed for this study under BioProject [PRJCA-ACCESSION]; browsable at https://alphaindex.zju.edu.cn/urgg.
- **Numbers behind every figure:** the paper's Source Data file.
- **Archived version of this repository:** Zenodo, doi:[ZENODO-DOI].

## Repository layout

```
pipeline/
  01_catalog/      QC, species (95% ANI) and strain (99% ANI) dereplication, URGG IDs, accumulation
  02_annotation/   GTDB-Tk, gene prediction, gene catalog, functional annotation, AMR, plasmids
  03_virome/       viral detection, CheckV, vOTUs, taxonomy, host prediction, lifestyle, comparisons
analysis/          statistics behind the reported results (host range, novelty, functional tests,
                   gain/loss, thiamine, AMR/VFDB, RUC lineage)
  lineage_iqtree/  RUC lineage tree and iTOL annotation
release/           tables and files for the website and the OMIX deposit
figures/           panel scripts (one per published panel or panel group)
config/            path templates for the figure scripts
```

All parameters are those stated in the paper's Methods. Pipeline and analysis scripts were run on a Linux server (SLURM-free, bash) and contain that server's paths, set at the top of each script or in `pipeline/02_annotation/config.sh` and `pipeline/03_virome/config.sh`; edit these before running elsewhere. Tool versions are listed in the Methods and are collected by `pipeline/02_annotation/71_collect_versions.sh`.

## Pipeline: run order

| Order | Script | Paper item |
|---|---|---|
| 1 | `01_catalog/01_batch_checkm2.sh` | CheckM2 quality (≥70% completeness, ≤5% contamination); Fig. 1b, 1d |
| 2 | `01_catalog/02_1_standardize_species_names.sh`, `03_prep_basalt_for_drep.sh` | genome naming and input lists |
| 3 | `01_catalog/05a_drep_95_stage1_per_host.sh` → `05b_drep_95_stage2_final.sh` | two-stage species clustering (95% ANI): 25,451 species |
| 4 | `01_catalog/04b_drep_99_strain.sh` | strain clustering (99% ANI): 23,999 strains |
| 5 | `01_catalog/14a_build_urgg_mapping.py`, `14b_create_release_copies.sh` | URGG / URGGS identifiers |
| 6 | `01_catalog/86a_build_genome_to_species.py`, `86b_species_accumulation.py` | species accumulation, Chao1 (Supplementary) |
| 7 | `02_annotation/00_install_tools.sh`, `00b_download_databases.sh`, `config.sh` | tools and databases |
| 8 | `02_annotation/01_gtdbtk.sh` | GTDB-Tk taxonomy |
| 9 | `02_annotation/02_prodigal.sh` | 46,918,358 predicted genes |
| 10 | `02_annotation/03_mmseqs2.sh` | 37,129,697 non-redundant genes |
| 11 | `02_annotation/04a_diamond_nr.sh`, `04b_diamond_vfdb.sh`, `04c_eggnog.sh`, `04d_dbcan.sh`, `53_dbcan_restore_gene_ids.sh` | functional annotation (KEGG, COG, CAZy, VFDB) |
| 12 | `02_annotation/05_amrfinderplus.sh`, `06_platon.sh`, `06b_mobsuite.sh` | AMR genes, plasmids |
| 13 | `02_annotation/13a_build_gene_to_genome_map.sh`, `13b_build_per_host_tables.py`, `13c_remap_hosts_to_species.py` | gene → genome → host tables |
| 14 | `03_virome/00_install_and_databases_virome.sh`, `config.sh` | viral tools and databases |
| 15 | `03_virome/01_virsorter2.sh`, `02_vibrant.sh`, `03_merge_checkv.sh` | viral sequences and quality |
| 16 | `03_virome/04_dramv.sh`, `05_prophage_flagging.sh` | viral annotation, prophage status |
| 17 | `03_virome/06_votu_clustering.sh`, `07_phabox2.sh`, `07b_explode_votu_membership.sh` | 8,980 vOTUs; taxonomy, lifestyle, host (PhaBOX2/CHERRY) |
| 18 | `03_virome/08a_generate_results_tables.sh`, `08f_host_sharing_tiers_final.sh` | vOTU tables and host-sharing tiers |
| 19 | `03_virome/12a_extract_terminase_marker.sh`, `12b_align_trim_iqtree.sh` | terminase phylogeny |
| 20 | `03_virome/14_coverm_virus_abundance.sh` | vOTU abundance |
| 21 | `03_virome/33_rvd_virome_comparison.sh`, `34_rgvc_virome_comparison.sh` | comparison with RVD and RGVC; Fig. 6f |
| 22 | `03_virome/47_cross_host_sankey_and_jaccard.R` | cross-host vOTUs and Jaccard; Fig. 6c, 6d, 6h |

## Analyses: script for each result

| Result | Script(s) in `analysis/` |
|---|---|
| Host-sharing tiers (20,186 host-specific, 4,811 shared, 454 core) | `33_cluster_host_sharing.R` |
| Taxonomic novelty | `12_taxonomic_novelty_from_gtdbtk.py`, `32a_taxonomic_novelty.R` |
| Matches to ten reference collections (skani); Fig. 1c | `32b_catalogue_comparison.sh` |
| Read mapping, abundance, diversity, PERMANOVA; Fig. 2b | `35b_coverm_abundance.sh`, `35c_diversity_analysis.R`, `43_fig5_diversity_composite_biplot_genus.R` |
| Read-based host range of novel vs described species; Fig. 2c, 2d | `80_readbased_host_range_novelty.R` |
| CAZy families by tier, propagated calls; Fig. 3a, 3d | `54_cazyme_propagated_tier_analysis.R` |
| Trait presence models within lineages | `49_trait_presence_taxonomy_glmm.R`, `62_run_adjustment_models.sh` |
| Metabolic sub-pathways; Fig. 3b, Supplementary Table 6 | `50_subpathway_completeness_propagated.R` |
| Cooperative genes; Fig. 3c | `47_cooperative_taxonomy_glmm.R`, `48_cooperative_mcmcglmm_phylo.R`, `73_coop_table8_primary_or.R` |
| Gain/loss rates and correlated evolution; Fig. 3d | `56_castor_cazy_propagated.R`, `79_castor_coop_categories.R`, `57_pagel_cazy_rerun.R` |
| Thiamine guilds; Fig. 4 | `63_thiamine_guild_typing_propagated.R`, `66_thiamine_by_host.R` |
| Virulence factors (VFDB) | `64_vfdb_propagated_tier_analysis.R` |
| AMR genes and plasmids; Fig. 6b | `65_amr_per_genome_check.R`, `78_para126_plasmid_amr.R` |
| RUC lineage tree and rings; Fig. 5a | `lineage_iqtree/01`–`08`, `81_target_lineage_ml_tree.sh` |
| RUC lineage functional profile; Fig. 5b, 5c | `58_ruc_lineage_propagated.R` |
| Unclassified sheep genome: sequence identity | `83_target_sequence_identity_panels.R` |
| Numbers quoted in the text | `61_pull_manuscript_numbers.R` |

## Release tables (website and OMIX)

`release/build_web_tables.py` builds the website tables (`web_genomes.tsv`, `web_votus.tsv`, …); `release/fix_vfdb_columns.py` splits VFDB results by level; `release/build_release.py` writes the download files (species and strain genome archives, gene catalog with URGGG identifiers, vOTU sequences, matrices) and `MD5SUMS.txt`.

## Figures

Figure scripts contain no machine-specific paths. Set two folders once per machine:

| Setting | Meaning |
|---|---|
| `URGG_PROJECT` | folder with the inputs (pipeline outputs laid out as on the server, the release tables under `website_release/web_data/`, and copies of the tables used on the PC under `figure_inputs/FigN/...`) |
| `URGG_FIGURES` | where figures are written (one subfolder per figure) |

Each is read from an environment variable, else from `config/local_paths.R` / `config/local_paths.py` (copy the `.example` files; git-ignored), else defaults to `data/` and `figures_out/` in the repository. Run with `Rscript figures/<fig>/<script>.R` or `python figures/<fig>/<script>.py`.

| Panel | Script | Inputs (relative to `URGG_PROJECT`) |
|---|---|---|
| 1a | `fig1/Fig1a_host_sunburst.py` | `website_release/web_data/web_genomes.tsv` |
| 1e | `fig1/Fig1e_quality_density.py` | `website_release/web_data/web_genomes.tsv` |
| 1f | `fig1/Fig1f_phylum_by_host.py` | `website_release/web_data/web_genomes.tsv` |
| 2b | `analysis/43_fig5_diversity_composite_biplot_genus.R` | CoverM abundance table (server) |
| 3a | `fig3/Fig3a_cazy_phylum_tier.py` | `figure_inputs/Fig2/supp/representative_host_range.csv`, `figure_inputs/Fig3/supp/54_cazy_family_presence_propagated.tsv`, `figure_inputs/Fig3/supp/gtdbtk.{bac120,ar53}.summary.tsv` |
| 3b | `fig3/Fig3b_subpathway_phylum.py` | `figure_inputs/Fig3/supp/subpathway_presence_propagated_50.tsv` and the files of 3a |
| 3c | `fig3/Fig3c_cooperative_ridge_bars.R` | `annotation/results_phase3/s80_figC_tables/figC_c_genome_props.csv`, `figC_c_cooperative.csv` |
| 3d | `fig3/Fig3d_lossgain_scatter.R` | `annotation/results_phase3/s80_figC_tables/figC_d_castor.csv` |
| 4a–c | `fig4/Fig4abc_thiamine_guilds.py` | `figure_inputs/Fig4/supp/63_genome_thiamine_guilds.tsv`, `table_s49_thiamine_guild_presence_models.csv`; the files of 3a; optional `figure_inputs/Fig5/supp/lineage_full_classification_table.csv` |
| 4d–f | `fig4/Fig4d_guild_by_host.R`, `Fig4e_hostpair_thiC_share.R`, `Fig4f_host_network.R` (shared `Fig4_common_setup.R`) | `annotation/results_phase3/s63_…`, `s66_thiamine_by_host/66_host_guilds_all.csv`, `66_hostpair_shared_producers.csv` |
| 5a | iTOL (https://itol.embl.de), tree and datasets from `analysis/lineage_iqtree/` | |
| 6b | `fig6/Fig6b_amr_heatmap.py` | `figure_inputs/Fig6/result_files/amr_heatmap_matrix_reps.tsv`, `amr_family_meta.tsv` |
| 6c, 6d, 6h | `pipeline/03_virome/47_cross_host_sankey_and_jaccard.R` | virome results tables |
| 6e | `fig6/Fig6e_viral_class_and_rank.R` | `website_release/web_data/web_votus.tsv` |

Final figure layout (panel arrangement, silhouettes in Fig. 1a, the Fig. 1b drawing) was assembled outside these scripts.

## Requirements

- **Pipeline:** the tools and versions listed in the Methods (CheckM2, dRep, skani, GTDB-Tk, Prodigal, MMseqs2, DIAMOND, eggNOG-mapper, run_dbcan, AMRFinderPlus, Platon, MOB-suite, VirSorter2, VIBRANT, CheckV, DRAM-v, PhaBOX2, IQ-TREE, CoverM).
- **R ≥ 4.3:** data.table, ggplot2, patchwork, ggforce, vegan, lme4, MCMCglmm, castor, phytools; optional svglite, ragg, systemfonts.
- **Python ≥ 3.9:** pandas, numpy, scipy, matplotlib.

## Licence

Code: MIT (see `LICENSE`). Data are released under CC BY 4.0 through OMIX.

## Citation

Raza Q. et al. Ruminant gut genome catalog links microbial novelty to narrow host range and function to lineage. *Nature Communications* (submitted). Code: doi:[ZENODO-DOI].
