# URGG: Unified Ruminant Gut Genome catalog — code

Code for *Ruminant gut genome catalog links microbial novelty to narrow host range and function to lineage* (Raza et al., submitted to *Nature Communications*).

The catalog holds 25,451 species-level genomes (URGG00001–URGG25451) and 23,999 strain-level genomes from 117,145 quality-filtered MAGs of 17 ruminant and camelid host species, a 37,129,697-gene non-redundant catalog and 8,980 vOTUs.

- **Data:** NGDC OMIX [OMIX-ACCESSION]; raw reads reprocessed for this study under BioProject [PRJCA-ACCESSION].
- **Numbers behind every figure:** the paper's Source Data file.
- **Archived version of this repository:** Zenodo, doi:[ZENODO-DOI].

## Layout

```
pipeline/   1p–13p   building the catalog: reads → MAGs → species/strains → taxonomy, genes, annotation → virome
analysis/   1a–33a   statistics behind the reported results
figures/    1f–14f   figure panels (_setup/ holds shared helpers)
```

Scripts are numbered in run order. Scripts that combine several steps run all of them by default, or only the steps named on the command line:

```bash
bash pipeline/8p_virus_detection.sh                  # batches, virsorter2, vibrant, checkv
bash pipeline/8p_virus_detection.sh checkv           # one step
Rscript analysis/24a_thiamine.R hosts                # one part
```

All parameters are those stated in the paper's Methods. Every script starts with a header giving its purpose, inputs and outputs (paths relative to `URGG_PROJECT`). Tools are expected on `PATH`; versions are listed in the Methods.

## Paths

No script contains a machine-specific path. Locations are read from environment variables:

| Variable | Location |
|---|---|
| `URGG_PROJECT` | project folder (`drep_analysis/`, `annotation/`, `phase_5_virome/`, …) |
| `URGG_READS_ROOT` | read-processing folder (1p) |
| `URGG_BINS` | MAG/bin folders of each BioProject (2p) |
| `URGG_HOCORT` | clean read pairs used for read mapping (7a) |
| `URGG_DB` | reference databases (GTDB R232, NCBI nr, eggNOG, dbCAN, VFDB, AMRFinderPlus, Platon, MOB-suite, CheckV, VirSorter2, VIBRANT, DRAM-v, PhaBOX2, RVD, RGVC, RGMC, …) |
| `URGG_EXTDATA`, `URGG_UPGG` | external genome collections used in 4a (RGMGC, Hungate1000; pig gut catalog) |
| `CATALOGUE_R_DIR` | plotting helpers sourced by 8a (default `analysis/helpers`) |

```bash
export URGG_PROJECT=/data/urgg URGG_DB=/data/databases
bash pipeline/5p_taxonomy_and_gene_catalog.sh
```

Figure scripts read `URGG_PROJECT` and `URGG_FIGURES` (output folder) from the environment, or from `config/local_paths.R` / `config/local_paths.py` (copy the `.example` files).

## Pipeline

| Script | Steps | Paper item |
|---|---|---|
| `config_reads.sh`, `config_annotation.sh`, `config_virome.sh` | shared paths and settings | |
| `1p_read_processing.sh` | qc (fastp v1.0.1) → hostremoval (HoCoRT) → assembly (MEGAHIT v1.2.9) → binning (BASALT v1.2.2) | reprocessed BioProjects |
| `2p_mag_quality_and_inputs.sh` | checkm2 → names → inputs | 117,145 QC genomes; Fig. 1b, 1d |
| `3p_dereplication.sh` | species_stage1 → species_stage2 (95% ANI) → strains (99% ANI) | 25,451 species; 23,999 strains |
| `4p_genome_ids.py` | URGG identifiers | |
| `5p_taxonomy_and_gene_catalog.sh` | gtdbtk (R232) → prodigal → catalog (MMseqs2) → gene_map | taxonomy; 37,129,697 genes |
| `6p_functional_annotation.sh` | eggnog → nr → vfdb → dbcan → dbcan_ids | KEGG, CAZy, VFDB |
| `7p_amr_and_plasmids.sh` | amr → platon → mobsuite | AMR, plasmids |
| `8p_virus_detection.sh` | batches → virsorter2 → vibrant → checkv | viral sequences |
| `9p_virus_annotation.sh` | dramv → prophage | prophage status |
| `10p_votus.sh` | cluster → phabox2 → membership → tables → tiers | 8,980 vOTUs; Fig. 6a, 6d |
| `11p_virus_phylogeny.sh` | terminase → tree | TerL tree |
| `12p_virus_abundance.sh` | coverm | vOTU abundance |
| `13p_virome_comparison.sh` (+ `13p_classify_matches.py`) | rvd → rgvc | Fig. 6f |

## Analyses

| Script | Result |
|---|---|
| `1a_host_sharing_tiers.R` | genome provenance; host-sharing tiers (20,186 / 4,811 / 454) |
| `2a_taxonomic_novelty.R`, `3a_gtdbtk_novelty_summary.py` | taxonomic novelty |
| `4a_reference_comparison.sh` | matches to ten reference collections; Fig. 1c |
| `5a_genome_to_species.py`, `6a_species_accumulation.py` | species accumulation, Chao1; Suppl. Fig. 5 |
| `7a_read_mapping_abundance.sh` | sample manifest (also used by 12p); genome abundance per sample |
| `8a_diversity.R` (+ `helpers/`) | diversity, PERMANOVA, envfit; Fig. 2b |
| `9a_readbased_host_range.R` | read-based host range; Fig. 2c, 2d |
| `10a_cooperative_genes_and_pathways.R` | cooperative-gene shares and categories, gene-level KO join, pathway completeness |
| `11a_cazy_genome_matrix.R`, `12a_cazy_join_tiers.R` | genome × CAZy family table with tiers |
| `13a_cazyme_tiers.R` | CAZy families by tier; Fig. 3a |
| `14a_subpathways.R` | metabolic sub-pathways; Fig. 3b |
| `15a_trait_presence_glmm.R`, `16a_run_trait_models.sh` | within-lineage models for every trait family |
| `17a_cooperative_glmm.R` | cooperative genes by tier; Fig. 3c |
| `18a_select_tree_genomes.R`, `19a_build_subset_tree.sh` | subset tree for correlated evolution |
| `20a_cooperative_phylo_mcmcglmm.R`, `21a_cooperative_table.R` | phylogenetic sensitivity; cooperative-category table |
| `22a_gain_loss_castor.R`, `23a_pagel_correlated_evolution.R` | gain/loss rates, correlated evolution; Fig. 3d |
| `24a_thiamine.R` | thiamine guilds; Fig. 4 |
| `25a_virulence_amr_plasmids.R` | virulence factors, AMR, plasmids; Fig. 6b |
| `26a_ruc_lineage_table.R`, `27a_ruc_lineage_tree.sh`, `28a_ruc_tree_annotation.R` | RUC lineage genomes and tree |
| `29a_ruc_lineage_traits.R` | RUC lineage functional profile; Fig. 5b, 5c |
| `30a_target_genome_tree.sh`, `31a_target_tree_itol.R`, `32a_target_sequence_identity.R` | placement of the unclassified sheep genome; Fig. 5a |
| `33a_manuscript_numbers.R` | numbers quoted in the text |

## Figures

| Script | Panel |
|---|---|
| `1f`–`3f` | Fig. 1a, 1e, 1f |
| `4f`–`7f` | Fig. 3a, 3b, 3c, 3d |
| `8f`–`11f` | Fig. 4a–c, 4d, 4e, 4f |
| `12f`, `13f`, `14f` | Fig. 6b; 6c and 6h; 6e |

Fig. 1a, 1e, 1f and 6e read the catalog summary tables (genome and vOTU tables deposited in OMIX), placed under `URGG_PROJECT/catalog_tables/`. Fig. 5a was drawn in iTOL (https://itol.embl.de) from the tree of 30a and the datasets of 31a. Final panel arrangement was done outside these scripts.

## Requirements

- **Tools:** as listed in the Methods (CheckM2, dRep, skani, GTDB-Tk, Prodigal, MMseqs2, DIAMOND, eggNOG-mapper, run_dbcan, AMRFinderPlus, Platon, MOB-suite, VirSorter2, VIBRANT, CheckV, DRAM-v, PhaBOX2, MAFFT, trimAl, IQ-TREE, CoverM, fastp, HoCoRT, MEGAHIT, BASALT).
- **R ≥ 4.3:** data.table, dplyr, readr, tidyr, ggplot2, patchwork, vegan, lme4, MCMCglmm, ape, castor, phytools, ggalluvial, ggforce.
- **Python ≥ 3.9:** pandas, numpy, scipy, matplotlib.

## Licence

Code: MIT (see `LICENSE`). Data are released under CC BY 4.0 through OMIX.

## Citation

Raza Q. et al. Ruminant gut genome catalog links microbial novelty to narrow host range and function to lineage. *Nature Communications* (submitted). Code: doi:[ZENODO-DOI].
