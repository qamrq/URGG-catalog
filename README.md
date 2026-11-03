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

| Script | Process |
|---|---|
| `config_*.sh` | Shared settings |
| `1p` | Read processing, assembly and binning |
| `2p` | MAG quality control |
| `3p` | Species and strain dereplication |
| `4p` | Genome identifiers |
| `5p` | Taxonomy and gene catalog |
| `6p` | Functional annotation |
| `7p` | Antimicrobial resistance and plasmids |
| `8p` | Viral sequence detection |
| `9p` | Viral annotation and prophages |
| `10p` | vOTU clustering, taxonomy and host prediction |
| `11p` | Viral phylogeny |
| `12p` | Viral abundance |
| `13p` | Comparison with viral catalogs |

## Analyses

| Script | Analysis |
|---|---|
| `1a` | Host-sharing tiers |
| `2a`–`3a` | Taxonomic novelty |
| `4a` | Comparison with reference catalogs |
| `5a`–`6a` | Species accumulation |
| `7a` | Read mapping and abundance |
| `8a` | Microbial diversity |
| `9a` | Read-based host range |
| `10a` | Cooperative genes and metabolic pathways |
| `11a`–`13a` | Carbohydrate-active enzymes |
| `14a` | Metabolic sub-pathways |
| `15a`–`16a` | Within-lineage trait models |
| `17a`, `20a`–`21a` | Cooperative traits |
| `18a`–`19a` | Phylogenetic tree |
| `22a`–`23a` | Trait gain, loss and correlated evolution |
| `24a` | Thiamine guilds |
| `25a` | Virulence factors, resistance genes and plasmids |
| `26a`–`29a` | RUC lineage |
| `30a`–`32a` | Placement of the unclassified sheep genome |
| `33a` | Summary statistics |

## Figures

| Script | Figure |
|---|---|
| `1f`–`3f` | Fig. 1 |
| `4f`–`7f` | Fig. 3 |
| `8f`–`11f` | Fig. 4 |
| `12f`–`14f` | Fig. 6 |

Fig. 1a, 1e, 1f and 6e read the catalog summary tables (genome and vOTU tables deposited in OMIX), placed under `URGG_PROJECT/catalog_tables/`. Fig. 5a was drawn in iTOL (https://itol.embl.de). Final panel arrangement was done outside these scripts.

## Requirements

- **Tools:** as listed in the Methods (CheckM2, dRep, skani, GTDB-Tk, Prodigal, MMseqs2, DIAMOND, eggNOG-mapper, run_dbcan, AMRFinderPlus, Platon, MOB-suite, VirSorter2, VIBRANT, CheckV, DRAM-v, PhaBOX2, MAFFT, trimAl, IQ-TREE, CoverM, fastp, HoCoRT, MEGAHIT, BASALT).
- **R ≥ 4.3:** data.table, dplyr, readr, tidyr, ggplot2, patchwork, vegan, lme4, MCMCglmm, ape, castor, phytools, ggalluvial, ggforce.
- **Python ≥ 3.9:** pandas, numpy, scipy, matplotlib.

## Licence

Code: MIT (see `LICENSE`). Data are released under CC BY 4.0 through OMIX.

## Citation

Raza Q. et al. Ruminant gut genome catalog links microbial novelty to narrow host range and function to lineage. *Nature Communications* (submitted). Code: doi:[ZENODO-DOI].
