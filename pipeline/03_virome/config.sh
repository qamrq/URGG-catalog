#!/usr/bin/env bash
###############################################################################
#
#                           #                #
#                          ###              ###
#                         #####            #####
#                        #######          #######
#                       #########        #########
#                      ###########      ###########
#                     #############    #############
#                   ###############    ###############
#                  #################  #################
#                 ######################################
#                ########        #######         ########
#               ###########      #######       ###########
#              ##############    #######     ##############
#             ################## #######  ##################
#            ################################################
#            ################################################
#               ##########################################
#                   ##################################
#                       ##########################
#                          ####################
#                           ##################
#                            ################
#                             ##############
#                              ############
#                              ############
#                               ##########
#                                ########
#                                 ######
#
#        P H A S E   5   ::   V I R O M E   A N A L Y S I S
#             rumen MAGs  |  cross-host comparison
#
###############################################################################
# config.sh — Phase 5: Virome Analysis Pipeline (rumen MAGs, cross-host)
#
# Mirrors the structure of Phase 3's config.sh. Edit the paths in this file
# once; every 0X_*.sh script sources this file.
#
# Pipeline design recap (per Yan & Yu 2024, Microbiome, adapted for
# cross-host comparison):
#   1. VirSorter2 (replaces PhageBoost) on dRep 99% ANI MAGs -> viral calls
#   2. CheckV -> quality filter (keep VirSorter2 category 1 & 2)
#   3. VIBRANT (-virome) -> second independent confirmation
#   4. DRAM-v -> functional annotation of confirmed viral sequences
#   5. Prophage flagging (CheckV integrated calls + VOGDB integrase/excisionase/
#      repressor gene hits)
#   6. Cluster confirmed viral seqs into vOTUs at 95% ANI / 85% shorter-contig
#      coverage (CheckV scripts) -- UNIFIED catalog across all host species
#   7. PhaGCN2.1 -> taxonomic classification of vOTUs
#   8. GTDB-Tk host taxonomy (already running/complete from Phase 3/4 -- reused)
#   9. PropagAtE -> active vs cryptic prophage (needs QC'd reads, per-sample)
#  10. MinCED (CRISPR arrays on MAGs) + BLASTn spacer matching -> host linkage
#  11. Prodigal + CARD -> ARGs carried by prophages
###############################################################################
set -euo pipefail

### ---- TOP-LEVEL DIRECTORIES (EDIT THESE) ---------------------------------
BASE_DIR="/disk252/qamar/metagenome_analysis"
PHASE5_DIR="${BASE_DIR}/phase_5_virome"

# Extension of MAG genome files (matches Phase 3 config.sh convention)
export GENOME_EXT=fa

# Input MAG sets (from dRep) -- EDIT to actual paths
# *** VIRAL CALLING (VirSorter2/CheckV/VIBRANT/DRAM-v) USES THE 99% ANI STRAIN
# CATALOG — prophages are strain-specific, and the 95% species reps would
# collapse exactly the resolution this step needs. This is the OPPOSITE of
# Phase 3 annotation, which correctly uses the 95% catalog. Do not swap. ***
#
# The ":-" defaulting lets you override from the command line for test runs
# without editing this file, e.g.:
#   MAGS_99ANI_DIR=/tmp/mag_test_batch bash 01_virsorter2.sh
MAGS_99ANI_DIR="${MAGS_99ANI_DIR:-/disk252/qamar/metagenome_analysis/drep_analysis/drep_99_strain/dereplicated_genomes}"

# 95% ANI species catalog -- used ONLY for clustering confirmed viral seqs
# into vOTUs (step 6), matching the paper's threshold for that specific step.
MAGS_95ANI_DIR="${MAGS_95ANI_DIR:-/disk252/qamar/metagenome_analysis/drep_analysis/drep_95_species_final/dereplicated_genomes}"

# MAG -> host-species / BioProject metadata table (tab-separated).
# Reusing your existing Phase 3 metadata.tsv directly (genome|bioproject|
# host|species|sample_site|...) -- no new table needs to be built.
MAG_METADATA_TSV="${BASE_DIR}/drep_analysis/metadata.tsv"

# QC'd, host-depleted paired-end reads manifest (reused from Phase 4 / SV
# analysis -- same samples_manifest.tsv format: sample_id  fq1  fq2  bioproject  host)
SAMPLES_MANIFEST="${BASE_DIR}/phase4_sv/samples_manifest.tsv"

# GTDB-Tk output (already run in Phase 3/4 -- reused here for host taxonomy)
GTDBTK_OUTPUT_DIR="${BASE_DIR}/annotation/gtdbtk_output_v232"

### ---- PHASE 5 OUTPUT SUBDIRECTORIES ---------------------------------------
VS2_OUT="${PHASE5_DIR}/01_virsorter2"
CHECKV_OUT="${PHASE5_DIR}/02_checkv"
VIBRANT_OUT="${PHASE5_DIR}/03_vibrant"
CONFIRMED_VIRAL_OUT="${PHASE5_DIR}/04_confirmed_viral"
DRAMV_OUT="${PHASE5_DIR}/05_dramv"
PROPHAGE_FLAG_OUT="${PHASE5_DIR}/06_prophage_flagging"
VOTU_OUT="${PHASE5_DIR}/07_votu_clustering"
PHAGCN_OUT="${PHASE5_DIR}/08_phagcn2"
PROPAGATE_OUT="${PHASE5_DIR}/09_propagate"
CRISPR_OUT="${PHASE5_DIR}/10_crispr_host_linkage"
ARG_OUT="${PHASE5_DIR}/11_prophage_args"
LOG_DIR="${PHASE5_DIR}/logs"

### ---- DATABASES ------------------------------------------------------------
# NOTE: all databases live under /disk252/qamar/databases, NOT under
# BASE_DIR (metagenome_analysis) -- confirmed from actual install paths used
# throughout setup (CheckV, DRAM-v, VOGDB, CARD all installed here).
DATABASES_DIR="/disk252/qamar/databases"

CHECKV_DB="${DATABASES_DIR}/checkv/checkv-db-v1.5"                       # confirmed 6.4G, working (checkv end_to_end tested successfully)
VOGDB_HMM="${DATABASES_DIR}/dramv_config/vog_latest_hmms.txt"            # reused from DRAM-v setup, already hmmpressed -- no separate download
CARD_DB="${DATABASES_DIR}/card/card.json"                                # extracted from card-data.tar.bz2
DRAMV_CONFIG="${DATABASES_DIR}/dramv_config"                             # fully configured -- print_config confirms all paths populated
PHAGCN_DB=""   # not used -- PhaGCN2.1 dropped in favor of PhaBOX2 (confirmed installed/working separately)

### ---- CONDA ENVIRONMENTS (one per tool -- see note below) ----------------
# VirSorter2 and CheckV are commonly bundled since CheckV is VirSorter2's
# companion QC tool and shares compatible dependency versions.
ENV_VIRSORTER2="virsorter2_env"
ENV_CHECKV="virsorter2_env"
ENV_VIBRANT="vibrant_env"
ENV_DRAMV="dramv_env"
ENV_PHAGCN="phagcn_env"
ENV_MINCED="vibrant_env"      # MinCED is a single jar -- no need for its own env
ENV_PROPAGATE="propagate_env"
ENV_BLAST="vibrant_env"       # blastn is a light dependency, fine to reuse

### ---- BATCHING (see METHODS_NOTES.md) ------------------------------------
# VirSorter2/VIBRANT score each CONTIG independently -- they have no concept
# of genome boundaries -- so running N MAGs concatenated into one FASTA is
# equivalent to N separate runs, provided contig headers stay unique and
# traceable. Per-MAG invocation wasted ~17 min of snakemake/conda orchestration
# on a 34-contig genome; batching amortises that across ~BATCH_SIZE genomes.
BATCH_SIZE="${BATCH_SIZE:-200}"     # MAGs per batch FASTA
BATCH_JOBS="${BATCH_JOBS:-6}"       # concurrent batches
BATCH_THREADS="${BATCH_THREADS:-8}" # threads per batch (BATCH_JOBS*BATCH_THREADS <= THREADS)

BATCH_DIR="${PHASE5_DIR}/00_batches"
MAG_INDEX_TSV="${BATCH_DIR}/mag_index.tsv"   # mag_index | mag_id | batch_id | source_path

### ---- RESOURCES (matches Phase 3 THREADS/PARALLEL_JOBS convention) -------
THREADS="${THREADS:-48}"            # threads for single big multi-threaded tools (GTDB-Tk, MMseqs2, DIAMOND)
PARALLEL_JOBS=24      # concurrent MAGs for per-genome tools (MinCED, Prodigal)
PARALLEL_THREADS=2    # threads each per-genome job gets

### ---- Helper: activate a conda env non-interactively (same as Phase 3) ---
activate_env () {
    local env_name="$1"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

### ---- Helper: timestamped logging (same as Phase 3) -----------------------
log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

### ---- SETUP ------------------------------------------------------------
mkdir -p "${BATCH_DIR}" "${VS2_OUT}" "${CHECKV_OUT}" "${VIBRANT_OUT}" "${CONFIRMED_VIRAL_OUT}" \
         "${DRAMV_OUT}" "${PROPHAGE_FLAG_OUT}" "${VOTU_OUT}" "${PHAGCN_OUT}" \
         "${PROPAGATE_OUT}" "${CRISPR_OUT}" "${ARG_OUT}" "${LOG_DIR}"