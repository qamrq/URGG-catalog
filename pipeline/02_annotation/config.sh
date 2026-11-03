#!/usr/bin/env bash
# ============================================================
# config.sh  —  Shared configuration for Phase 3 annotation
# Source this at the top of every step script:  source ./config.sh
# Edit paths / thread counts here ONCE; all scripts inherit them.
# ============================================================

# ---- Core project paths ----
export PROJECT=/disk252/qamar/metagenome_analysis
export DREP_DIR=${PROJECT}/drep_analysis
export METADATA=${DREP_DIR}/metadata.tsv                            # genome|bioproject|host|species|sample_site|...
export GENOME_EXT=fa                                               # extension of genome files

# *** PHASE 3 ANNOTATION USES THE 95% ANI SPECIES CATALOG — NOT the 99% strain catalog ***
# Per the dRep export doc (Part 5): 25,473 species reps are the primary input for
# taxonomy / gene prediction / functional annotation / abundance mapping.
# The 99% ANI strain catalog (23,999 genomes) is reserved ONLY for the SGVFinder2
# SV reference database (separate pipeline, not part of these Phase 3 scripts).
export GENOME_DIR=${DREP_DIR}/drep_95_species_final/dereplicated_genomes   # 25,473 *.fa species reps

# Kept for reference / used only by the future SV pipeline — do NOT point Phase 3 here.
export SV_STRAIN_DIR=${DREP_DIR}/drep_99_strain/dereplicated_genomes       # 23,999 *.fa strain reps (SV reference ONLY)

# ---- Output roots ----
export ANNOT=${PROJECT}/annotation
export GENES=${ANNOT}/genes                 # Prodigal .faa/.gff per genome
export CATALOG=${ANNOT}/gene_catalog        # MMseqs2 non-redundant catalog
export LOGS=${ANNOT}/logs                   # per-step logs
export DB=/disk252/qamar/databases          # all reference databases live here

# ---- Reference database paths (populated by 00b_download_databases.sh) ----
# GTDB-Tk DB: using EXISTING R226 reference package on disk212 — confirmed
# compatible with the installed GTDB-Tk v2.5.2. Do NOT re-download via
# 00b_download_databases.sh (that block is now redundant for this tool).
export GTDBTK_DB=/disk252/qamar/databases/gtdbtk_database/release232/
export NR_DMND=${DB}/ncbi_nr/nr.dmnd
export VFDB_DMND=${DB}/vfdb/vfdb.dmnd
export EGGNOG_DB=${DB}/eggnog
export DBCAN_DB=${DB}/dbcan
# AMRFinderPlus DB: `amrfinder -u` IGNORES the AMRFINDER_DB env var and always
# installs to its own default location inside the conda env
# (<env>/share/amrfinderplus/data/<version>/). Point this at that path so
# 05_amrfinderplus.sh's per-genome calls use it explicitly and don't rely on
# whatever the tool picks as "latest" at run time.
# AMRFinderPlus DB: `amrfinder -u` IGNORES the AMRFINDER_DB env var and always
# installs to its own default location inside the conda env
# (<env>/share/amrfinderplus/data/<version>/), and maintains a `latest`
# symlink pointing at whatever version is current. Point at that symlink
# rather than a hardcoded date string so future `amrfinder -u` runs don't
# require editing this file again.
export AMRFINDER_DB=/home/qamar/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest
export PLATON_DB=${DB}/platon/db
export CHECKV_DB=${DB}/checkv/checkv-db-v1.5

# ---- Compute settings ----
# Both dRep catalogs (95% species, 99% strain) are complete as of this export —
# no competing dRep job should be running. Still worth an idle `top -u qamar`
# check before launching GTDB-Tk (heaviest single step) in case anything else
# is queued on the box.
export THREADS=48            # threads for single big multi-threaded tools (GTDB-Tk, MMseqs2, DIAMOND)
export PARALLEL_JOBS=24      # concurrent genomes for per-genome tools (AMRFinder, Platon, PhageBoost)
export PARALLEL_THREADS=2    # threads each per-genome job gets (PARALLEL_JOBS * PARALLEL_THREADS <= THREADS-ish)

# ---- Conda envs (created by 00_install_tools.sh) ----
export ENV_DREP=drep         # holds prodigal, mmseqs2, diamond, parallel, pigz (lightweight standalone)

# GTDB-Tk: using your EXISTING install on disk192 (v2.5.2, matched to R226 DB
# on disk212) — NOT created by 00_install_tools.sh. Activation differs from
# other envs (different miniconda base), so GTDB-Tk is activated directly in
# 01_gtdbtk.sh via its own `source .../activate <path>` line, bypassing the
# activate_env() helper used for everything else.
#export GTDBTK_CONDA_BASE=/disk192/qamar/miniconda3
export GTDBTK_ENV_PATH=/disk192/qamar/miniconda3/envs/gtdbtk-2.5.2

export ENV_EGGNOG=eggnog
export ENV_DBCAN=dbcan
export ENV_AMR=amrfinder
export ENV_PLATON=platon
export ENV_CHECKV=checkv
export ENV_PHAGEBOOST=phageboost

# ---- Helper: activate a conda env non-interactively ----
activate_env () {
    local env_name="$1"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

# ---- Helper: timestamped logging ----
log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---- Create output dirs up front ----
mkdir -p "${ANNOT}" "${GENES}" "${CATALOG}" "${LOGS}" "${DB}"