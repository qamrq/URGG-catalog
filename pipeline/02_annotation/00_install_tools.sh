#!/usr/bin/env bash
# ============================================================
# 00_install_tools.sh  —  One-shot installation of all Phase 3 tools
#
# STRATEGY (important — read the comment):
#   Lightweight standalone tools (Prodigal, MMseqs2, DIAMOND, GNU parallel,
#   pigz) are added to your EXISTING `drep` env, as requested.
#
#   Tools with heavy / conflicting dependencies get their OWN env. Putting
#   GTDB-Tk, dbCAN, eggNOG-mapper, Platon, CheckV and (especially) PhageBoost
#   into `drep` would downgrade numpy/hmmer/prodigal/scikit-learn and break
#   dRep + each other. Isolated envs is the standard, reproducible approach.
#
# Idempotent: re-running skips envs that already exist.
# Run from the scripts directory:  bash 00_install_tools.sh
# ============================================================
set -euo pipefail
source ./config.sh

source "$(conda info --base)/etc/profile.d/conda.sh"

env_exists () { conda env list | awk '{print $1}' | grep -qx "$1"; }

log "=== Phase 3 tool installation starting ==="

# ------------------------------------------------------------
# 1) Lightweight standalone tools -> existing `drep` env
# ------------------------------------------------------------
log "[1/8] Adding Prodigal, MMseqs2, DIAMOND, parallel, pigz to '${ENV_DREP}'"
conda activate "${ENV_DREP}"
conda install -y -c bioconda -c conda-forge \
    prodigal \
    "mmseqs2>=15" \
    "diamond>=2.1" \
    parallel \
    pigz
conda deactivate

# ------------------------------------------------------------
# 2) GTDB-Tk — SKIPPED: reusing existing install
#    v2.5.2 at /disk192/qamar/miniconda3/envs/gtdbtk, paired with the R226
#    reference DB at /disk212/qamar/databases/release226/. Confirmed
#    compatible version pairing. See config.sh (GTDBTK_CONDA_BASE,
#    GTDBTK_ENV_PATH, GTDBTK_DB) and 01_gtdbtk.sh.
# ------------------------------------------------------------
log "[2/8] GTDB-Tk — skipping install, reusing existing env at ${GTDBTK_ENV_PATH:-<not set>}"

# ------------------------------------------------------------
# 3) eggNOG-mapper (isolated) — provides KO/KEGG, COG, GO, CAZy, etc.
# ------------------------------------------------------------
if env_exists "${ENV_EGGNOG}"; then
    log "[3/8] '${ENV_EGGNOG}' already exists — skipping"
else
    log "[3/8] Creating '${ENV_EGGNOG}' (eggNOG-mapper)"
    conda create -y -n "${ENV_EGGNOG}" -c bioconda -c conda-forge "eggnog-mapper>=2.1.12"
fi

# ------------------------------------------------------------
# 4) dbCAN (isolated) — CAZyme annotation
# ------------------------------------------------------------
if env_exists "${ENV_DBCAN}"; then
    log "[4/8] '${ENV_DBCAN}' already exists — skipping"
else
    log "[4/8] Creating '${ENV_DBCAN}' (dbCAN / run_dbcan)"
    conda create -y -n "${ENV_DBCAN}" -c bioconda -c conda-forge "dbcan>=4.1"
fi

# ------------------------------------------------------------
# 5) AMRFinderPlus (isolated)
# ------------------------------------------------------------
if env_exists "${ENV_AMR}"; then
    log "[5/8] '${ENV_AMR}' already exists — skipping"
else
    log "[5/8] Creating '${ENV_AMR}' (NCBI AMRFinderPlus)"
    conda create -y -n "${ENV_AMR}" -c bioconda -c conda-forge ncbi-amrfinderplus
fi

# ------------------------------------------------------------
# 6) Platon (isolated) — plasmid detection
# ------------------------------------------------------------
if env_exists "${ENV_PLATON}"; then
    log "[6/8] '${ENV_PLATON}' already exists — skipping"
else
    log "[6/8] Creating '${ENV_PLATON}' (Platon)"
    conda create -y -n "${ENV_PLATON}" -c bioconda -c conda-forge platon
fi

# ------------------------------------------------------------
# 7) CheckV (isolated) — viral genome quality
# ------------------------------------------------------------
if env_exists "${ENV_CHECKV}"; then
    log "[7/8] '${ENV_CHECKV}' already exists — skipping"
else
    log "[7/8] Creating '${ENV_CHECKV}' (CheckV)"
    conda create -y -n "${ENV_CHECKV}" -c bioconda -c conda-forge checkv
fi

# ------------------------------------------------------------
# 8) PhageBoost (isolated — pins old numpy/scikit-learn; MUST be separate)
# ------------------------------------------------------------
if env_exists "${ENV_PHAGEBOOST}"; then
    log "[8/8] '${ENV_PHAGEBOOST}' already exists — skipping"
else
    log "[8/8] Creating '${ENV_PHAGEBOOST}' (PhageBoost)"
    # PhageBoost is a pip package with strict deps; give it a clean py3.7 env.
    conda create -y -n "${ENV_PHAGEBOOST}" -c conda-forge python=3.7 pip
    conda activate "${ENV_PHAGEBOOST}"

    # --- Known packaging bug workaround ---
    # PhageBoost pulls in the latest pyrodigal (currently 3.7.1) unpinned.
    # pyrodigal 3.7.1's pyproject.toml declares `license = "GPL-3.0-or-later"`
    # as a plain string, which newer scikit-build-core (>=0.10) rejects under
    # PEP 639 (expects a table, not a string) — build fails with:
    #   "Field 'project.license' has an invalid type... (got 'GPL-3.0-or-later')"
    # This is a pyrodigal/scikit-build-core version mismatch, not an env issue.
    # Fix: pre-install an older pyrodigal (pre-dates the strict license field
    # requirement, and still has prebuilt wheels for cp37) so pip reuses it
    # instead of building the newest one from source.
    pip install "pyrodigal<3.0" || pip install "pyrodigal==2.3.0"

    pip install phageboost
    conda deactivate
fi

log "=== Tool installation complete ==="
log "Next: run 00b_download_databases.sh to fetch reference databases (large, one-time)."
log "Verify with:  conda env list"