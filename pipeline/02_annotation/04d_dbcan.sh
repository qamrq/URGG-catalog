#!/usr/bin/env bash
# ============================================================
# 04d_dbcan.sh  —  CAZyme annotation (run_dbcan)
# Env: dbcan | Input: CATALOG/nr_gene_catalog_rep_seq.fasta
#
# dbCAN gives dedicated, higher-confidence CAZyme calls (HMMER + DIAMOND +
# dbCAN_sub consensus) — important for your fibre-degradation / methane
# functional story. Run on the non-redundant catalog for efficiency.
# ============================================================
set -euo pipefail
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/04d_dbcan_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DBCAN}"

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/dbcan_output
mkdir -p "${OUT}"

log "Running run_dbcan CAZyme_annotation (dbCAN v5.x subcommand syntax)..."
# NOTE: --methods defaults to diamond,hmm,dbCANsub. If dbCAN-sub.hmm hasn't
# finished downloading yet (known server-side transfer issue — see download
# troubleshooting notes), this will fail on the dbCANsub step. Restricting to
# diamond,hmm here so annotation can proceed now; re-run with all three once
# dbCAN-sub.hmm is confirmed complete (grep -c "^NAME" on the file is a quick check).
#DBCAN_METHODS=${DBCAN_METHODS:-diamond,hmm}
run_dbcan CAZyme_annotation \
    --input_raw_data "${QUERY}" \
    --mode protein \
    --output_dir "${OUT}" \
    --db_dir "${DBCAN_DB}" \
    --threads "${THREADS}"

log "dbCAN done. Consensus CAZyme table: ${OUT}/overview.txt"