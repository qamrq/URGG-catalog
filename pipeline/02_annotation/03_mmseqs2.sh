#!/usr/bin/env bash
# ============================================================
# 03_mmseqs2.sh  —  Build non-redundant gene catalog
# Env: drep | Input: all GENES/*.faa | Output: CATALOG/nr_gene_catalog_rep_seq.fasta
# Clusters at 95% identity / 90% coverage (standard for gene catalogs).
# The rep_seq FASTA is the input for all DIAMOND / eggNOG / dbCAN steps.
# ============================================================
set -euo pipefail
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/03_mmseqs2_$(date +%Y%m%d_%H%M%S).log") 2>&1
activate_env "${ENV_DREP}"

ALL_PROT=${CATALOG}/all_proteins.faa
TMP=${CATALOG}/tmp
mkdir -p "${TMP}"

log "Concatenating all predicted proteins..."
# stream-concat to avoid a huge argv; strip nothing, keep headers unique per genome
: > "${ALL_PROT}"
for faa in "${GENES}"/*.faa; do
    base=$(basename "$faa" .faa)
    # prefix each header with its genome so catalog members trace back to source
    sed "s/^>/>${base}~/" "$faa" >> "${ALL_PROT}"
done
log "Total proteins: $(grep -c '^>' "${ALL_PROT}")"

log "Running MMseqs2 easy-linclust (95% id, 90% cov)..."
mmseqs easy-linclust \
    "${ALL_PROT}" \
    "${CATALOG}/nr_gene_catalog" \
    "${TMP}" \
    --min-seq-id 0.95 -c 0.9 --cov-mode 1 \
    --threads "${THREADS}"

log "Non-redundant catalog: ${CATALOG}/nr_gene_catalog_rep_seq.fasta"
log "Representatives: $(grep -c '^>' "${CATALOG}/nr_gene_catalog_rep_seq.fasta")"
log "Cluster membership: ${CATALOG}/nr_gene_catalog_cluster.tsv"
rm -rf "${TMP}"
