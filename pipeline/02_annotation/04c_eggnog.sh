#!/usr/bin/env bash
# ============================================================
# 04c_eggnog.sh  —  eggNOG-mapper functional annotation
# Env: eggnog | Input: CATALOG/nr_gene_catalog_rep_seq.fasta
#
# This is your KEGG/KO route. eggNOG-mapper outputs KO numbers, KEGG
# pathways/modules, COG categories, GO terms, EC numbers and CAZy families
# in one pass. (CowPI from the proposal is a 16S-amplicon PICRUSt tool, not
# valid on MAG proteins — eggNOG replaces it for KO/KEGG. If you want a
# dedicated KEGG-ortholog caller, KofamScan is the free alternative.)
# ============================================================
# Temporarily disable strict unbound variable checking for Conda
set +u
activate_env "${ENV_EGGNOG}"
set -u
source "/disk252/qamar/metagenome_analysis/phase_3_analysis/scripts/config.sh"
exec > >(tee -a "${LOGS}/04c_eggnog_$(date +%Y%m%d_%H%M%S).log") 2>&1
#activate_env "${ENV_EGGNOG}"

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/eggnog_output
mkdir -p "${OUT}"

log "Running eggNOG-mapper (DIAMOND mode)..."
emapper.py \
    -i "${QUERY}" \
    --itype proteins \
    -m diamond \
    --data_dir "${EGGNOG_DB}" \
    -o nr_catalog \
    --output_dir "${OUT}" \
    --cpu "${THREADS}" \
    --override

log "eggNOG-mapper done."
log "Annotations: ${OUT}/nr_catalog.emapper.annotations"
log "  -> KO column = KEGG orthologs; KEGG_Pathway / KEGG_Module = pathways; CAZy column present too."
