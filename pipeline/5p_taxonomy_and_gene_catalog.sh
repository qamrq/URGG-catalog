#!/usr/bin/env bash
#                 ############################
#               ##                            ##
#             ##  ####                    ####  ##
#           ##   ######                  ######   ##
#         ##    ########                ########    ##
#       ##      ##########            ##########      ##
#     ##        ####  ####################  ####        ##
#   ##          #####  ##################  #####          ##
# ##            ################################            ##
#   ##            ############################            ##
#     ##              ####################              ##
#       ##             ##  ##########  ##             ##
#         ##           ##   ########   ##           ##
#           ##          #    ######    #          ##
#             ##              ####              ##
#               ##                            ##
#                 ############################
#
# ==========================================================
#  PROJECT : URGG catalog
#  SCRIPT  : 5p_taxonomy_and_gene_catalog.sh  -  GTDB-Tk taxonomy, gene prediction and non-redundant gene catalog
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_95_species_final/dereplicated_genomes/*.fa
#  OUTPUT  : annotation/gtdbtk_output_v232/
#            annotation/genes/
#            annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#            annotation/gene_to_genome.tsv
#  USAGE   : bash 5p_taxonomy_and_gene_catalog.sh [gtdbtk prodigal catalog gene_map]   (no argument: all steps in this order)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

step_gtdbtk() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/01_gtdbtk_$(date +%Y%m%d_%H%M%S).log") 2>&1

export GTDBTK_DATA_PATH="${GTDBTK_DB}"

log "Using GTDB-Tk: $(gtdbtk --version | head -n1)"
log "Using GTDB reference DB: ${GTDBTK_DATA_PATH}"

OUT=${ANNOT}/gtdbtk_output_v232
mkdir -p "${OUT}"

log "Starting GTDB-Tk classify_wf on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes"

gtdbtk classify_wf \
    --genome_dir "${GENOME_DIR}" \
    --out_dir "${OUT}" \
    --extension "${GENOME_EXT}" \
    --cpus "${THREADS}" \
    --pplacer_cpus 8

log "GTDB-Tk complete."

MERGED=${OUT}/gtdbtk.all.summary.tsv
head -n1 "${OUT}/gtdbtk.bac120.summary.tsv" > "${MERGED}"
for f in "${OUT}"/gtdbtk.bac120.summary.tsv "${OUT}"/gtdbtk.ar53.summary.tsv; do
    [[ -f "$f" ]] && tail -n +2 "$f" >> "${MERGED}"
done
log "Combined taxonomy written to ${MERGED}"
)

step_prodigal() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/02_prodigal_$(date +%Y%m%d_%H%M%S).log") 2>&1

predict_one () {
    local f="$1" GENES="$2"
    local base; base=$(basename "$f" ."${GENOME_EXT:-fa}")
    local faa="${GENES}/${base}.faa"
    local gff="${GENES}/${base}.gff"
    [[ -s "$faa" ]] && return 0
    if ! prodigal -i "$f" -a "$faa" -o "$gff" -f gff -p single -q 2>/dev/null; then
        prodigal -i "$f" -a "$faa" -o "$gff" -f gff -p meta -q
    fi
}
export -f predict_one
export GENOME_EXT

log "Predicting genes on $(find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | wc -l) genomes with ${PARALLEL_JOBS} jobs"

find "${GENOME_DIR}" -maxdepth 1 -name "*.${GENOME_EXT}" | \
    parallel -j "${PARALLEL_JOBS}" --joblog "${LOGS}/02_prodigal.joblog" --bar \
    predict_one {} "${GENES}"

log "Prodigal complete. $(ls "${GENES}"/*.faa | wc -l) protein files in ${GENES}"
)

step_catalog() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/03_mmseqs2_$(date +%Y%m%d_%H%M%S).log") 2>&1

ALL_PROT=${CATALOG}/all_proteins.faa
TMP=${CATALOG}/tmp
mkdir -p "${TMP}"

log "Concatenating all predicted proteins..."
: > "${ALL_PROT}"
for faa in "${GENES}"/*.faa; do
    base=$(basename "$faa" .faa)
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
)

step_gene_map() (
set -euo pipefail
GENES_DIR=${URGG_PROJECT}/annotation/genes
OUT=${URGG_PROJECT}/annotation/gene_to_genome.tsv
THREADS=24

echo "Building gene->genome map from $(find "$GENES_DIR" -maxdepth 1 -name "*.faa" | wc -l) .faa files..."

extract_one () {
    local f="$1"
    local genome; genome=$(basename "$f" .faa)
    grep '^>' "$f" | sed -e 's/^>//' -e 's/ .*//' | \
        awk -v g="$genome" -v OFS='\t' '{print $0, g}'
}
export -f extract_one

find "$GENES_DIR" -maxdepth 1 -name "*.faa" | \
    parallel -j "${THREADS}" extract_one {} > "${OUT}"

echo "Done: $(wc -l < "${OUT}") gene->genome rows written to ${OUT}"

n=$(wc -l < "${OUT}")
if [[ "$n" -ne 46918358 ]]; then
    echo "WARNING: expected 46,918,358 rows (from Prodigal QC earlier), got ${n}."
    echo "  Check for genomes with 0 predicted genes or a partial/interrupted run."
fi

echo ""
echo "=== Header format check ==="
head -3 ${URGG_PROJECT}/annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta | grep '^>'
echo ""
echo "Compare the above to a sample from gene_to_genome.tsv:"
shuf -n 3 "${OUT}"
echo ""
echo "If the gene catalog headers look IDENTICAL in format to the first column"
echo "of gene_to_genome.tsv (same exact strings), downstream joins will work by"
echo "direct exact match. If the catalog headers look different (e.g. have a"
echo "genome name prefixed that gene_to_genome.tsv's IDs don't have, or vice"
echo "versa), note the difference before running the annotation steps."
)

for s in ${@:-gtdbtk prodigal catalog gene_map}; do
  echo "==> 5p_taxonomy_and_gene_catalog.sh: $s"
  "step_$s"
done
