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
#  SCRIPT  : 13a_build_gene_to_genome_map.sh  -  Gene-to-genome table for the gene catalog
#  AUTHOR  : QRQ
#  INPUT   : annotation/genes/*.faa
#  OUTPUT  : annotation/gene_to_genome.tsv
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

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
echo "=== Header format check (informs 13b's join strategy) ==="
head -3 ${URGG_PROJECT}/annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta | grep '^>'
echo ""
echo "Compare the above to a sample from gene_to_genome.tsv:"
shuf -n 3 "${OUT}"
echo ""
echo "If the gene catalog headers look IDENTICAL in format to the first column"
echo "of gene_to_genome.tsv (same exact strings), the join in 13b will work by"
echo "direct exact match. If the catalog headers look different (e.g. have a"
echo "genome name prefixed that gene_to_genome.tsv's IDs don't have, or vice"
echo "versa), note the difference — 13b has a fallback but flag this now."
