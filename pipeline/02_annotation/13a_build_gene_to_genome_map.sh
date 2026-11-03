#!/usr/bin/env bash
# ============================================================
# 13a_build_gene_to_genome_map.sh
#
# eggNOG/dbCAN/VFDB were run on the 37.1M-representative non-redundant gene
# catalog (MMseqs2 output), not directly per genome. To attribute a given
# genome's functional annotation, we need to know which MMseqs2 cluster
# representative each of ITS OWN genes belongs to.
#
# Rather than assume a header-naming convention for the concatenated gene
# catalog (risky — I can't verify whether genome names were prefixed during
# concatenation), this rebuilds the mapping from the SOURCE OF TRUTH: each
# genome's own Prodigal .faa file, where we always know which genome we're
# reading regardless of what the header itself says.
#
# Output: gene_to_genome.tsv  (gene_id <tab> genome)   ~46.9M rows
# ============================================================
set -euo pipefail
GENES_DIR=/disk252/qamar/metagenome_analysis/annotation/genes
OUT=/disk252/qamar/metagenome_analysis/annotation/gene_to_genome.tsv
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

# ---- Sanity check: confirm this matches the known total gene count ----
n=$(wc -l < "${OUT}")
if [[ "$n" -ne 46918358 ]]; then
    echo "WARNING: expected 46,918,358 rows (from Prodigal QC earlier), got ${n}."
    echo "  Check for genomes with 0 predicted genes or a partial/interrupted run."
fi

# ---- Diagnostic: check whether the MMseqs2 catalog's own headers already
# carry genome-identifying info, so 13b can pick the right join strategy ----
echo ""
echo "=== Header format check (informs 13b's join strategy) ==="
head -3 /disk252/qamar/metagenome_analysis/annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta | grep '^>'
echo ""
echo "Compare the above to a sample from gene_to_genome.tsv:"
shuf -n 3 "${OUT}"
echo ""
echo "If the gene catalog headers look IDENTICAL in format to the first column"
echo "of gene_to_genome.tsv (same exact strings), the join in 13b will work by"
echo "direct exact match. If the catalog headers look different (e.g. have a"
echo "genome name prefixed that gene_to_genome.tsv's IDs don't have, or vice"
echo "versa), note the difference — 13b has a fallback but flag this now."
