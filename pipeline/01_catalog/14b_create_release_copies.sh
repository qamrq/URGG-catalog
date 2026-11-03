#!/usr/bin/env bash
# ============================================================
# 14b_create_release_copies.sh
#
# Creates a SEPARATE copy of the genome catalog with clean URGG filenames,
# for the website/public release. The original catalog under
# drep_95_species_final/dereplicated_genomes/ is never touched — every
# existing script and joined table in this pipeline depends on those
# original filenames.
#
# Also relabels each genome's FASTA headers to prefix the public ID onto
# the original contig name (>URGG00001|<original_contig_header>), so a
# downloaded genome is self-describing even in isolation, while the
# original contig identity is still preserved for traceability back to
# your internal records if ever needed.
# ============================================================
set -euo pipefail

DREP=/disk252/qamar/metagenome_analysis/drep_analysis/drep_95_species_final
SRC_DIR="${DREP}/dereplicated_genomes"
MAPPING="${DREP}/genome_id_mapping.tsv"
RELEASE_DIR=/disk252/qamar/metagenome_analysis/website_release/genomes
THREADS=24

mkdir -p "${RELEASE_DIR}"

[[ -f "${MAPPING}" ]] || { echo "ERROR: run 14a_build_urgg_mapping.py first"; exit 1; }

rename_one () {
    local original="$1" public_id="$2" src_dir="$3" out_dir="$4"
    local src="${src_dir}/${original}.fa"
    local dst="${out_dir}/${public_id}.fa"
    if [[ ! -f "$src" ]]; then
        echo "WARNING: source genome not found for ${original} (expected ${public_id})" >&2
        return 1
    fi
    # Prefix each contig header with the public ID, keep original for traceability
    sed "s/^>/>${public_id}|/" "$src" > "$dst"
}
export -f rename_one

echo "Creating renamed release copies..."
tail -n +2 "${MAPPING}" | cut -f1,2 | \
    parallel -j "${THREADS}" --colsep '\t' \
    rename_one {1} {2} "${SRC_DIR}" "${RELEASE_DIR}"

n_release=$(find "${RELEASE_DIR}" -maxdepth 1 -name "*.fa" | wc -l)
n_expected=$(($(wc -l < "${MAPPING}") - 1))
echo ""
echo "Release copies created: ${n_release} (expected ${n_expected})"
if [[ "${n_release}" -ne "${n_expected}" ]]; then
    echo "WARNING: count mismatch — check the WARNING lines above for missing sources"
fi
echo "Release genomes at: ${RELEASE_DIR}"
echo "Original catalog UNCHANGED at: ${SRC_DIR}"
