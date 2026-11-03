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
#  SCRIPT  : 14b_create_release_copies.sh  -  Copy species genomes under URGG file names
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis/drep_95_species_final/dereplicated_genomes, genome_id_mapping.tsv
#  OUTPUT  : website_release/genomes/URGG*.fa
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

DREP=${URGG_PROJECT}/drep_analysis/drep_95_species_final
SRC_DIR="${DREP}/dereplicated_genomes"
MAPPING="${DREP}/genome_id_mapping.tsv"
RELEASE_DIR=${URGG_PROJECT}/website_release/genomes
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
