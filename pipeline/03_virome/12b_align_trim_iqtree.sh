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
#  SCRIPT  : 12b_align_trim_iqtree.sh  -  TerL alignment, trimming and ML tree (MAFFT, trimAl, IQ-TREE)
#  AUTHOR  : QRQ
#  INPUT   : phase_5_virome/12_phylogeny/terminase_sequences_filtered.faa
#  OUTPUT  : phase_5_virome/12_phylogeny/terminase_trimmed.faa, .treefile
# ==========================================================

source "$(dirname "$0")/config.sh"

PHYLO_DIR="${PHASE5_DIR}/12_phylogeny"
INPUT_FAA="${PHYLO_DIR}/terminase_sequences_filtered.faa"

if [[ ! -s "${INPUT_FAA}" ]]; then
    log "ERROR: ${INPUT_FAA} not found or empty. Run 12a_extract_terminase_marker.sh first."
    exit 1
fi

n_seqs=$(grep -c '^>' "${INPUT_FAA}")
log "Input: ${n_seqs} terminase sequences"
if [[ "${n_seqs}" -lt 4 ]]; then
    log "ERROR: fewer than 4 sequences -- cannot build a meaningful tree."
    exit 1
fi

for tool in mafft trimal; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        log "ERROR: ${tool} not found. Install: conda install -y -c bioconda ${tool}"
        exit 1
    fi
done
IQTREE_BIN=""
for candidate in iqtree2 iqtree; do
    if command -v "${candidate}" >/dev/null 2>&1; then
        IQTREE_BIN="${candidate}"
        break
    fi
done
if [[ -z "${IQTREE_BIN}" ]]; then
    log "ERROR: iqtree/iqtree2 not found. Install: conda install -y -c bioconda iqtree"
    exit 1
fi
log "Using IQ-TREE binary: ${IQTREE_BIN}"

ALIGNED_FAA="${PHYLO_DIR}/terminase_aligned.faa"
if [[ ! -s "${ALIGNED_FAA}" ]]; then
    log "Running MAFFT alignment (this may take a while at ${n_seqs} sequences)..."
    mafft --auto --thread "${THREADS}" "${INPUT_FAA}" > "${ALIGNED_FAA}" 2> "${LOG_DIR}/mafft.log"
else
    log "Existing alignment found -- skipping (delete ${ALIGNED_FAA} to force a rerun)"
fi
log "  alignment: $(grep -c '^>' "${ALIGNED_FAA}") sequences"

TRIMMED_FAA="${PHYLO_DIR}/terminase_trimmed.faa"
if [[ ! -s "${TRIMMED_FAA}" ]]; then
    log "Running trimAl (default parameters)..."
    trimal -in "${ALIGNED_FAA}" -out "${TRIMMED_FAA}" -automated1 \
        > "${LOG_DIR}/trimal.log" 2>&1
else
    log "Existing trimmed alignment found -- skipping"
fi

log "Running IQ-TREE (MFP model selection + 1000 ultrafast bootstraps)..."
log "  This is the long-running step -- consider tmux/screen."
"${IQTREE_BIN}" -s "${TRIMMED_FAA}" -m MFP -bb 1000 -nt "${THREADS}" -redo \
    > "${LOG_DIR}/iqtree.log" 2>&1

TREE_FILE="${TRIMMED_FAA}.treefile"
if [[ -s "${TREE_FILE}" ]]; then
    log "IQ-TREE finished -> ${TREE_FILE}"
    log "Next: fig_phylo_tree_real.R (visualization with function rings)"
else
    log "ERROR: IQ-TREE did not produce a .treefile -- check ${LOG_DIR}/iqtree.log"
    exit 1
fi
