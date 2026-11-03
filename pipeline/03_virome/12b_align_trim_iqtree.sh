#!/usr/bin/env bash
###############################################################################
# 12b_align_trim_iqtree.sh — Phase 5: align, trim, and build ML tree
#
# Follows the exact methodology from the reference paper's phylogenetics
# section, substituting MAFFT for GTDB-Tk's alignment step (GTDB-Tk itself
# is bacterial/archaeal-marker-specific and doesn't apply to viral markers,
# but the downstream trimAl + IQ-TREE steps are identical):
#   1. MAFFT alignment of terminase large subunit protein sequences
#   2. trimAl (v1.4, default parameters) -- EXACT match to reference methods
#   3. IQ-TREE (-m MFP -bb 1000 -redo -mredo) -- EXACT match to reference
#      methods: MFP = auto model selection, -bb 1000 = ultrafast bootstrap
#      with 1000 replicates
#
# *** VERIFY BEFORE RUNNING ***
#   which mafft trimal iqtree2   # or `iqtree` depending on install
# Install if missing:
#   conda install -y -c bioconda mafft trimal iqtree
#
# Runtime note: ~9,000 vOTUs is a real tree size for IQ-TREE with 1000
# ultrafast bootstraps -- this can take hours depending on how many vOTUs
# actually have a terminase hit (likely a few thousand, not all 8,980,
# since only Caudoviricetes-type phages carry this specific marker). Check
# terminase_sequences.faa's sequence count from 12a before committing to a
# long unattended run.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

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

###############################################################################
# 1. MAFFT alignment
###############################################################################
ALIGNED_FAA="${PHYLO_DIR}/terminase_aligned.faa"
if [[ ! -s "${ALIGNED_FAA}" ]]; then
    log "Running MAFFT alignment (this may take a while at ${n_seqs} sequences)..."
    mafft --auto --thread "${THREADS}" "${INPUT_FAA}" > "${ALIGNED_FAA}" 2> "${LOG_DIR}/mafft.log"
else
    log "Existing alignment found -- skipping (delete ${ALIGNED_FAA} to force a rerun)"
fi
log "  alignment: $(grep -c '^>' "${ALIGNED_FAA}") sequences"

###############################################################################
# 2. trimAl (v1.4, default parameters -- exact match to reference methods)
###############################################################################
TRIMMED_FAA="${PHYLO_DIR}/terminase_trimmed.faa"
if [[ ! -s "${TRIMMED_FAA}" ]]; then
    log "Running trimAl (default parameters)..."
    trimal -in "${ALIGNED_FAA}" -out "${TRIMMED_FAA}" -automated1 \
        > "${LOG_DIR}/trimal.log" 2>&1
else
    log "Existing trimmed alignment found -- skipping"
fi

###############################################################################
# 3. IQ-TREE (-m MFP -bb 1000 -redo -mredo -- exact match to reference methods)
###############################################################################
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
