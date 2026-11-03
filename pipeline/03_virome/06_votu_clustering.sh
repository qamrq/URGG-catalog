#!/usr/bin/env bash
###############################################################################
# 06_votu_clustering.sh — Phase 5: Virome pipeline, Step 6
#
# Clusters the confirmed viral set into vOTUs at 95% ANI / 85% AF (MIUViG
# standard, and the threshold the source protocol uses for this specific
# step -- see METHODS_NOTES.md §2 for why this differs from the 99% ANI
# MAG catalogue used for viral calling itself).
#
# METHOD: CheckV's own recommended pipeline (all-vs-all blastn -> anicalc.py
# -> aniclust.py), NOT a MAG-dereplication tool like dRep -- vOTU clustering
# in the virus literature is conventionally done this way because ANI between
# partial/fragmented viral contigs needs alignment-fraction-aware clustering,
# which dRep's genome-completeness assumptions aren't built for.
#
# Input : CONFIRMED_VIRAL_OUT/confirmed_viral.fa  (15,342 seqs, from 03)
# Output: VOTU_OUT/votu_representatives.fa         one representative per vOTU
#         VOTU_OUT/votu_clusters.tsv               cluster_id | representative | members
#         VOTU_OUT/vOTU_index.tsv                  votu_id | representative | n_members | member_mag_indices
#
# This is the UNIFIED catalog decision from earlier discussion: one vOTU
# clustering across ALL host species combined, not per-species -- required
# for the cross-host presence/absence comparison to have a shared coordinate
# system at all.
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"
activate_env "${ENV_VIRSORTER2}"

log "Step 6: vOTU clustering (95% ANI / 85% AF) on confirmed viral set"

CONFIRMED_FA="${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa"
if [[ ! -s "${CONFIRMED_FA}" ]]; then
    log "ERROR: ${CONFIRMED_FA} not found or empty. Run 03_merge_checkv.sh first."
    exit 1
fi
n_input=$(grep -c '^>' "${CONFIRMED_FA}")
log "Input: ${n_input} confirmed viral sequences"

###############################################################################
# 1. Locate (or fetch) CheckV's clustering helper scripts.
#    These ship in the CheckV GitHub repo under scripts/, not always bundled
#    with the conda package -- search the env first, download only if missing.
###############################################################################
ANICALC=$(find "$(conda info --base)/envs/${ENV_VIRSORTER2}" -name "anicalc.py" 2>/dev/null | head -1)
ANICLUST=$(find "$(conda info --base)/envs/${ENV_VIRSORTER2}" -name "aniclust.py" 2>/dev/null | head -1)

if [[ -z "${ANICALC}" || -z "${ANICLUST}" ]]; then
    log "anicalc.py/aniclust.py not found in ${ENV_VIRSORTER2} -- fetching from CheckV GitHub"
    mkdir -p "${VOTU_OUT}/_checkv_scripts"
    if ! command -v git >/dev/null 2>&1; then
        log "  git not in this env -- installing"
        conda install -y -c conda-forge git
    fi
    if [[ ! -f "${VOTU_OUT}/_checkv_scripts/anicalc.py" ]]; then
        git clone --depth 1 https://github.com/snayfach/CheckV.git "${VOTU_OUT}/_checkv_scripts/CheckV_repo" \
            > "${LOG_DIR}/votu_checkv_scripts_clone.log" 2>&1
        cp "${VOTU_OUT}/_checkv_scripts/CheckV_repo/scripts/anicalc.py"  "${VOTU_OUT}/_checkv_scripts/"
        cp "${VOTU_OUT}/_checkv_scripts/CheckV_repo/scripts/aniclust.py" "${VOTU_OUT}/_checkv_scripts/"
    fi
    ANICALC="${VOTU_OUT}/_checkv_scripts/anicalc.py"
    ANICLUST="${VOTU_OUT}/_checkv_scripts/aniclust.py"
fi

if [[ ! -f "${ANICALC}" || ! -f "${ANICLUST}" ]]; then
    log "ERROR: could not locate or fetch anicalc.py/aniclust.py."
    log "       Manually place them at ${VOTU_OUT}/_checkv_scripts/ and rerun."
    log "       Source: https://github.com/snayfach/CheckV/tree/master/scripts"
    exit 1
fi
log "Using anicalc: ${ANICALC}"
log "Using aniclust: ${ANICLUST}"

###############################################################################
# 2. All-vs-all blastn (megablast) -- CheckV's recommended parameters
###############################################################################
BLAST_DB="${VOTU_OUT}/confirmed_viral_blastdb"
BLAST_OUT="${VOTU_OUT}/blast_all_vs_all.tsv"

if [[ ! -s "${BLAST_OUT}" ]]; then
    log "Building BLAST db..."
    makeblastdb -in "${CONFIRMED_FA}" -dbtype nucl -out "${BLAST_DB}" \
        > "${LOG_DIR}/votu_makeblastdb.log" 2>&1

    log "Running all-vs-all megablast (this is O(n^2)-ish; ${n_input} seqs -- watch runtime)..."
    blastn -query "${CONFIRMED_FA}" -db "${BLAST_DB}" \
        -outfmt '6 std qlen slen' \
        -max_target_seqs 10000 \
        -perc_identity 90 \
        -num_threads "${THREADS}" \
        -out "${BLAST_OUT}" \
        > "${LOG_DIR}/votu_blastn.log" 2>&1
else
    log "Existing BLAST output found -- skipping (delete ${BLAST_OUT} to force a rerun)"
fi

n_pairs=$(wc -l < "${BLAST_OUT}")
log "  ${n_pairs} pairwise hits"

###############################################################################
# 3. anicalc.py -- converts raw blast hits into per-genome-pair ANI + AF,
#    correctly handling multiple HSPs per pair (a single blastn line per
#    alignment segment is NOT the same as ANI between the two full sequences).
###############################################################################
log "Computing pairwise ANI/AF (anicalc.py)..."
ANI_OUT="${VOTU_OUT}/ani.tsv"
python "${ANICALC}" -i "${BLAST_OUT}" -o "${ANI_OUT}" \
    > "${LOG_DIR}/votu_anicalc.log" 2>&1

if [[ ! -s "${ANI_OUT}" ]]; then
    log "ERROR: anicalc.py produced no output -- see ${LOG_DIR}/votu_anicalc.log"
    exit 1
fi

###############################################################################
# 4. aniclust.py -- greedy centroid clustering at 95% ANI / 85% AF (the
#    coverage of the SHORTER/target sequence, per CheckV's own vOTU
#    recommendation and MIUViG standards for viral species-level clusters).
###############################################################################
log "Clustering into vOTUs (95% ANI, 85% AF)..."
CLUSTERS_OUT="${VOTU_OUT}/votu_clusters.tsv"
python "${ANICLUST}" \
    --fna "${CONFIRMED_FA}" \
    --ani "${ANI_OUT}" \
    --out "${CLUSTERS_OUT}" \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0 \
    > "${LOG_DIR}/votu_aniclust.log" 2>&1

if [[ ! -s "${CLUSTERS_OUT}" ]]; then
    log "ERROR: aniclust.py produced no output -- see ${LOG_DIR}/votu_aniclust.log"
    exit 1
fi

n_votus=$(wc -l < "${CLUSTERS_OUT}")
log "  ${n_votus} vOTUs from ${n_input} input sequences"

###############################################################################
# 5. Extract representative sequences (first column of each cluster line)
###############################################################################
log "Extracting vOTU representative sequences..."
cut -f1 "${CLUSTERS_OUT}" | sort -u > "${VOTU_OUT}/.rep_ids.txt"
sed 's/^/^/' "${VOTU_OUT}/.rep_ids.txt" > "${VOTU_OUT}/.rep_ids_regex.txt"

if ! command -v seqkit >/dev/null 2>&1; then
    log "ERROR: seqkit required. conda install -n ${ENV_VIRSORTER2} -c bioconda seqkit"
    exit 1
fi
seqkit grep -n -r -f "${VOTU_OUT}/.rep_ids_regex.txt" "${CONFIRMED_FA}" \
    > "${VOTU_OUT}/votu_representatives.fa"

n_reps=$(grep -c '^>' "${VOTU_OUT}/votu_representatives.fa" 2>/dev/null || echo 0)
log "  ${n_reps} representative sequences extracted"
if [[ "${n_reps}" -lt "${n_votus}" ]]; then
    log "  WARNING: extracted ${n_reps} but expected ${n_votus} -- check for CheckV suffix"
    log "  mismatches (same class of issue documented in 03_merge_checkv.sh)"
fi

###############################################################################
# 6. Build vOTU index: votu_id | representative | n_members | member MAG indices
#    (member_mag_indices is what makes this queryable for the cross-host
#    core/shared/host-specific breakdown downstream -- join against
#    mag_index.tsv's host-species metadata by the M###### prefix.)
###############################################################################
log "Building vOTU index with member MAG provenance..."
awk -F'\t' '
    BEGIN { OFS="\t"; print "votu_id","representative","n_members","member_ids" }
    {
        rep = $1
        n = split($2, members, ",")
        votu_id = sprintf("vOTU_%06d", NR)
        print votu_id, rep, n, $2
    }
' "${CLUSTERS_OUT}" > "${VOTU_OUT}/vOTU_index.tsv"

rm -f "${VOTU_OUT}/.rep_ids.txt" "${VOTU_OUT}/.rep_ids_regex.txt"

log "Step 6 finished."
log "  ${n_votus} vOTUs -> ${VOTU_OUT}/votu_representatives.fa"
log "  cluster membership -> ${VOTU_OUT}/vOTU_index.tsv"
log "  source MAG / host species: join member_ids (M###### prefix) against ${MAG_INDEX_TSV}"
log "Next: 07 (PhaBOX2 taxonomy on vOTU representatives), 09 (PropagAtE), 10b (CRISPR matching)"
