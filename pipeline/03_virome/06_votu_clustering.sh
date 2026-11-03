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
#  SCRIPT  : 06_votu_clustering.sh  -  vOTU clustering (95% ANI, 85% AF)
#  AUTHOR  : QRQ
#  INPUT   : confirmed_viral.fa
#  OUTPUT  : phase_5_virome/07_votu_clustering/votu_representatives.fa, vOTU_index.tsv
# ==========================================================

source "$(dirname "$0")/config.sh"
activate_env "${ENV_VIRSORTER2}"

log "Step 6: vOTU clustering (95% ANI / 85% AF) on confirmed viral set"

CONFIRMED_FA="${CONFIRMED_VIRAL_OUT}/confirmed_viral.fa"
if [[ ! -s "${CONFIRMED_FA}" ]]; then
    log "ERROR: ${CONFIRMED_FA} not found or empty. Run 03_merge_checkv.sh first."
    exit 1
fi
n_input=$(grep -c '^>' "${CONFIRMED_FA}")
log "Input: ${n_input} confirmed viral sequences"

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

log "Computing pairwise ANI/AF (anicalc.py)..."
ANI_OUT="${VOTU_OUT}/ani.tsv"
python "${ANICALC}" -i "${BLAST_OUT}" -o "${ANI_OUT}" \
    > "${LOG_DIR}/votu_anicalc.log" 2>&1

if [[ ! -s "${ANI_OUT}" ]]; then
    log "ERROR: anicalc.py produced no output -- see ${LOG_DIR}/votu_anicalc.log"
    exit 1
fi

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
