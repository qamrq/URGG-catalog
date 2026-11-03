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
#  SCRIPT  : 04a_diamond_nr.sh  -  DIAMOND search of genes without eggNOG annotation against NCBI nr
#  AUTHOR  : QRQ
#  INPUT   : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#            annotation/eggnog_output/nr_catalog.emapper.annotations
#            $URGG_DB/ncbi_nr/nr.dmnd
#  OUTPUT  : annotation/diamond_nr_output/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/04a_diamond_nr_$(date +%Y%m%d_%H%M%S).log") 2>&1

activate_env "${ENV_DIAMOND:-diamond}"

# ---- Paths (edit PROT_FASTA if auto-detection below fails) ----
NR_DB=${URGG_DB}/ncbi_nr/nr.dmnd
EGGNOG_ANN=${ANNOT}/eggnog_output/nr_catalog.emapper.annotations
CATALOG_DIR=${ANNOT}/gene_catalog
OUT=${ANNOT}/diamond_nr_output
WORK=${OUT}/work
mkdir -p "${OUT}" "${WORK}"

# ---- Tunables ----
MAX_TARGETS=25
EVALUE=1e-5
BLOCK_SIZE=8
INDEX_CHUNKS=1
DMND_TMP=/tmp/diamond_tmp
mkdir -p "${DMND_TMP}"

log "Validating DIAMOND nr database..."
if [[ ! -s "${NR_DB}" ]]; then
    log "ERROR: ${NR_DB} not found. Has 'diamond makedb' finished?"
    exit 1
fi
diamond dbinfo -d "${NR_DB}" || { log "ERROR: dbinfo failed — DB is incomplete or corrupt."; exit 1; }

log "Locating protein FASTA for the non-redundant gene catalog..."
PROT_FASTA="${URGG_PROJECT}/annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta"
for cand in \
    "${CATALOG_DIR}/nr_gene_catalog_rep_seq.faa" \
    "${CATALOG_DIR}/nr_gene_catalog_rep_seq.fasta" \
    "${CATALOG_DIR}"/*.faa ; do
    [[ -s "${cand}" ]] && { PROT_FASTA="${cand}"; break; }
done

if [[ -z "${PROT_FASTA}" ]]; then
    log "ERROR: no candidate protein FASTA found in ${CATALOG_DIR}"
    log "Contents:"; find "${CATALOG_DIR}" -maxdepth 1 -type f -printf '  %f\n'
    exit 1
fi
log "Candidate protein FASTA: ${PROT_FASTA}"

FIRST_SEQ=$(awk '/^>/{if(n++)exit;next}{printf "%s",$0}' "${PROT_FASTA}" | head -c 400)
NON_ACGTN=$(printf '%s' "${FIRST_SEQ}" | tr -d 'ACGTNacgtn' | wc -c)
if [[ "${NON_ACGTN}" -lt 20 ]]; then
    log "ERROR: ${PROT_FASTA} looks like NUCLEOTIDE sequence, not protein."
    log "       blastp requires amino acids. Use the Prodigal .faa output"
    log "       (protein translations), or switch to 'diamond blastx'"
    log "       (far slower). Aborting rather than producing garbage."
    exit 1
fi
log "Confirmed protein (amino acid) sequences."

N_TOTAL=$(grep -c '^>' "${PROT_FASTA}")
log "Total proteins in catalog: ${N_TOTAL}"

ANNOTATED_IDS=${WORK}/eggnog_annotated_ids.txt
if [[ -s "${ANNOTATED_IDS}" ]]; then
    log "Reusing existing annotated-ID list: ${ANNOTATED_IDS}"
else
    log "Extracting eggNOG-annotated gene IDs..."
    grep -v '^#' "${EGGNOG_ANN}" | cut -f1 | sort -u > "${ANNOTATED_IDS}"
fi
N_ANNOT=$(wc -l < "${ANNOTATED_IDS}")
log "eggNOG-annotated genes: ${N_ANNOT}"

UNANNOT_FASTA=${WORK}/unannotated_proteins.faa
if [[ -s "${UNANNOT_FASTA}" ]]; then
    log "Reusing existing unannotated subset: ${UNANNOT_FASTA}"
else
    log "Extracting proteins with NO eggNOG annotation..."
    awk 'NR==FNR { seen[$0]=1; next }
         /^>/ { id=substr($1,2); keep=!(id in seen) }
         keep { print }' \
        "${ANNOTATED_IDS}" "${PROT_FASTA}" > "${UNANNOT_FASTA}"
fi
N_UNANNOT=$(grep -c '^>' "${UNANNOT_FASTA}")
log "Unannotated proteins to search against nr: ${N_UNANNOT}"

# Consistency check — this should be close to N_TOTAL - N_ANNOT.
EXPECTED=$(( N_TOTAL - N_ANNOT ))
log "Expected ~${EXPECTED} from arithmetic; got ${N_UNANNOT}"
if [[ "${N_UNANNOT}" -eq 0 ]]; then
    log "ERROR: subset is empty. Likely an ID-format mismatch between the"
    log "       eggNOG query column and the FASTA headers. Compare:"
    log "       head -1 ${ANNOTATED_IDS}"
    log "       grep -m1 '^>' ${PROT_FASTA}"
    exit 1
fi

HITS=${OUT}/nr_hits.tsv
log "Starting DIAMOND blastp: ${N_UNANNOT} queries vs nr"
log "  max-target-seqs=${MAX_TARGETS}  evalue=${EVALUE}  block-size=${BLOCK_SIZE}"

diamond blastp \
    --query "${UNANNOT_FASTA}" \
    --db "${NR_DB}" \
    --out "${HITS}" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle \
    --evalue "${EVALUE}" \
    --max-target-seqs "${MAX_TARGETS}" \
    --threads "${THREADS}" \
    --block-size "${BLOCK_SIZE}" \
    --index-chunks "${INDEX_CHUNKS}" \
    --tmpdir "${DMND_TMP}"

log "DIAMOND blastp complete."

N_HIT_GENES=$(cut -f1 "${HITS}" | sort -u | wc -l)
N_NO_HIT=$(( N_UNANNOT - N_HIT_GENES ))
log "-------------------------------------------------------"
log "Catalog proteins:                 ${N_TOTAL}"
log "  with eggNOG annotation:         ${N_ANNOT}"
log "  without eggNOG annotation:      ${N_UNANNOT}"
log "    of those, with an nr hit:     ${N_HIT_GENES}"
log "    of those, NO hit anywhere:    ${N_NO_HIT}   <-- candidate truly novel genes"
log "-------------------------------------------------------"
log "Hits table: ${HITS}"

# The no-hit set is the headline novelty result — write it out.
NOHIT_IDS=${OUT}/no_nr_hit_gene_ids.txt
grep '^>' "${UNANNOT_FASTA}" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "${WORK}/unannot_ids.txt"
cut -f1 "${HITS}" | sort -u > "${WORK}/nr_hit_ids.txt"
comm -23 "${WORK}/unannot_ids.txt" "${WORK}/nr_hit_ids.txt" > "${NOHIT_IDS}"
log "Gene IDs with no homolog in eggNOG or nr: ${NOHIT_IDS}"
