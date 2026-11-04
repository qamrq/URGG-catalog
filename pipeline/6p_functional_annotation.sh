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
#  SCRIPT  : 6p_functional_annotation.sh  -  eggNOG, NCBI nr (unannotated genes), VFDB and dbCAN annotation of the gene catalog
#  AUTHOR  : QRQ
#  INPUT   : annotation/gene_catalog/nr_gene_catalog_rep_seq.fasta
#  OUTPUT  : annotation/eggnog_output/, diamond_nr_output/, diamond_vfdb/, dbcan_output_restored/
#  USAGE   : bash 6p_functional_annotation.sh [eggnog nr vfdb dbcan dbcan_ids]   (no argument: all steps in this order)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

step_eggnog() (
set +u
set -u
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/04c_eggnog_$(date +%Y%m%d_%H%M%S).log") 2>&1

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/eggnog_output
mkdir -p "${OUT}"

log "Running eggNOG-mapper (DIAMOND mode)..."
emapper.py \
    -i "${QUERY}" \
    --itype proteins \
    -m diamond \
    --data_dir "${EGGNOG_DB}" \
    -o nr_catalog \
    --output_dir "${OUT}" \
    --cpu "${THREADS}" \
    --override

log "eggNOG-mapper done."
log "Annotations: ${OUT}/nr_catalog.emapper.annotations"
log "  -> KO column = KEGG orthologs; KEGG_Pathway / KEGG_Module = pathways; CAZy column present too."
)

step_nr() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/04a_diamond_nr_$(date +%Y%m%d_%H%M%S).log") 2>&1

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
DMND_TMP="${DMND_TMP:-${TMPDIR:-/tmp}/diamond_tmp}"
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
)

step_vfdb() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/04b_diamond_vfdb_$(date +%Y%m%d_%H%M%S).log") 2>&1

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/diamond_vfdb
mkdir -p "${OUT}"
RAW=${OUT}/diamond_vfdb_raw.tsv
FILT=${OUT}/diamond_vfdb_filtered.tsv

log "DIAMOND blastp vs VFDB..."
diamond blastp \
    -q "${QUERY}" \
    -d "${VFDB_DMND}" \
    -o "${RAW}" \
    --threads "${THREADS}" \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --id 80 --query-cover 70 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle

awk -F'\t' '$3>=80 && $5<=1e-5' "${RAW}" > "${FILT}"

log "VFDB hits (raw): $(wc -l < "${RAW}")"
log "VFDB hits (>=80% id): $(wc -l < "${FILT}")"
log "Filtered virulence hits -> ${FILT}"
)

step_dbcan() (
set -euo pipefail
source "$(dirname "$0")/config_annotation.sh"
exec > >(tee -a "${LOGS}/04d_dbcan_$(date +%Y%m%d_%H%M%S).log") 2>&1

QUERY=${CATALOG}/nr_gene_catalog_rep_seq.fasta
OUT=${ANNOT}/dbcan_output
mkdir -p "${OUT}"

log "Running run_dbcan CAZyme_annotation (dbCAN v5.x subcommand syntax)..."
run_dbcan CAZyme_annotation \
    --input_raw_data "${QUERY}" \
    --mode protein \
    --output_dir "${OUT}" \
    --db_dir "${DBCAN_DB}" \
    --threads "${THREADS}"

log "dbCAN done. Consensus CAZyme table: ${OUT}/overview.txt"
)

step_dbcan_ids() (
set -euo pipefail

A=${A:-${URGG_PROJECT}/annotation}
REP=${REP:-${A}/gene_catalog/nr_gene_catalog_rep_seq.fasta}
UNI=${UNI:-${A}/dbcan_output/uniInput.faa}
OV=${OV:-${A}/dbcan_output/overview.tsv}
OUT=${OUT:-${A}/dbcan_output_restored}
TRUNC=${TRUNC:-80}

mkdir -p "${OUT}"
exec > >(tee -a "${OUT}/53_restore_report.txt") 2>&1
log() { echo "[$(date '+%F %T')] $*"; }
for f in "${REP}" "${UNI}" "${OV}"; do [[ -s "$f" ]] || { log "ERROR: missing $f"; exit 1; }; done
[[ -s "${OUT}/overview.tsv" ]] && { log "ERROR: ${OUT}/overview.tsv exists — move it aside first"; exit 1; }

fp() {
  awk 'BEGIN{OFS="\t"}
       function emit(){ gsub(/\*/,"",s); L=length(s);
                        print n, id, L ":" substr(s,1,30) ":" substr(s,(L>30?L-29:1)) }
       /^>/ { if (id!="") emit(); n++; id=substr($1,2); s=""; next }
       { s = s $0 }
       END { if (id!="") emit() }' "$1"
}
log "Fingerprinting catalogue FASTA (${REP}) ..."
fp "${REP}" > "${OUT}/rep.fp.tsv"
log "Fingerprinting dbCAN input copy (${UNI}) ..."
fp "${UNI}" > "${OUT}/uni.fp.tsv"

n_rep=$(wc -l < "${OUT}/rep.fp.tsv"); n_uni=$(wc -l < "${OUT}/uni.fp.tsv")
log "Records: catalogue ${n_rep} | dbCAN input ${n_uni}"
[[ "${n_rep}" -eq "${n_uni}" ]] || { log "ERROR: record counts differ — positional pairing is not valid. Stop."; exit 1; }

# ---------------- 2. Pair by position and verify ----------------
log "Pairing records by position and verifying sequence + truncation rule ..."
paste "${OUT}/rep.fp.tsv" "${OUT}/uni.fp.tsv" | awk -v T="${TRUNC}" -v MAP="${OUT}/id_map.tsv" '
  BEGIN { FS=OFS="\t" }
  {
    rep_id=$2; rep_fp=$3; uni_id=$5; uni_fp=$6
    if (rep_fp != uni_fp) { seq_bad++; if (seq_bad<=5) print "SEQ MISMATCH", NR, rep_id, uni_id > "/dev/stderr" }
    if (uni_id == rep_id) { same++ }
    else if (length(rep_id) > T && uni_id == substr(rep_id,1,T)) { trunc_plain++ }   # cut to T chars, no counter
    else {
      pre = uni_id; sub(/_[0-9]+$/, "", pre)                                         # cut to T chars + "_<n>" counter
      if (length(rep_id) > T && length(pre) == T && substr(rep_id,1,T) == pre) trunc_ctr++
      else { rule_bad++; if (rule_bad<=5) print "RULE MISMATCH", NR, rep_id, uni_id > "/dev/stderr" }
    }
    print uni_id, rep_id > MAP
  }
  END {
    printf "  identical IDs (<=%d chars): %d\n  truncated to %d chars, no counter: %d\n  truncated to %d chars + _<n> counter: %d\n  sequence mismatches: %d\n  truncation-rule violations: %d\n",
           T, same, T, trunc_plain+0, T, trunc_ctr+0, seq_bad+0, rule_bad+0
    if (seq_bad+rule_bad > 0) exit 1
  }' || { log "ERROR: verification failed — see mismatches above. Map deleted."; rm -f "${OUT}/id_map.tsv"; exit 1; }

n_dup=$(cut -f1 "${OUT}/id_map.tsv" | sort | uniq -d | wc -l)
log "Duplicate dbCAN IDs in map: ${n_dup} (must be 0)"
[[ "${n_dup}" -eq 0 ]] || { log "ERROR: dbCAN IDs are not unique — mapping is ambiguous."; exit 1; }
rm -f "${OUT}/rep.fp.tsv" "${OUT}/uni.fp.tsv"

log "Restoring Gene IDs in overview.tsv ..."
awk -v MAPF="${OUT}/id_map.tsv" '
  BEGIN { FS=OFS="\t" }
  NR==FNR { if (FNR>1) need[$1]=1; next }                       # pass 1: overview IDs
  FILENAME==MAPF { if ($1 in need) m[$1]=$2; next }              # pass 2: map (only needed rows kept)
  FNR==1 { print; next }                                         # pass 3: header
  { if ($1 in m) { $1 = m[$1]; ok++ } else { miss++; if (miss<=5) print "UNMAPPED: " $1 > "/dev/stderr" }; print }
  END { printf "  overview rows restored: %d | unmapped: %d\n", ok, miss+0 > "/dev/stderr"; if (miss>0) exit 1 }
' "${OV}" "${OUT}/id_map.tsv" "${OV}" > "${OUT}/overview.tsv" \
  || { log "ERROR: some overview IDs could not be mapped."; exit 1; }

# ---------------- 4. Summary ----------------
log "Genome coverage before vs after restoring IDs:"
Rscript - "${OV}" "${OUT}/overview.tsv" <<'EOF'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
old <- fread(a[1], select = 1, quote = "")[[1]]
new <- fread(a[2], select = 1, quote = "")[[1]]
g_old <- ifelse(grepl("~", old, fixed = TRUE), sub("~.*$", "", old), sub("_[0-9]+$", "", old))   # old 34_ logic
g_new <- sub("~.*$", "", new)
cat(sprintf("  genes: %s | with '~' after restore: %.4f%%\n", format(length(new), big.mark=","), 100*mean(grepl("~", new, fixed=TRUE))))
cat(sprintf("  distinct genomes, old parsing: %s | restored: %s\n",
            format(uniqueN(g_old), big.mark=","), format(uniqueN(g_new), big.mark=",")))
cat(sprintf("  genome-name length > %d chars among restored genomes: %s\n", 80,
            format(sum(nchar(unique(g_new)) > 80), big.mark=",")))
EOF
log "Done. Restored overview: ${OUT}/overview.tsv"
)

for s in ${@:-eggnog nr vfdb dbcan dbcan_ids}; do
  echo "==> 6p_functional_annotation.sh: $s"
  "step_$s"
done
