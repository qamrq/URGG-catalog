#!/usr/bin/env bash
# ============================================================
# 53_dbcan_restore_gene_ids.sh
#
# Restores full catalogue gene IDs in dbCAN's overview.tsv.
#
# FINDING: dbCAN annotated the correct catalogue
# (nr_gene_catalog_rep_seq.fasta), but gene IDs longer than 80 characters
# were cut to 80 characters and given a "_<n>" counter to keep them
# unique (e.g. "...bin.135~TP20200812" + "_1015"). IDs of <=80 characters
# were left unchanged. This truncation — not stray BASALT bins — is the
# source of the tilde-less IDs and the "unresolved genome" names found
# earlier (34c/34e): long genome names were cut before the '~'.
#
# METHOD: dbCAN's input copy (uniInput.faa) lists the same proteins in the
# same order as the catalogue FASTA, only with shortened IDs. Records are
# paired by position, and every pair is VERIFIED by sequence (length +
# first/last 30 residues, stop codons removed) and by the truncation rule
# (short ID unchanged; or long ID cut to exactly 80 characters; or cut to
# 80 characters plus a "_<n>" counter). v2: the first run wrongly rejected
# 5,685 IDs cut to 80 characters WITHOUT a counter whose last characters
# happened to be "_<digits>" (e.g. "..._length_86924"), because the counter
# was stripped before comparing.
# The script stops if either check fails.
#
# OUTPUTS (OUT, a new folder — nothing existing is modified):
#   id_map.tsv            dbCAN_id <tab> catalogue_rep_id   (all 37.1M)
#   overview.tsv          dbCAN overview with 'Gene ID' restored to catalogue IDs
#   53_restore_report.txt
# ============================================================
set -euo pipefail

A=${A:-/disk252/qamar/metagenome_analysis/annotation}
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

# ---------------- 1. Fingerprint every record in both FASTAs ----------------
# index <tab> id <tab> length:first30:last30   (stop codons '*' removed)
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

# ---------------- 3. Rewrite overview.tsv with catalogue IDs ----------------
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