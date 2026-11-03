#!/bin/bash
# 33_rvd_virome_comparison.sh
# Compare this project's 8,980 species-level vOTUs against the published
# Rumen Virome Database (RVD; Yan & Yu 2024, Microbiome 12:82).
# Uses CheckV's anicalc/aniclust (the field-standard tool for viral vOTU
# comparison), NOT skani, at 95% ANI / 85% coverage-of-shorter-contig,
# matching RVD's own published construction method.
#
# One-directional comparison: for each of our 8,980 vOTUs, is there an
# RVD vOTU that represents "the same virus" at 95%/85%? We do NOT run
# a full redundancy-merge clustering (that's what aniclust is designed
# for when combining two sets into one non-redundant set) — we just want
# matched vs novel classification, so we read anicalc's pairwise ANI
# table directly. If you also want RVD-style representative selection
# downstream, aniclust can be run on the combined table (see bottom).

set -euo pipefail


# ---- paths -------------------------------------------------------------
PROJ_DIR="/disk252/qamar/metagenome_analysis/phase_5_virome"
WORK_DIR="${PROJ_DIR}/09_rvd_comparison"
QUERY_FA="${PROJ_DIR}/07_votu_clustering/votu_representatives.fa"
RVD_FA="/disk252/qamar/databases/rvd/RVD.fa"
RVD_DB="${WORK_DIR}/rvd_blastdb/rvd_db"

THREADS=64

mkdir -p "${WORK_DIR}/rvd_blastdb"
cd "${WORK_DIR}"

# ---- 1. Build BLAST db on RVD (build once, reuse across runs) ----------
if [ ! -f "${RVD_DB}.nin" ] && [ ! -f "${RVD_DB}.00.nin" ]; then
    echo "[1/4] Building BLAST database from RVD.fa ..."
    makeblastdb -in "${RVD_FA}" -dbtype nucl -out "${RVD_DB}"
else
    echo "[1/4] RVD BLAST database already exists, skipping build."
fi

# ---- 2. All-vs-RVD BLASTN (query is small: 8,980 seqs; no need to ------
#         chunk the query -- chunking matters when the QUERY is huge,
#         not here) ------------------------------------------------------
BLAST_OUT="${WORK_DIR}/votu_vs_rvd_blast.tsv"

if [ ! -s "${BLAST_OUT}" ]; then
    echo "[2/4] Running BLASTN (query vOTUs vs RVD) ..."
    # outfmt '6 std qlen slen' is REQUIRED by anicalc (needs query/target
    # lengths appended after the 12 standard columns).
    # -max_target_seqs kept modest (25) since we only need each query's
    # best RVD hit(s), not an exhaustive all-pairs graph (unlike aniclust's
    # own redundancy-clustering use case, which wants -max_target_seqs 10000).
    blastn -query "${QUERY_FA}" \
           -db "${RVD_DB}" \
           -outfmt '6 std qlen slen' \
           -max_target_seqs 25 \
           -perc_identity 90 \
           -num_threads "${THREADS}" \
           -out "${BLAST_OUT}"
else
    echo "[2/4] BLAST output already exists, skipping."
fi

# ---- 3. anicalc: pairwise ANI + alignment fraction ----------------------
ANI_OUT="${WORK_DIR}/votu_vs_rvd_ani.tsv"

echo "[3/4] Running anicalc ..."
anicalc -i "${BLAST_OUT}" -o "${ANI_OUT}"

# ---- 4. Classify matched vs novel at 95% ANI / 85% coverage-of-shorter -
echo "[4/4] Classifying matched vs novel vOTUs ..."
python3 "$(dirname "$0")/classify_rvd_matches.py" \
    --query-fa "${QUERY_FA}" \
    --ani "${ANI_OUT}" \
    --min-ani 95 \
    --min-cov 85 \
    --out-table "${WORK_DIR}/votu_rvd_classification.tsv" \
    --out-summary "${WORK_DIR}/votu_rvd_summary.txt"

echo "Done. See:"
echo "  ${WORK_DIR}/votu_rvd_classification.tsv  (per-vOTU matched/novel calls)"
echo "  ${WORK_DIR}/votu_rvd_summary.txt          (matched/novel % summary)"

# ---------------------------------------------------------------------
# OPTIONAL: if you later want RVD-style redundancy clustering (combined
# representative selection across query+RVD, not just matched/novel
# classification), run aniclust on a combined FASTA + combined ani.tsv:
#
#   cat "${QUERY_FA}" "${RVD_FA}" > combined.fna
#   aniclust --fna combined.fna --ani "${ANI_OUT}" \
#            --out combined_clusters.tsv \
#            --min_ani 95 --min_tcov 85 --min_qcov 0
#
# Not run by default -- the classification step above already answers
# "similarity and difference vs RVD" without needing a full merge.
# ---------------------------------------------------------------------
