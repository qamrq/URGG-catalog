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
#  SCRIPT  : 30a_target_genome_tree.sh  -  ML tree of the RUC lineage, GTDB relatives and the unclassified sheep genome (Fig. 5a)
#  AUTHOR  : QRQ
#  INPUT   : annotation/gtdbtk_output_v232/
#            $URGG_DB/gtdbtk_database/release232
#  OUTPUT  : annotation/results_phase3/s32_taxonomy_novelty/target_ml_tree/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

set -euo pipefail

ANNOT=${ANNOT:-${URGG_PROJECT}/annotation}
GTDBTK_OUT=${GTDBTK_OUT:-${ANNOT}/gtdbtk_output_v232}
SUMMARY=${SUMMARY:-${GTDBTK_OUT}/gtdbtk.all.summary.tsv}
GTDB_DIR=${GTDB_DIR:-${URGG_DB}/gtdbtk_database/release232}
LINEAGE_BUILD_DIR=${LINEAGE_BUILD_DIR:-${ANNOT}/iqtree_output/tree_build_lineage_only_fresh}
LABEL_CSV=${LABEL_CSV:-${URGG_PROJECT}/phase_3_analysis/scripts/results/lineage_full_classification_table.csv}
WORK=${WORK:-${ANNOT}/results_phase3/s32_taxonomy_novelty/target_ml_tree}
LOGS=${LOGS:-${ANNOT}/logs}

TARGET=${TARGET:-CNP0001390__Sheep__Feces__T1900044680_TS20190146_bin.8}
LINEAGE_REGEX=${LINEAGE_REGEX:-'UMGS1840|HGM11507'}
EXPECT_LINEAGE=${EXPECT_LINEAGE:-76}

REF_IDS=${REF_IDS:-"RS_GCF_018918095.1 GB_GCA_017414505.1 RS_GCF_002252565.1 GB_GCA_030156665.1 GB_GCA_022360955.1 GB_GCA_012735465.1 GB_GCA_045655245.1 GB_GCA_017620795.1"}
REF_MSA=${REF_MSA:-}
GTDB_TAXONOMY=${GTDB_TAXONOMY:-${GTDB_DIR}/taxonomy/gtdb_taxonomy.tsv}
NEIGH_REGEX=${NEIGH_REGEX:-'o__SK-Y3;|o__RGIG8785;'}
N_NEIGH_CAT=${N_NEIGH_CAT:-6}
OUT_REGEX=${OUT_REGEX:-'c__Bacilli;'}
N_OUT=${N_OUT:-4}
MIN_MSA_PERCENT_NEIGH=${MIN_MSA_PERCENT_NEIGH:-70}

MSA_FILES=${MSA_FILES:-}
MIN_SEQ_COVER=${MIN_SEQ_COVER:-0.25}
MIN_TRIM_COLS=${MIN_TRIM_COLS:-1000}
THREADS=${THREADS:-32}
SEED=${SEED:-42}
UFBOOT=${UFBOOT:-1000}
STOP_AFTER=${STOP_AFTER:-}
FORCE=${FORCE:-0}
ALLOW_MISSING=${ALLOW_MISSING:-0}

mkdir -p "${WORK}" "${LOGS}"
exec > >(tee -a "${LOGS}/81_target_lineage_ml_tree_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] CHECKPOINT FAILED: $*" >&2; exit 1; }
stop_if() { [[ "${STOP_AFTER}" == "$1" ]] && { log "STOP_AFTER=$1: stopping here as requested."; exit 0; }; return 0; }
[[ "${FORCE}" == "1" ]] && { log "FORCE=1: clearing earlier outputs in ${WORK}"; rm -f "${WORK}"/*.tsv "${WORK}"/*.fasta "${WORK}"/*.txt "${WORK}"/*.nwk "${WORK}"/ml_tree.* "${WORK}"/.done_* ; }

log "=== C0: tools and inputs ==="
[[ -s "${SUMMARY}" ]] || die "GTDB-Tk summary not found: ${SUMMARY}"
command -v python3 >/dev/null || die "python3 not found"
command -v Rscript >/dev/null || die "Rscript not found (needed for the final tree checks)"
Rscript -e 'suppressMessages(library(ape))' 2>/dev/null || die "R package 'ape' not installed"
IQTREE=${IQTREE:-$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)}
TRIMAL=${TRIMAL:-$(command -v trimal || true)}
if [[ -n "${MSA_FILES}" ]]; then
  [[ -n "${IQTREE}" ]] || die "IQ-TREE not on PATH (looked for iqtree3, iqtree2, iqtree)."
  [[ -n "${TRIMAL}" ]] || die "trimal not on PATH."
  log "IQ-TREE: ${IQTREE} ($(${IQTREE} --version 2>&1 | head -n1))"
  log "trimAl : ${TRIMAL} ($(${TRIMAL} --version 2>&1 | head -n1))   <- fill into Methods"
fi
log "C0 passed."

IDS_TSV="${WORK}/wanted_ids.tsv"
if [[ ! -s "${IDS_TSV}" ]]; then
log "=== C1: selecting genomes ==="
python3 - "${SUMMARY}" "${TARGET}" "${LINEAGE_REGEX}" "${EXPECT_LINEAGE}" "${NEIGH_REGEX}" "${N_NEIGH_CAT}" \
          "${OUT_REGEX}" "${N_OUT}" "${MIN_MSA_PERCENT_NEIGH}" "${REF_IDS}" "${IDS_TSV}" "${GTDB_TAXONOMY}" <<'PY'
import sys, csv, re
(summary, target, lin_rx, expect, neigh_rx, n_neigh, out_rx, n_out, min_msa, ref_ids, out_tsv, tax_file) = sys.argv[1:13]
expect, n_neigh, n_out, min_msa = int(expect), int(n_neigh), int(n_out), float(min_msa)

rows = list(csv.DictReader(open(summary), delimiter="\t"))
by_id = {r["user_genome"]: r for r in rows}
def fail(m): sys.exit("CHECKPOINT FAILED: " + m)

if target not in by_id: fail(f"target {target} is not in the GTDB-Tk summary")
tcls = by_id[target]["classification"]
print(f"target classification: {tcls}")
if not re.search(r"c__Clostridia;o__;", tcls):
    fail("target is not 'c__Clostridia;o__;...' as expected (order empty). Wrong TARGET, or the summary changed.")
if re.search(lin_rx, tcls): fail("target unexpectedly matches the lineage regex")

lineage = [r["user_genome"] for r in rows if re.search(lin_rx, r["classification"])]
print(f"lineage genomes ({lin_rx}): {len(lineage)}")
if len(lineage) != expect: fail(f"expected {expect} lineage genomes, found {len(lineage)}")

def msa_pct(r):
    try: return float(r.get("msa_percent", ""))
    except ValueError: return None          # blank/NA = not in the alignment (e.g. ANI pre-screen)

def pick(rx, n, min_pct):
    pool = [r for r in rows if re.search(rx, r["classification"]) and r["user_genome"] not in set(lineage)
            and r["user_genome"] != target and msa_pct(r) is not None and msa_pct(r) >= min_pct]
    groups = {}
    for r in sorted(pool, key=lambda r: -msa_pct(r)):
        groups.setdefault(r["classification"].split(";")[2], []).append(r)   # round-robin across orders
    out, i = [], 0
    while len(out) < n and any(groups.values()):
        for g in list(groups):
            if groups[g] and len(out) < n: out.append(groups[g].pop(0))
        i += 1
    return out

import os
tax = {}
if os.path.exists(tax_file):
    for l in open(tax_file):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2: tax[a[0]] = a[1]
else:
    print(f"NOTE: {tax_file} not found; references will not be grouped by order in the report.")
refs = list(dict.fromkeys(ref_ids.split()))          # de-duplicated, order kept
refs = [g for g in refs if g not in by_id]
print(f"GTDB reference genomes requested: {len(refs)} ({sum(1 for g in refs if g in tax or re.sub(r'^(RS_|GB_)', '', g) in tax)} with a GTDB taxonomy entry)")
neigh = pick(neigh_rx, n_neigh, min_msa)
outg  = pick(out_rx, n_out, 90)
print(f"catalogue neighbours chosen ({neigh_rx}): {len(neigh)}  | outgroup chosen ({out_rx}): {len(outg)}")
if len(outg) < 2: fail(f"fewer than 2 outgroup genomes found for {out_rx!r}; set OUT_REGEX to another well-covered class")
if not neigh: print("NOTE: no catalogue genomes of the neighbour orders; neighbours will come from the GTDB reference genomes only.")

with open(out_tsv, "w") as fh:
    fh.write("id\trole\tclassification\n")
    fh.write(f"{target}\ttarget\t{tcls}\n")
    for g in lineage: fh.write(f"{g}\tlineage\t{by_id[g]['classification']}\n")
    for g in refs: fh.write(f"{g}\tref_neighbour\t{tax.get(g, tax.get(re.sub(r'^(RS_|GB_)', '', g), 'GTDB reference (taxonomy not found)'))}\n")
    for r in neigh: fh.write(f"{r['user_genome']}\tcat_neighbour\t{r['classification']}\n")
    for r in outg:  fh.write(f"{r['user_genome']}\toutgroup\t{r['classification']}\n")
print("\nOutgroup genomes:"); [print("  ", r["user_genome"], r["classification"]) for r in outg]
print("Catalogue neighbours:"); [print("  ", r["user_genome"], r["classification"]) for r in neigh]
PY
fi
log "role counts:"; cut -f2 "${IDS_TSV}" | tail -n +2 | sort | uniq -c
log "C1 passed. Table: ${IDS_TSV}"
stop_if select

if [[ -z "${MSA_FILES}" ]]; then
  log "=== C2: MSA_FILES not set -> SCAN MODE (nothing is built) ==="
  log "The IQ-TREE run that made the 76-tip tree recorded its input alignment:"
  grep -h -E "Input file name|Command:" "${LINEAGE_BUILD_DIR}"/*.iqtree 2>/dev/null | head -n 5 || log "  (no .iqtree file found in ${LINEAGE_BUILD_DIR})"
  CAND="${WORK}/candidate_alignments.txt"
  { find "${GTDBTK_OUT}" "${LINEAGE_BUILD_DIR}" -maxdepth 4 -type f \( -iname "*msa*" -o -iname "*align*" -o -iname "*aln*" -o -iname "*trim*" \) 2>/dev/null
    find "${LINEAGE_BUILD_DIR}" -maxdepth 2 -type f \( -iname "*.fa" -o -iname "*.fa.gz" -o -iname "*.fasta" -o -iname "*.fasta.gz" -o -iname "*.faa" -o -iname "*.faa.gz" \) 2>/dev/null
    grep -h "Input file name" "${LINEAGE_BUILD_DIR}"/*.iqtree 2>/dev/null | sed 's/.*: //' | while read -r f; do [[ -f "$f" ]] && echo "$f"; done
    [[ -d "${GTDB_DIR}/msa" ]] && find "${GTDB_DIR}/msa" -maxdepth 2 -type f 2>/dev/null
  } | grep -v -E "\.(log|treefile|tree|tsv|csv|txt|ckp\.gz|bionj|mldist|iqtree|contree|splits\.nex|uniqueseq\.phy)$" | sort -u > "${CAND}" || true
  log "candidate files: $(wc -l < "${CAND}")  (listed in ${CAND})"
  python3 - "${IDS_TSV}" "${CAND}" <<'PY'
import sys, gzip, os, re
ids_tsv, cand = sys.argv[1:3]
want, role = set(), {}
for l in list(open(ids_tsv))[1:]:
    i, r, _ = l.rstrip("\n").split("\t"); want.add(re.sub(r"^(RS_|GB_)", "", i)); role[re.sub(r"^(RS_|GB_)", "", i)] = r
def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
print(f"\n{'file':<90} {'MB':>7} {'len':>6} {'target':>6} {'lineage':>8} {'refs':>5} {'cat_nb+out':>10}")
for p in [l.strip() for l in open(cand) if l.strip()]:
    try:
        n = {"target": 0, "lineage": 0, "ref_neighbour": 0, "cat_neighbour": 0, "outgroup": 0}; length = None; first = True; nseq = 0
        with opener(p) as fh:
            seq_len = 0; head = None
            for line in fh:
                if line.startswith(">"):
                    if head is not None and length is None: length = seq_len
                    head = re.sub(r"^(RS_|GB_)", "", line[1:].split()[0]); seq_len = 0; nseq += 1
                    if head in want: n[role[head]] += 1
                else: seq_len += len(line.strip())
            if length is None: length = seq_len
        if nseq == 0: print(f"{p:<90} (not FASTA / empty)"); continue
        print(f"{p[-90:]:<90} {os.path.getsize(p)/1e6:>7.0f} {length:>6} {n['target']:>6} {n['lineage']:>8} {n['ref_neighbour']:>5} {n['cat_neighbour']+n['outgroup']:>10}")
    except Exception as e:
        print(f"{p[-90:]:<90} (unreadable: {e})")
PY
  cat <<EOF

Choose the files and run again, for example:
  * a file with the lineage genomes AND the target (look for 76 / 1);
  * a file with the reference genomes (refs column > 0): set it as REF_MSA;
  * outgroup / catalogue neighbours are usually in gtdbtk.bac120.user_msa.fasta.gz.
  All chosen files must have the SAME length (len column), normally 5040 for
  GTDB-Tk's bac120 alignment. A file with a shorter length is already trimmed
  and cannot be combined with untrimmed ones.
  MSA_FILES="fileA fileB" REF_MSA=fileC THREADS=32 bash $0
EOF
  exit 0
fi

RAW="${WORK}/combined_raw.fasta"; EXREP="${WORK}/extraction_report.tsv"
if [[ ! -s "${RAW}" ]]; then
log "=== C3: extracting sequences ==="
python3 - "${IDS_TSV}" "${RAW}" "${EXREP}" "${EXPECT_LINEAGE}" "${ALLOW_MISSING}" "${REF_MSA}" ${MSA_FILES} <<'PY'
import sys, gzip, re, collections
ids_tsv, out_fa, rep, expect, allow, ref_msa = sys.argv[1:7]; files = sys.argv[7:]
expect, allow = int(expect), int(allow)
if ref_msa: files = files + [ref_msa]
norm = lambda s: re.sub(r"^(RS_|GB_)", "", s)
want, role = {}, {}
for l in list(open(ids_tsv))[1:]:
    i, r, _ = l.rstrip("\n").split("\t"); want[norm(i)] = i; role[i] = r
def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
found, lens = {}, collections.defaultdict(collections.Counter)
for p in files:
    left = {k for k in want if want[k] not in found}
    if not left: break
    head, buf = None, []
    def flush():
        if head is not None and norm(head) in want and want[norm(head)] not in found:
            seq = "".join(buf).upper(); found[want[norm(head)]] = (seq, p); lens[p][len(seq)] += 1
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(); head = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
        flush()
    print(f"read {p}: total found so far {len(found)} of {len(want)}")
def fail(m): sys.exit("CHECKPOINT FAILED: " + m)
byrole = collections.defaultdict(lambda: [0, 0])
for i, r in role.items():
    byrole[r][1] += 1; byrole[r][0] += i in found
with open(rep, "w") as fh:
    fh.write("id\trole\tfound\tsource\tlength\n")
    for i, r in role.items():
        s = found.get(i); fh.write(f"{i}\t{r}\t{int(s is not None)}\t{s[1] if s else ''}\t{len(s[0]) if s else ''}\n")
print("\nfound / wanted by role:"); [print(f"  {r:<14} {a} / {b}") for r, (a, b) in byrole.items()]
if "target" not in {role[i] for i in found}: fail("the target genome was not found in any alignment file")
if byrole["lineage"][0] < expect and not allow:
    miss = [i for i, r in role.items() if r == "lineage" and i not in found]
    fail(f"only {byrole['lineage'][0]} of {expect} lineage genomes found; missing e.g. {miss[:5]} (full list: {rep}). "
         "Add the alignment that holds them to MSA_FILES, or ALLOW_MISSING=1 to continue (then say so in the Methods).")
if byrole["outgroup"][0] < 2: fail("fewer than 2 outgroup sequences found")
if byrole["ref_neighbour"][0] == 0: print("WARNING: no GTDB reference neighbours found (REF_MSA missing or IDs differ); the tree will have catalogue neighbours only.")
all_len = collections.Counter(len(s) for s, _ in found.values())
if len(all_len) != 1:
    print("alignment lengths:", dict(all_len))
    for p, c in lens.items(): print("  ", p, dict(c))
    fail("sequences have different lengths: files are not the same alignment (one is probably already trimmed). Use untrimmed GTDB-Tk alignments for all files.")
with open(out_fa, "w") as fh:
    for i in role:
        if i in found: fh.write(f">{i}\n{found[i][0]}\n")
print(f"\nwrote {len(found)} sequences of length {list(all_len)[0]} -> {out_fa}")
PY
fi
log "C3 passed. ($(grep -c '>' "${RAW}") sequences) Report: ${EXREP}"
stop_if extract

FILT="${WORK}/combined_filtered.fasta"; TRIM="${WORK}/combined_trimmed.fasta"
if [[ ! -s "${FILT}" ]]; then
log "=== C4: filtering very gappy sequences (target is never dropped) ==="
python3 - "${RAW}" "${FILT}" "${IDS_TSV}" "${MIN_SEQ_COVER}" "${WORK}/filter_report.tsv" <<'PY'
import sys
raw, out, ids_tsv, thr, rep = sys.argv[1:6]; thr = float(thr)
role = {l.split("\t")[0]: l.split("\t")[1] for l in list(open(ids_tsv))[1:]}
recs, h, b = [], None, []
for line in open(raw):
    if line.startswith(">"):
        if h: recs.append((h, "".join(b)))
        h, b = line[1:].strip(), []
    else: b.append(line.strip())
if h: recs.append((h, "".join(b)))
kept, dropped = [], []
with open(rep, "w") as fr:
    fr.write("id\trole\tnon_gap_fraction\tkept\n")
    for h, s in recs:
        cov = sum(c not in "-.X" for c in s) / len(s)
        keep = cov >= thr or role[h] == "target"
        fr.write(f"{h}\t{role[h]}\t{cov:.3f}\t{int(keep)}\n")
        (kept if keep else dropped).append((h, s))
for h, s in dropped: print(f"  dropped (coverage < {thr}): {h} [{role[h]}]")
tc = [sum(c not in '-.X' for c in s)/len(s) for h, s in recs if role[h] == "target"][0]
print(f"target non-gap fraction: {tc:.3f}  (GTDB-Tk reported msa_percent 75.19)")
if sum(role[h] == "lineage" for h, _ in dropped) > 5: sys.exit("CHECKPOINT FAILED: >5 lineage genomes dropped for low coverage; check the alignment")
if sum(role[h] == "outgroup" for h, _ in kept) < 2: sys.exit("CHECKPOINT FAILED: fewer than 2 outgroup sequences left")
with open(out, "w") as fh:
    for h, s in kept: fh.write(f">{h}\n{s}\n")
print(f"kept {len(kept)} sequences, dropped {len(dropped)}; report {rep}")
PY
fi
log "C4 passed."
if [[ ! -s "${TRIM}" ]]; then
  log "=== C5: trimAl -automated1 ==="
  "${TRIMAL}" -in "${FILT}" -out "${TRIM}" -automated1
fi
NCOL=$(awk '/^>/{if(n++)exit;next}{c+=length($0)}END{print c}' "${TRIM}")
NSEQ=$(grep -c '>' "${TRIM}")
log "trimmed alignment: ${NSEQ} sequences x ${NCOL} columns"
[[ "${NCOL}" -ge "${MIN_TRIM_COLS}" ]] || die "only ${NCOL} columns left after trimAl (< ${MIN_TRIM_COLS}); inspect the alignment"
grep -q "^>${TARGET}$" "${TRIM}" || die "target missing from the trimmed alignment"
log "C5 passed."
stop_if trim

PREFIX="${WORK}/ml_tree"
OUTG=$(awk -F'\t' '$2=="outgroup"{print $1}' "${IDS_TSV}" | while read -r g; do grep -q "^>${g}$" "${TRIM}" && echo "$g"; done | paste -sd, -)
[[ -n "${OUTG}" ]] || die "no outgroup sequence in the trimmed alignment"
if [[ ! -s "${PREFIX}.treefile" ]]; then
  log "=== C6: IQ-TREE (ModelFinder + ${UFBOOT} UFBoot, seed ${SEED}) ==="
  log "outgroup: ${OUTG}"
  "${IQTREE}" -s "${TRIM}" -m MFP -B "${UFBOOT}" -bnni -T "${THREADS}" --seed "${SEED}" -o "${OUTG}" --prefix "${PREFIX}"
fi
[[ -s "${PREFIX}.treefile" ]] || die "IQ-TREE produced no treefile; see ${PREFIX}.log"
MODEL=$(grep -m1 "Best-fit model" "${PREFIX}.iqtree" 2>/dev/null | sed 's/.*: //' || true)
log "C6 passed. Model: ${MODEL:-'(see .iqtree)'}   <- fill into Methods"
stop_if tree

log "=== C7: tree checks ==="
REPORT="${WORK}/target_placement_report.txt"
Rscript - "${PREFIX}.treefile" "${IDS_TSV}" "${WORK}" "${LABEL_CSV}" "${UFBOOT}" <<'RS' | tee "${REPORT}"
suppressMessages(library(ape))
a <- commandArgs(TRUE); tf <- a[1]; ids <- read.delim(a[2], stringsAsFactors = FALSE); wd <- a[3]; lab <- a[4]
tr <- read.tree(tf)
role <- setNames(ids$role, ids$id); role <- role[names(role) %in% tr$tip.label]
cat("tips in tree:", Ntip(tr), "| roles:", paste(names(table(role)), table(role), collapse = "; "), "\n")
stopifnot("tree tips do not match the alignment ids" = setequal(tr$tip.label, names(role)))
tgt  <- names(role)[role == "target"]; lin <- names(role)[role == "lineage"]; og <- names(role)[role == "outgroup"]
cat("\n[outgroup] monophyletic in the unrooted tree:", is.monophyletic(tr, og), "\n")
if (!is.monophyletic(tr, og)) cat("  WARNING: outgroup is not monophyletic.\n")
if (is.monophyletic(tr, og)) { rt <- root(tr, outgroup = og, resolve.root = TRUE)
} else { cat("  Rooting on the single outgroup tip", og[1], "instead. Treat the rooted topology with caution.\n")
         rt <- root(tr, outgroup = og[1], resolve.root = TRUE) }
write.tree(rt, file.path(wd, "ml_tree_rooted_rawIDs.nwk"))
sup <- function(node) { l <- rt$node.label[node - Ntip(rt)]; if (is.null(l) || is.na(l) || l == "") NA else suppressWarnings(as.numeric(l)) }
cat("\n[lineage] the", length(lin), "UMGS1840/HGM11507 genomes monophyletic (excluding target):", is.monophyletic(rt, lin), "\n")
cat("[lineage + target] monophyletic:", is.monophyletic(rt, c(lin, tgt)), "\n")
m <- getMRCA(rt, c(lin, tgt)); tips_m <- extract.clade(rt, m)$tip.label
intr <- setdiff(tips_m, c(lin, tgt))
cat("MRCA clade of lineage + target: ", length(tips_m), "tips; support", sup(m), "; non-lineage tips inside:", length(intr), "\n")
if (length(intr)) cat("   ", paste(intr, "[", role[intr], "]", collapse = "\n    "), "\n")
p  <- rt$edge[rt$edge[, 2] == which(rt$tip.label == tgt), 1]
sis <- setdiff(extract.clade(rt, p)$tip.label, tgt)
cat("\n[target] parent node support:", sup(p), "| sister clade (", length(sis), "tips) role counts:\n"); print(table(role[sis]))
if (length(sis) <= 12) cat("   sister tips:", paste(sis, collapse = ", "), "\n")
d <- cophenetic.phylo(rt)[tgt, ]; d <- sort(d[names(d) != tgt])[1:12]
cat("\n[target] 12 nearest tips by ML patristic distance:\n")
print(data.frame(tip = names(d), role = role[names(d)], dist = round(d, 4)), row.names = FALSE)
ids$order <- ifelse(grepl(";o__", ids$classification), sub(".*;o__([^;]*);.*", "\\1", ids$classification), NA)
ord <- setNames(ids$order, ids$id)[names(role)]
cat("\n[order-level groups] (lineage genomes count as their GTDB order, UMGS1840):\n")
grp <- do.call(rbind, lapply(setdiff(unique(ord[!is.na(ord) & ord != ""]), NA), function(o) {
  tips <- names(ord)[which(ord == o)]; if (length(tips) < 2) return(NULL)
  n <- getMRCA(rt, tips); cl <- extract.clade(rt, n)$tip.label
  data.frame(order = o, tips = length(tips), monophyletic = is.monophyletic(rt, tips), support = sup(n),
             clade_size = length(cl), foreign_tips_inside = length(setdiff(cl, tips)))
}))
if (!is.null(grp)) print(grp[order(-grp$tips), ], row.names = FALSE)
cat("\n[target + each order]: smallest clade containing the target and that order (small foreign_tips = target nests close to it):\n")
tg <- do.call(rbind, lapply(setdiff(unique(ord[!is.na(ord) & ord != ""]), NA), function(o) {
  tips <- names(ord)[which(ord == o)]; if (length(tips) < 2) return(NULL)
  n <- getMRCA(rt, c(tips, tgt)); cl <- extract.clade(rt, n)$tip.label
  data.frame(order = o, clade_size = length(cl), support = sup(n), foreign_tips_inside = length(setdiff(cl, c(tips, tgt))))
}))
if (!is.null(tg)) print(tg[order(tg$clade_size), ], row.names = FALSE)
cat("\n[target] sister clade by order:\n"); print(table(ifelse(is.na(ord[sis]) | ord[sis] == "", "(unknown)", ord[sis])))
cat("\nHOW TO READ THIS (guidance, not a result):\n",
    " - a statement about the target is supported only where the deciding node has UFBoot >= 95; below that, call it unresolved;\n",
    " - the order table shows whether each order is monophyletic and how many foreign tips sit inside it;\n",
    " - the target+order table shows which order the target nests with most closely;\n",
    " - patristic distances above are NOT evidence of closeness (long branches distort them); use the topology and support.\n")
if (file.exists(lab)) {
  L <- read.csv(lab, stringsAsFactors = FALSE)
  if (all(c("genome", "public_id") %in% names(L))) {
    nm <- rt$tip.label; i <- match(nm, L$genome); nm[!is.na(i)] <- L$public_id[i[!is.na(i)]]; nm[nm == tgt] <- "TARGET"
    r2 <- rt; r2$tip.label <- nm; write.tree(r2, file.path(wd, "ml_tree_rooted_publicID.nwk"))
    cat("\nwrote ml_tree_rooted_publicID.nwk (lineage = public_id, target = 'TARGET', others raw)\n")
  }
} else cat("\n(LABEL_CSV not found; public-ID tree not written)\n")
write.table(data.frame(id = names(role), role = role), file.path(wd, "tip_roles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
RS
log "C7 done. Report: ${REPORT}"
log "-------------------------------------------------------"
log "Outputs in ${WORK}:"
log "  ml_tree.treefile / .iqtree          ML tree with UFBoot, model, log"
log "  ml_tree_rooted_rawIDs.nwk           rooted on the outgroup"
log "  ml_tree_rooted_publicID.nwk         lineage tips relabelled (public_id), target = TARGET"
log "  target_placement_report.txt         placement checks"
log "  tip_roles.tsv                       tip -> role (for an iTOL colour strip)"
log "  extraction_report.tsv, filter_report.tsv"
log "Methods facts to record: trimAl version, IQ-TREE version, model (${MODEL:-see .iqtree}), ${UFBOOT} UFBoot, seed ${SEED}, N tips, N columns (${NCOL})."
log "-------------------------------------------------------"
