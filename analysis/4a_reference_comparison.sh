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
#  SCRIPT  : 4a_reference_comparison.sh  -  Matches to ten reference collections (skani, 95% ANI, 50% AF)
#  AUTHOR  : QRQ
#  INPUT   : species genomes
#            reference collections ($URGG_DB, $URGG_EXTDATA, $URGG_UPGG)
#  OUTPUT  : annotation/results_phase3/s32_taxonomy_novelty/skani/
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"
export URGG_EXTDATA="${URGG_EXTDATA:-/path/to/external_catalogues}"
export URGG_UPGG="${URGG_UPGG:-/path/to/UPGG_dRep99}"

set -euo pipefail

QUERY_DIR=${URGG_PROJECT}/drep_analysis/drep_95_species_final/dereplicated_genomes
RGMGC_DIR=${URGG_EXTDATA}/RGMGC_10000_MAG
RGMC_DIR=${URGG_DB}/rgmc/dereplicated_genomes
HUNGATE_DIR=${URGG_EXTDATA}/Hungate1000_Collection/HungateCollection
GTDB_DIR=${URGG_DB}/gtdbtk_database/release232
UPGG_DIR=${UPGG_DIR:-${URGG_UPGG}/dereplicated_genomes}

OUT=${URGG_PROJECT}/annotation/results_phase3/s32_taxonomy_novelty/skani
TEST_N=${TEST_N:-0}
ONLY_REF=${ONLY_REF:-}
HRGM2_ONLY=${HRGM2_ONLY:-0}
[[ "${HRGM2_ONLY}" -eq 1 ]] && ONLY_REF=hrgm2
if [[ "${TEST_N}" -gt 0 ]]; then
  ONLY_REF=${ONLY_REF:-hrgm2}
  OUT=${OUT}_${ONLY_REF}_test
fi
WORK=${OUT}/work
LOGS=${URGG_PROJECT}/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${LOGS}"

exec > >(tee -a "${LOGS}/32b_catalogue_comparison_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-48}
MIN_ANI=95
MIN_AF=50
HRGM2_MIN_AF=60
SCREEN_ANI=80

command -v skani >/dev/null 2>&1 || {
  log "ERROR: skani not on PATH."
  exit 1
}
log "skani: $(skani --version 2>&1 | head -n1)"

log "Building genome lists..."

find "${QUERY_DIR}" -maxdepth 1 -name "*.fa" > "${WORK}/query_list.txt"
N_QUERY=$(wc -l < "${WORK}/query_list.txt")
log "Query genomes (95% ANI representatives): ${N_QUERY}"
[[ "${N_QUERY}" -gt 0 ]] || { log "ERROR: no query genomes found in ${QUERY_DIR}"; exit 1; }
if [[ "${TEST_N}" -gt 0 ]]; then
  shuf -n "${TEST_N}" --random-source=<(yes 42) "${WORK}/query_list.txt" > "${WORK}/query_list.sub"
  mv "${WORK}/query_list.sub" "${WORK}/query_list.txt"
  N_QUERY=$(wc -l < "${WORK}/query_list.txt")
  log "TEST MODE: using ${N_QUERY} random query genomes; only ${ONLY_REF} will be run"
fi

find "${RGMGC_DIR}" -name "*_genomic.fna" -o -name "*_genomic.fna.gz" \
  | sort -u > "${WORK}/rgmgc_list.txt"
N_RGMGC=$(wc -l < "${WORK}/rgmgc_list.txt")
log "RGMGC reference genomes: ${N_RGMGC}"

find "${RGMC_DIR}" -maxdepth 1 -name "*.fa" \
  | sort -u > "${WORK}/rgmc_list.txt"
N_RGMC=$(wc -l < "${WORK}/rgmc_list.txt")
log "RGMC reference genomes: ${N_RGMC}"

find "${UPGG_DIR}" -maxdepth 1 \( -name "*.fa" -o -name "*.fa.gz" \) \
  | sort -u > "${WORK}/upgg_list.txt"
N_UPGG=$(wc -l < "${WORK}/upgg_list.txt")
log "UPGG reference genomes: ${N_UPGG}"

find "${HUNGATE_DIR}" -name "submission.assembly.fasta" \
  | sort -u > "${WORK}/hungate_list.txt"
N_HUNGATE=$(wc -l < "${WORK}/hungate_list.txt")
log "Hungate1000 reference genomes: ${N_HUNGATE}"

GTDB_SKETCH_DB="${GTDB_DIR}/skani/database"
if [[ ! -s "${GTDB_SKETCH_DB}/sketches.db" ]]; then
    log "WARNING: expected pre-built skani sketch DB not found at ${GTDB_SKETCH_DB}"
    log "         (looked for sketches.db). GTDB comparison will be SKIPPED."
    GTDB_SKETCH_DB=""
else
    log "GTDB R232 comparison will reuse the pre-built skani sketch DB at ${GTDB_SKETCH_DB}"
    log "  ($(du -sh "${GTDB_SKETCH_DB}" 2>/dev/null | cut -f1) - no sketching needed for this reference set)"
fi
: > "${WORK}/gtdb_list.txt"

run_comparison() {
    local name="$1" ref_list="$2" prebuilt_sketch="${3:-}"
    local out_tsv="${OUT}/skani_vs_${name}.tsv"

    if [[ -n "${ONLY_REF}" && "${name}" != "${ONLY_REF}" ]]; then
        return 0
    fi

    if [[ -s "${out_tsv}" ]]; then
        log "SKIP ${name}: ${out_tsv} already exists (delete it to force a re-run)."
        return 0
    fi

    local sketch_dir
    if [[ -n "${prebuilt_sketch}" ]]; then
        sketch_dir="${prebuilt_sketch}"
        log "=== ${name}: using pre-built sketch DB at ${sketch_dir} (no sketching needed) ==="
    else
        local n_ref; n_ref=$(wc -l < "${ref_list}")
        if [[ "${n_ref}" -eq 0 ]]; then
            log "SKIP ${name}: reference list is empty."
            return 0
        fi
        sketch_dir="${WORK}/${name}_sketch"
        log "=== ${name}: sketching ${n_ref} reference genomes ==="
        rm -rf "${sketch_dir}"
        skani sketch -l "${ref_list}" -o "${sketch_dir}" -t "${THREADS}"
    fi

    log "=== ${name}: searching ${N_QUERY} queries against the sketch ==="
    skani search -d "${sketch_dir}" \
        --ql "${WORK}/query_list.txt" \
        -o "${out_tsv}" \
        -t "${THREADS}" \
        --min-af "${MIN_AF}" \
        -s "${SCREEN_ANI}"

    log "${name} complete: $(wc -l < "${out_tsv}") hit rows -> ${out_tsv}"
}

run_comparison "rgmgc"   "${WORK}/rgmgc_list.txt"
run_comparison "rgmc"    "${WORK}/rgmc_list.txt"
run_comparison "upgg"    "${WORK}/upgg_list.txt"
run_comparison "hungate" "${WORK}/hungate_list.txt"
run_comparison "gtdb"    "${WORK}/gtdb_list.txt" "${GTDB_SKETCH_DB}"

MGNIFY_ROOT=${URGG_DB}/mgnify_datasets

for biome_version in "sheep-rumen/v1.0" "cow-rumen/v1.0.1" \
                     "human-gut/v2.0.2" "pig-gut/v1.0" \
                     "mouse-gut/v1.0" "chicken-gut/v1.0.1" ; do
    biome="${biome_version%%/*}"
    biome_dir="${MGNIFY_ROOT}/${biome_version}/species_catalogue"
    if [[ ! -d "${biome_dir}" ]]; then
        log "SKIP mgnify_${biome}: ${biome_dir} not found (not downloaded yet)"
        continue
    fi
    list_file="${WORK}/mgnify_${biome//-/_}_list.txt"
    find "${biome_dir}" -name "*.fna" -o -name "*.fna.gz" | sort -u > "${list_file}"
    n=$(wc -l < "${list_file}")
    if [[ "${n}" -eq 0 ]]; then
        log "SKIP mgnify_${biome}: 0 genomes found under ${biome_dir}"
        continue
    fi
    log "mgnify_${biome}: ${n} species representatives found"
    run_comparison "mgnify_${biome//-/_}" "${list_file}"
done

HRGM2_TAR=${URGG_DB}/urgm2/HRGMv2_Rep_Genome.tar.gz
HRGM2_DIR=${URGG_DB}/urgm2/HRGMv2_Rep_Genome

if [[ ! -d "${HRGM2_DIR}" && -s "${HRGM2_TAR}" ]]; then
    log "Extracting ${HRGM2_TAR} -> $(dirname "${HRGM2_DIR}")/"
    tar -xzf "${HRGM2_TAR}" -C "$(dirname "${HRGM2_DIR}")"
fi

if [[ -n "${ONLY_REF}" && "${ONLY_REF}" != "hrgm2" ]]; then
    log "SKIP hrgm2 block: ONLY_REF=${ONLY_REF}"
elif [[ -d "${HRGM2_DIR}" ]]; then
    find "${HRGM2_DIR}" \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \
         -o -name "*.fa.gz" -o -name "*.fna.gz" -o -name "*.fasta.gz" \) \
         | sort -u > "${WORK}/hrgm2_list.txt"
    N_HRGM2=$(wc -l < "${WORK}/hrgm2_list.txt")
    log "HRGM2 representative genomes: ${N_HRGM2} (paper: 4,824)"
    [[ "${N_HRGM2}" -eq 4824 ]] || log "WARNING: expected 4,824 HRGM2 representatives, found ${N_HRGM2}. Check the extraction."
    run_comparison "hrgm2" "${WORK}/hrgm2_list.txt"

    REV_OUT=${OUT}/hrgm2_reverse_vs_ours.tsv
    if [[ -s "${REV_OUT}" ]]; then
        log "SKIP reverse HRGM2: ${REV_OUT} already exists (delete it to force a re-run)."
    else
        OURS_SKETCH="${WORK}/ours_sketch"
        if [[ ! -s "${OURS_SKETCH}/sketches.db" ]]; then
            log "=== sketching ${N_QUERY} of our representatives (for reverse HRGM2 search) ==="
            rm -rf "${OURS_SKETCH}"
            skani sketch -l "${WORK}/query_list.txt" -o "${OURS_SKETCH}" -t "${THREADS}"
        fi
        log "=== reverse: ${N_HRGM2} HRGM2 representatives vs our ${N_QUERY} ==="
        skani search -d "${OURS_SKETCH}" \
            --ql "${WORK}/hrgm2_list.txt" \
            -o "${REV_OUT}" \
            -t "${THREADS}" \
            --min-af "${MIN_AF}" \
            -s "${SCREEN_ANI}"
        log "reverse HRGM2 complete: $(wc -l < "${REV_OUT}") hit rows -> ${REV_OUT}"
    fi
else
    log "SKIP hrgm2: ${HRGM2_DIR} not found (extract HRGMv2_Rep_Genome.tar.gz and set HRGM2_DIR)"
fi

log "Summarising best hits..."

python3 - "${OUT}" "${WORK}/query_list.txt" "${MIN_ANI}" "${MIN_AF}" "${HRGM2_MIN_AF}" <<'PY'
import os, sys, csv, glob

out_dir, query_list, min_ani, min_af = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
hrgm2_af = float(sys.argv[5])   # HRGM2's own coverage cut-off (0.6 in its paper)

def basename_noext(p):
    b = os.path.basename(p)
    for ext in (".fa", ".fna", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b

queries = [basename_noext(l.strip()) for l in open(query_list) if l.strip()]
summary = {q: {} for q in queries}
ref_sets = []

for path in sorted(glob.glob(os.path.join(out_dir, "skani_vs_*.tsv"))):
    ref = os.path.basename(path)[len("skani_vs_"):-len(".tsv")]
    ref_sets.append(ref)
    best = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        def pick(*cands):
            for c in cands:
                if c in cols: return c
            return None
        c_q   = pick("Query_file", "Query_name", "Query")
        c_ani = pick("ANI")
        c_af  = pick("Align_fraction_query", "Align_fraction_ref", "Align_fraction")
        c_r   = pick("Ref_file", "Ref_name", "Ref")
        if not all([c_q, c_ani, c_af]):
            sys.exit(f"ERROR: unexpected skani columns in {path}: {cols}")
        for row in rd:
            q = basename_noext(row[c_q])
            try:
                ani = float(row[c_ani]); af = float(row[c_af])
            except (ValueError, TypeError):
                continue
            if q not in best or ani > best[q][0]:
                best[q] = (ani, af, row.get(c_r, ""))
    for q in queries:
        if q in best:
            ani, af, rname = best[q]
            af_cut = hrgm2_af if ref == "hrgm2" else min_af
            summary[q][ref] = (ani, af, rname, (ani >= min_ani and af >= af_cut))
        else:
            summary[q][ref] = (0.0, 0.0, "", False)

out_csv = os.path.join(out_dir, "novelty_vs_reference_catalogues.csv")
with open(out_csv, "w", newline="") as fh:
    cols = ["genome"]
    for r in ref_sets:
        cols += [f"{r}_best_ani", f"{r}_best_af", f"{r}_best_hit", f"{r}_matched"]
    cols.append("novel_vs_all")
    w = csv.writer(fh); w.writerow(cols)
    for q in queries:
        row = [q]; matched_any = False
        for r in ref_sets:
            ani, af, rname, matched = summary[q][r]
            row += [f"{ani:.2f}", f"{af:.2f}", rname, int(matched)]
            matched_any = matched_any or matched
        row.append(int(not matched_any))
        w.writerow(row)

print(f"\n=== Novelty vs. reference catalogues (>= {min_ani}% ANI and AF >= {min_af}%; HRGM2 AF >= {hrgm2_af}% = match) ===")
print(f"Total query genomes: {len(queries):,}\n")
for r in ref_sets:
    n_match = sum(1 for q in queries if summary[q][r][3])
    n_novel = len(queries) - n_match
    print(f"{r:>8}: {n_match:>6,} matched ({100*n_match/len(queries):5.1f}%) | "
          f"{n_novel:>6,} NOVEL ({100*n_novel/len(queries):5.1f}%)")

n_novel_all = sum(1 for q in queries if not any(summary[q][r][3] for r in ref_sets))
print(f"\nNovel vs. ALL reference sets simultaneously: {n_novel_all:,} "
      f"({100*n_novel_all/len(queries):.1f}%)")
print(f"\nPer-genome table: {out_csv}")
PY

if [[ -s "${OUT}/skani_vs_hrgm2.tsv" && -s "${OUT}/hrgm2_reverse_vs_ours.tsv" ]]; then
python3 - "${OUT}" "${WORK}/query_list.txt" "${WORK}/hrgm2_list.txt" "${MIN_ANI}" "${HRGM2_MIN_AF}" <<'PY'
import os, sys, csv

out_dir, q_list, h_list = sys.argv[1], sys.argv[2], sys.argv[3]
min_ani, min_af = float(sys.argv[4]), float(sys.argv[5])

def bn(p):
    b = os.path.basename(p.strip())
    for ext in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if b.endswith(ext):
            return b[: -len(ext)]
    return b

def best_hits(path):
    best = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q, r = bn(row["Query_file"]), bn(row["Ref_file"])
            ani, af = float(row["ANI"]), float(row["Align_fraction_query"])
            if ani >= min_ani and af >= min_af and (q not in best or ani > best[q][0]):
                best[q] = (ani, af, r)
    return best

ours = [bn(l) for l in open(q_list) if l.strip()]
hrgm2 = [bn(l) for l in open(h_list) if l.strip()]
fwd = best_hits(os.path.join(out_dir, "skani_vs_hrgm2.tsv"))        # ours -> HRGM2
rev = best_hits(os.path.join(out_dir, "hrgm2_reverse_vs_ours.tsv"))  # HRGM2 -> ours

shared_ours = [q for q in ours if q in fwd]
shared_hrgm2 = [h for h in hrgm2 if h in rev]

unique_hrgm2 = [h for h in hrgm2 if h not in rev]
with open(os.path.join(out_dir, "hrgm2_species_not_in_our_catalogue.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["hrgm2_representative"])
    for h in unique_hrgm2:
        w.writerow([h])

print("\n=== HRGM2 overlap (ANI >= %s, coverage >= %s%%, both directions) ===" % (min_ani, min_af))
print(f"Our representatives with an HRGM2 species match: {len(shared_ours):,} of {len(ours):,} "
      f"({100*len(shared_ours)/len(ours):.1f}%); novel to HRGM2: {len(ours)-len(shared_ours):,}")
print(f"HRGM2 representatives with a match in our catalogue: {len(shared_hrgm2):,} of {len(hrgm2):,} "
      f"({100*len(shared_hrgm2)/len(hrgm2):.1f}%); HRGM2-only: {len(unique_hrgm2):,}")
print("The two sides differ because species clusters are not one-to-one between catalogues.")
print(f"HRGM2-only list: {os.path.join(out_dir, 'hrgm2_species_not_in_our_catalogue.csv')}")
PY
fi

if [[ "${TEST_N}" -gt 0 && "${ONLY_REF}" == "hrgm2" ]]; then
  log "TEST MODE: ignore the 'HRGM2 representatives with a match' and 'HRGM2-only' numbers above;"
  log "           the reverse search only saw ${N_QUERY} of our genomes. The forward numbers are valid for this subset."
fi
log "-------------------------------------------------------"
log "Done. Outputs in ${OUT}"
log "  skani_vs_<ref>.tsv                     raw skani hits per reference set"
log "  novelty_vs_reference_catalogues.csv    per-genome best hit + novelty calls"
log "-------------------------------------------------------"
