#!/usr/bin/env bash
# ============================================================
# 32b_catalogue_comparison.sh — RESULTS SECTION 3.2
#
# Compares the 25,451 species-level representatives against multiple
# reference sets with skani, giving several DIFFERENT novelty claims:
#
#   RGMGC (10,372 MAGs)  -> novel vs. the RGMGC gene/MAG-catalogue paper
#   RGMC (~20-40k MAGs)  -> novel vs. the largest, most recent global
#                           ruminant MAG catalogue (Zhang et al. 2026,
#                           Animal Research and One Health) -- NOT the
#                           same resource as RGMGC above, despite the
#                           similar name
#   Hungate1000 (~500)   -> novel vs. cultured rumen isolates
#   GTDB R232 reps       -> novel vs. all described prokaryotic diversity
#   HRGM2 (4,824 NC reps)-> novel vs. the HRGM2 human gut catalogue
#                           (Nature Microbiology 2026); run in both directions
#                           so shared/unique species can be reported from each
#                           side, as in the HRGM2 paper's overlap analysis
#   UPGG (pig gut)       -> novel vs. the UPGG pig gut genome catalogue
#                           (dRep 99% ANI representatives; forward direction
#                           only: it is a cross-host, non-ruminant reference)
#   MGnify (6 biomes)    -> novel vs. cow-rumen/sheep-rumen (same
#                           methodology as this catalogue) plus four
#                           cross-physiology gut catalogues
#
# This is a separate novelty axis from 32a's GTDB-Tk analysis: GTDB-Tk
# asks "does this genome fall inside a described taxon", skani asks "is
# there any genome in these collections within 95% ANI of it". Reviewers
# of catalogue papers generally expect the rumen-specific comparisons
# (RGMGC/RGMC/Hungate/MGnify rumen catalogues), which GTDB-Tk alone
# cannot provide.
#
# METHOD: skani sketch (reference) + skani search (queries), for every
# reference set except GTDB. For GTDB R232, this reuses GTDB-Tk's OWN
# pre-built skani sketch database (release232/skani/database/) directly
# — the same index GTDB-Tk itself queries during ani_screening — so no
# FASTA download or sketching is needed for that comparison at all.
#
# NOTE: uses `find ... > list` + skani's -l/--ql flags throughout, never
# shell globs. Genome counts here (25k queries, 113k refs) are well past
# the point where `ls dir/*.fa` fails with "Argument list too long" — the
# recurring bug pattern in this project.
# ============================================================
set -euo pipefail

QUERY_DIR=/disk252/qamar/metagenome_analysis/drep_analysis/drep_95_species_final/dereplicated_genomes
RGMGC_DIR=/disk252/qamar/disk2_data/RGMGC_10000_MAG
RGMC_DIR=/disk252/qamar/databases/rgmc/dereplicated_genomes
HUNGATE_DIR=/disk252/qamar/disk2_data/Hungate1000_Collection/HungateCollection
GTDB_DIR=/disk252/qamar/databases/gtdbtk_database/release232
# UPGG pig gut catalogue: dRep 99% ANI dereplicated genomes, flat folder of *.fa
UPGG_DIR=${UPGG_DIR:-/disk202/lius/microbiome/11.dRep/dRep99_all/dereplicated_genomes}

OUT=/disk252/qamar/metagenome_analysis/annotation/results_phase3/s32_taxonomy_novelty/skani
# TEST MODE:  TEST_N=300 bash 32b_catalogue_comparison.sh
#   Runs only the HRGM2 comparison, on a fixed random subset of TEST_N of our
#   representatives, into a separate folder (${OUT}_hrgm2_test), so the real
#   outputs are never touched. Forward results (our subset vs all 4,824 HRGM2
#   reps) are valid; the reverse/two-sided numbers are NOT (reverse searches
#   against the subset only). Use it to check paths, columns and runtime.
# ONLY_REF=upgg  runs just one reference set (names: rgmgc rgmc hungate gtdb
#   hrgm2 upgg) into the normal output folder, leaving finished results alone.
#   With TEST_N it goes to ${OUT}_<ref>_test instead. TEST_N alone = hrgm2.
TEST_N=${TEST_N:-0}
ONLY_REF=${ONLY_REF:-}
HRGM2_ONLY=${HRGM2_ONLY:-0}
[[ "${HRGM2_ONLY}" -eq 1 ]] && ONLY_REF=hrgm2
if [[ "${TEST_N}" -gt 0 ]]; then
  ONLY_REF=${ONLY_REF:-hrgm2}
  OUT=${OUT}_${ONLY_REF}_test
fi
WORK=${OUT}/work
LOGS=/disk252/qamar/metagenome_analysis/annotation/logs
mkdir -p "${OUT}" "${WORK}" "${LOGS}"

exec > >(tee -a "${LOGS}/32b_catalogue_comparison_$(date +%Y%m%d_%H%M%S).log") 2>&1
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

THREADS=${THREADS:-48}
# Species boundary. A query with a >=95% ANI hit at >=50% aligned fraction
# to a reference is the SAME species as something already in that
# collection, i.e. not novel relative to it.
MIN_ANI=95
MIN_AF=50
# HRGM2 defines a conspecific pair as ANI > 95% and coverage > 0.6 (HRGM2
# Methods, "Identifying NC-genome species overlap"). Use 60 for that reference
# set only, so the overlap numbers line up with the HRGM2 paper. Every other
# reference set stays at MIN_AF.
HRGM2_MIN_AF=60
# skani screening floor: report hits down to 80% ANI so genus-level
# proximity can also be examined, not just the 95% species call.
SCREEN_ANI=80

command -v skani >/dev/null 2>&1 || {
  log "ERROR: skani not on PATH. Activate an env that provides it, e.g.:"
  log "  conda activate gtdbtk-2.7.2   # skani ships with GTDB-Tk"
  exit 1
}
log "skani: $(skani --version 2>&1 | head -n1)"

# ============================================================
# STEP 1 — build genome file lists
# ============================================================
log "Building genome lists..."

find "${QUERY_DIR}" -maxdepth 1 -name "*.fa" > "${WORK}/query_list.txt"
N_QUERY=$(wc -l < "${WORK}/query_list.txt")
log "Query genomes (95% ANI representatives): ${N_QUERY}"
[[ "${N_QUERY}" -gt 0 ]] || { log "ERROR: no query genomes found in ${QUERY_DIR}"; exit 1; }
if [[ "${TEST_N}" -gt 0 ]]; then
  # fixed seed, so the same subset is drawn on every rerun
  shuf -n "${TEST_N}" --random-source=<(yes 42) "${WORK}/query_list.txt" > "${WORK}/query_list.sub"
  mv "${WORK}/query_list.sub" "${WORK}/query_list.txt"
  N_QUERY=$(wc -l < "${WORK}/query_list.txt")
  log "TEST MODE: using ${N_QUERY} random query genomes; only ${ONLY_REF} will be run"
fi

# RGMGC: <GCA>/ncbi_dataset/data/<GCA>/*_genomic.fna
find "${RGMGC_DIR}" -name "*_genomic.fna" -o -name "*_genomic.fna.gz" \
  | sort -u > "${WORK}/rgmgc_list.txt"
N_RGMGC=$(wc -l < "${WORK}/rgmgc_list.txt")
log "RGMGC reference genomes: ${N_RGMGC}"

# RGMC (Zhang et al. 2026, Animal Research and One Health): a separate,
# much more recent MAG catalogue (40,812 strain-level / 20,530
# species-level non-redundant genomes, 18 hosts, 27 countries) -- NOT
# the same resource as RGMGC above. Flat directory, one genome per file,
# named RGMC.<id>.fa -- confirmed via direct testing, no nested
# structure like RGMGC's NCBI-datasets-tool layout.
find "${RGMC_DIR}" -maxdepth 1 -name "*.fa" \
  | sort -u > "${WORK}/rgmc_list.txt"
N_RGMC=$(wc -l < "${WORK}/rgmc_list.txt")
log "RGMC reference genomes: ${N_RGMC}"

# UPGG pig gut catalogue (dRep 99% ANI representatives, so strain-level:
# several genomes may fall in one 95% species, which is fine for a best-hit call).
find "${UPGG_DIR}" -maxdepth 1 \( -name "*.fa" -o -name "*.fa.gz" \) \
  | sort -u > "${WORK}/upgg_list.txt"
N_UPGG=$(wc -l < "${WORK}/upgg_list.txt")
log "UPGG reference genomes: ${N_UPGG}"

# Hungate1000: <STRAIN>/QC_and_Genome_Assembly/submission.assembly.fasta
find "${HUNGATE_DIR}" -name "submission.assembly.fasta" \
  | sort -u > "${WORK}/hungate_list.txt"
N_HUNGATE=$(wc -l < "${WORK}/hungate_list.txt")
log "Hungate1000 reference genomes: ${N_HUNGATE}"

# GTDB R232 comparison reuses GTDB-Tk's OWN pre-built skani sketch
# database at release232/skani/database/ (sketches.db, index.db,
# markers.bin) — the exact index GTDB-Tk itself queries during its
# ani_screening step. This needs no FASTA download and no sketching: it
# is already a valid `skani sketch` output directory, ready for
# `skani search -d <this dir>`. Confirmed present (80GB) on this server;
# using it is both faster AND more correct than falling back to raw R226
# genome FASTAs, since it matches the R232 release actually used for
# classification in 32a.
GTDB_SKETCH_DB="${GTDB_DIR}/skani/database"
if [[ ! -s "${GTDB_SKETCH_DB}/sketches.db" ]]; then
    log "WARNING: expected pre-built skani sketch DB not found at ${GTDB_SKETCH_DB}"
    log "         (looked for sketches.db). GTDB comparison will be SKIPPED."
    GTDB_SKETCH_DB=""
else
    log "GTDB R232 comparison will reuse the pre-built skani sketch DB at ${GTDB_SKETCH_DB}"
    log "  ($(du -sh "${GTDB_SKETCH_DB}" 2>/dev/null | cut -f1) — no sketching needed for this reference set)"
fi
: > "${WORK}/gtdb_list.txt"   # empty fallback ref list, only used if GTDB_SKETCH_DB is unset above

# ============================================================
# STEP 2 — sketch + search against each reference set
# ============================================================
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
        # Reuse an existing skani sketch database as-is (e.g. GTDB-Tk's
        # own R232 index) — no genome list needed, no sketching step.
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

# ------------------------------------------------------------------------
# MGnify biome-specific catalogues (species_catalogue/ downloads).
# cow-rumen and sheep-rumen are the PRIMARY additions: same rumen
# environment, same dRep/GTDB-Tk-based construction methodology as this
# catalogue, so a match/no-match result is directly comparable in a way
# GTDB or RGMGC (built differently) is not. human-gut/pig-gut/mouse-gut/
# chicken-gut are cross-physiology references — informative for the
# Discussion, not the core novelty claim, and only run if present.
# ------------------------------------------------------------------------
MGNIFY_ROOT=/disk252/qamar/databases/mgnify_datasets

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

# ------------------------------------------------------------------------
# HRGM2 (Nature Microbiology 2026): 4,824 species-level representatives of
# near-complete genomes (>=90% complete, <=5% contaminated), 41 countries.
# Extract HRGMv2_Rep_Genome.tar.gz first and point HRGM2_DIR at the folder.
# Forward: our reps vs HRGM2 reps (feeds the per-genome table below).
# Reverse: HRGM2 reps vs our reps (feeds the two-sided overlap summary).
# ------------------------------------------------------------------------
HRGM2_TAR=/disk252/qamar/databases/urgm2/HRGMv2_Rep_Genome.tar.gz
HRGM2_DIR=/disk252/qamar/databases/urgm2/HRGMv2_Rep_Genome

# Extract once. The tarball's top-level folder is HRGMv2_Rep_Genome/ (checked with
# tar -tzf), with .fna files in nested HRGMv2_17XX/HRGMv2_170X/ subfolders, so
# extract into the parent directory; find below searches recursively.
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

    # Reverse direction. Output name deliberately does NOT match skani_vs_*.tsv,
    # so the summary above does not pick it up as a reference set.
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

# ============================================================
# STEP 3 — per-query best hit + novelty call per reference set
# ============================================================
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
        # skani column names vary slightly by version; resolve them once.
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

# Per-genome table
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

# ------------------------------------------------------------------------
# HRGM2 two-sided overlap (same logic as the HRGM2 paper's Fig. 1h):
# a species is shared if its representative has a conspecific match
# (ANI >= 95, coverage >= 0.6) on the other side.
# ------------------------------------------------------------------------
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

# HRGM2 species absent from our catalogue: the "not shared" side, for review
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
