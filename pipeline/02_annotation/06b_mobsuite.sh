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
#  SCRIPT  : 06b_mobsuite.sh  -  Mobility typing of plasmid contigs (MOB-suite)
#  AUTHOR  : QRQ
#  INPUT   : annotation/platon_merged.tsv
#            species genomes (*.fa)
#  OUTPUT  : annotation/mobsuite_output/mobtyper_all.tsv
# ==========================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"
exec > >(tee -a "${LOGS}/06b_mobsuite_$(date +%Y%m%d_%H%M%S).log") 2>&1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mobsuite

PLATON=${ANNOT}/platon_merged.tsv
OUT=${ANNOT}/mobsuite_output
PERGENOME=${OUT}/per_genome
CONTIGS=${OUT}/plasmid_contigs
WORK=${OUT}/work
mkdir -p "${PERGENOME}" "${CONTIGS}" "${WORK}"

MOB_THREADS=2
JOBS=$(( THREADS / MOB_THREADS ))
[[ "${JOBS}" -lt 1 ]] && JOBS=1
[[ "${JOBS}" -gt 32 ]] && JOBS=32

PRIMARY_CLUSTER_DIST=0.06

log "MOB-suite: $(mob_typer -V 2>&1 | head -n1)"
log "Parallel genomes: ${JOBS}  (threads each: ${MOB_THREADS})"

[[ -s "${PLATON}" ]] || { log "ERROR: ${PLATON} missing or empty."; exit 1; }

MOB_DB=$(python -c "import mob_suite, os; print(os.path.join(os.path.dirname(mob_suite.__file__),'databases'))")
if [[ ! -s "${MOB_DB}/status.txt" ]]; then
    log "ERROR: MOB-suite databases not initialized (${MOB_DB}/status.txt absent)."
    log "       Run 00d_mobsuite_db_init_manual.py first."
    exit 1
fi
log "MOB-suite databases: ${MOB_DB}"

log "Parsing Platon calls..."
CONTIG_MAP=${WORK}/genome_contig_map.tsv

awk -F'\t' 'NR>1 && $1!="" && $2!="" { print $1"\t"$2 }' "${PLATON}" \
    | sort -u > "${CONTIG_MAP}.tmp"
mv "${CONTIG_MAP}.tmp" "${CONTIG_MAP}"

N_ROWS=$(wc -l < "${CONTIG_MAP}")
N_GENOMES_EXPECTED=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "${PLATON}" | sort -u | wc -l)
N_GENOMES=$(cut -f1 "${CONTIG_MAP}" | sort -u | wc -l)
log "Plasmid-flagged contigs: ${N_ROWS} across ${N_GENOMES} genomes"
[[ "${N_ROWS}" -gt 0 ]] || { log "ERROR: no plasmid contigs parsed from ${PLATON}"; exit 1; }

if [[ "${N_GENOMES}" -ne "${N_GENOMES_EXPECTED}" ]]; then
    log "ERROR: genome count mismatch (${N_GENOMES} vs expected ${N_GENOMES_EXPECTED})."
    log "       CONTIG_MAP may be corrupt or truncated. Aborting rather than"
    log "       silently proceeding on a short list. Delete ${WORK} and rerun."
    exit 1
fi

cut -f1 "${CONTIG_MAP}" | sort -u > "${WORK}/genomes.txt.tmp"
mv "${WORK}/genomes.txt.tmp" "${WORK}/genomes.txt"

N_GENOMES_LIST=$(wc -l < "${WORK}/genomes.txt")
if [[ "${N_GENOMES_LIST}" -ne "${N_GENOMES}" ]]; then
    log "ERROR: genomes.txt has ${N_GENOMES_LIST} lines, expected ${N_GENOMES}."
    log "       Aborting rather than running xargs on a short/corrupt list."
    exit 1
fi
log "Verified genomes.txt: ${N_GENOMES_LIST} genomes (matches expected)"

if awk '!/__/ || length($0)<20 {print; bad=1} END{exit !bad}' "${WORK}/genomes.txt"; then
    log "ERROR: malformed genome ID(s) found in genomes.txt (printed above)."
    log "       This indicates a truncated/corrupt write. Aborting."
    exit 1
fi

cat > "${WORK}/worker.sh" <<'WORKER'
#!/usr/bin/env bash
set -euo pipefail
G="$1"
GENOME_DIR="$2"; GENOME_EXT="$3"; CONTIG_MAP="$4"
CONTIGS="$5"; PERGENOME="$6"; MOB_THREADS="$7"

RESULT="${PERGENOME}/${G}.mobtyper.txt"
# Resume: a completed run has a header plus at least one data row.
if [[ -s "${RESULT}" ]] && [[ "$(wc -l < "${RESULT}")" -gt 1 ]]; then
    exit 0
fi

FASTA="${GENOME_DIR}/${G}.${GENOME_EXT}"
if [[ ! -s "${FASTA}" ]]; then
    echo "WARN: genome FASTA not found: ${FASTA}" >&2
    exit 0
fi

IDS="${CONTIGS}/${G}.ids"
SUB="${CONTIGS}/${G}.plasmids.fasta"
awk -F'\t' -v g="${G}" '$1==g { print $2 }' "${CONTIG_MAP}" | sort -u > "${IDS}"

# Pull only the plasmid-flagged contigs. Header IDs are taken up to the
# first whitespace, which is how Platon reports them.
#
# The header is REWRITTEN as '>' id, discarding anything after the first
# whitespace. Several assemblies in this dataset (PRJCA008889,
# PRJNA657455, PRJNA972320 and others) carry annotations in the contig
# header, e.g.
#     >15-271096 total_depth=5088.32 sample_depths=4.6,0.5,...
# MOB-suite splits such headers on whitespace and emits the extra tokens
# as additional leading columns, shifting every field after sample_id by
# two positions. That silently corrupts predicted_mobility, all four
# biomarker columns and mash_neighbor_distance without any error, so the
# header must be normalised here rather than detected downstream.
awk 'NR==FNR { want[$0]=1; next }
     /^>/ { id=substr($1,2); keep=(id in want); if (keep) print ">" id; next }
     keep { print }' "${IDS}" "${FASTA}" > "${SUB}"

N=$(grep -c '^>' "${SUB}" || true)
if [[ "${N}" -eq 0 ]]; then
    echo "WARN: no contigs extracted for ${G} (ID mismatch?)" >&2
    exit 0
fi

# --multi: treat every sequence as an independent plasmid. Without it,
# MOB-suite would merge all contigs into one hypothetical plasmid.
mob_typer --multi \
    --infile "${SUB}" \
    --out_file "${RESULT}" \
    --num_threads "${MOB_THREADS}" \
    >/dev/null 2>&1 || echo "WARN: mob_typer failed on ${G}" >&2
WORKER
chmod +x "${WORK}/worker.sh"

log "Running mob_typer across ${N_GENOMES} genomes..."
xargs -a "${WORK}/genomes.txt" -I{} -P "${JOBS}" \
    "${WORK}/worker.sh" {} \
    "${GENOME_DIR}" "${GENOME_EXT}" "${CONTIG_MAP}" \
    "${CONTIGS}" "${PERGENOME}" "${MOB_THREADS}"

log "mob_typer stage complete."

N_DONE=$(find "${PERGENOME}" -maxdepth 1 -name "*.mobtyper.txt" | wc -l)
log "Per-genome result files: ${N_DONE} / ${N_GENOMES} expected"

MIN_COVERAGE_FRAC="0.90"
COVERAGE=$(awk -v d="${N_DONE}" -v e="${N_GENOMES}" 'BEGIN{printf "%.4f", d/e}')
COVERAGE_OK=$(awk -v c="${COVERAGE}" -v m="${MIN_COVERAGE_FRAC}" 'BEGIN{print (c>=m)?1:0}')
if [[ "${COVERAGE_OK}" -ne 1 ]]; then
    log "ERROR: only ${N_DONE}/${N_GENOMES} genomes produced results (${COVERAGE})."
    log "       This is below the ${MIN_COVERAGE_FRAC} completeness threshold."
    log "       Refusing to merge a partial run as if it were complete."
    log "       Check for skipped genomes with:"
    log "         comm -23 <(sort ${WORK}/genomes.txt) <(find ${PERGENOME} -maxdepth 1 -name '*.mobtyper.txt' -printf '%f\\n' | sed 's/\\.mobtyper\\.txt\$//' | sort)"
    log "       Re-run this script — completed genomes are skipped automatically (resume-safe)."
    exit 1
fi

MERGED=${OUT}/mobtyper_all.tsv
python3 - "$PERGENOME" "$MERGED" "$PRIMARY_CLUSTER_DIST" <<'PY'
import os, sys, csv

pergenome, merged, dist_cut = sys.argv[1], sys.argv[2], float(sys.argv[3])

files = sorted(f for f in os.listdir(pergenome) if f.endswith(".mobtyper.txt"))
header, rows = None, []
shifted = []

BIOMARKERS = ["rep_type(s)", "relaxase_type(s)", "mpf_type", "orit_type(s)"]

for fn in files:
    genome = fn[:-len(".mobtyper.txt")]
    path = os.path.join(pergenome, fn)
    with open(path, newline="") as fh:
        rd = list(csv.DictReader(fh, delimiter="\t"))
    if not rd:
        continue
    if header is None:
        header = list(rd[0].keys())
    for r in rd:
        r["genome"] = genome

        # Guard against field shift. MOB-suite emits extra leading columns
        # when a FASTA header contains whitespace-separated annotations,
        # displacing every field after sample_id and silently corrupting
        # predicted_mobility, the biomarkers and mash_neighbor_distance.
        # num_contigs must always be an integer; anything else means the
        # row is misaligned and must not be interpreted.
        nc = (r.get("num_contigs", "") or "").strip()
        if not nc.isdigit():
            shifted.append((genome, r.get("sample_id", "?")))
            continue

        # A biomarker hit is positive evidence of mobility machinery.
        has_marker = any(
            (r.get(b, "-") or "-").strip() not in ("-", "", "NA")
            for b in BIOMARKERS
        )

        # A mash distance at/above the primary cluster cutoff means the
        # nearest reference plasmid is not a genuine match, so the cluster
        # ID and host attribution on that row are noise.
        try:
            d = float(r.get("mash_neighbor_distance", "1") or 1)
        except ValueError:
            d = 1.0
        near_ref = d < dist_cut

        if has_marker or near_ref:
            conf = "resolved"
        else:
            # Distinguish "no evidence found" from "confirmed immobile".
            # MOB-suite's reference DB is dominated by clinical and
            # Enterobacteriaceae plasmids; rumen anaerobes (Bacteroidota,
            # Clostridia, Lachnospiraceae, Oscillospiraceae) are sparsely
            # represented, so an empty result here is frequently a database
            # blind spot rather than a biological finding. Rows marked
            # unresolved must NOT be counted as non-mobilizable.
            conf = "unresolved_db_limited"
        r["mobility_confidence"] = conf
        rows.append(r)

if header is None:
    sys.exit("ERROR: no mob_typer results found to merge.")

if shifted:
    bad = sorted({g for g, _ in shifted})
    listing = os.path.join(os.path.dirname(merged), "shifted_genomes.txt")
    with open(listing, "w") as fh:
        fh.write("\n".join(bad) + "\n")
    sys.exit(
        f"ERROR: {len(shifted)} rows across {len(bad)} genomes have shifted "
        f"columns (num_contigs is not an integer).\n"
        f"  Cause: contig FASTA headers contained whitespace-separated "
        f"annotations, so mob_typer emitted extra leading columns.\n"
        f"  Affected genomes written to: {listing}\n"
        f"  Fix: ensure the extraction step rewrites headers as '>' id "
        f"(see worker.sh), delete those genomes' .mobtyper.txt files, and "
        f"re-run — the script is resume-safe and will redo only those.\n"
        f"  Refusing to merge: mobility, biomarker and mash-distance values "
        f"on these rows are misread, not missing."
    )

out_cols = ["genome"] + header + ["mobility_confidence"]
with open(merged, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

n = len(rows)
res = sum(1 for r in rows if r["mobility_confidence"] == "resolved")
print(f"MERGE: {n} typed plasmid contigs from {len(files)} genomes")
print(f"MERGE: resolved={res}  unresolved_db_limited={n-res}")
for cls in ("conjugative", "mobilizable", "non-mobilizable"):
    c = sum(1 for r in rows
            if (r.get("predicted_mobility", "") or "").lower() == cls)
    cr = sum(1 for r in rows
             if (r.get("predicted_mobility", "") or "").lower() == cls
             and r["mobility_confidence"] == "resolved")
    print(f"MERGE: {cls}: {c} total, {cr} resolved")
PY

log "-------------------------------------------------------"
log "Merged table: ${MERGED}"
log "Plasmid-flagged contigs submitted: ${N_ROWS}"
log "Genomes with results:              ${N_DONE} / ${N_GENOMES}"
log "-------------------------------------------------------"
log "NOTE: rows flagged 'unresolved_db_limited' have neither a mobility"
log "      biomarker nor a genuine reference match (mash distance >= ${PRIMARY_CLUSTER_DIST})."
log "      Report these as unresolved, NOT as non-mobilizable — their"
log "      cluster_id and host attribution columns are not meaningful."
log "Next: join ${MERGED} to metadata.tsv on 'species' for the cross-species"
log "      comparison, and to the AMR table to test plasmid-borne ARG mobility."
