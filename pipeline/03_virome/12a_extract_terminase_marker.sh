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
#  SCRIPT  : 12a_extract_terminase_marker.sh  -  Extract terminase large subunit (TerL) sequences
#  AUTHOR  : QRQ
#  INPUT   : DRAM-v annotations and proteins
#            $URGG_DB/dramv_config/database_download/vog.annotations.tsv.gz
#  OUTPUT  : phase_5_virome/12_phylogeny/terminase_sequences.faa
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_DB="${URGG_DB:-/path/to/databases}"

source "$(dirname "$0")/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
PHYLO_DIR="${PHASE5_DIR}/12_phylogeny"
mkdir -p "${PHYLO_DIR}"

VOG_ANNOT_GZ="${URGG_DB}/dramv_config/database_download/vog.annotations.tsv.gz"
if [[ ! -f "${VOG_ANNOT_GZ}" ]]; then
    log "ERROR: ${VOG_ANNOT_GZ} not found. Adjust path if your install differs."
    exit 1
fi

log "Discovering terminase large subunit VOGs from VOGDB annotations..."
MARKER_VOGS="${PHYLO_DIR}/terminase_marker_vogs.tsv"
zcat "${VOG_ANNOT_GZ}" | awk -F'\t' '
    NR>1 && tolower($5) ~ /terminase/ && tolower($5) ~ /large/ { print }
' > "${MARKER_VOGS}"

n_markers=$(wc -l < "${MARKER_VOGS}")
log "  ${n_markers} terminase-large-subunit VOGs discovered -> ${MARKER_VOGS}"
if [[ "${n_markers}" -eq 0 ]]; then
    log "ERROR: zero terminase-large-subunit VOGs found. Check the description"
    log "       search pattern, or inspect vog.annotations.tsv manually for the"
    log "       exact wording used for terminase entries in this VOGDB version."
    exit 1
fi
cut -f1,5 "${MARKER_VOGS}" | head -10

log "Locating DRAM-v protein FASTA files..."
PROTEIN_FILES=$(find "${DRAMV_OUT}" -maxdepth 3 -iname "*.faa" 2>/dev/null)
n_faa=$(echo "${PROTEIN_FILES}" | grep -c . || echo 0)
log "  ${n_faa} protein FASTA files found under ${DRAMV_OUT}"
if [[ "${n_faa}" -eq 0 ]]; then
    log "ERROR: no .faa files found. DRAM-v's protein output may be named/located"
    log "       differently than assumed. Find it manually with:"
    log "         find ${DRAMV_OUT} -iname '*.faa' -o -iname '*protein*'"
    log "       then adjust this script's PROTEIN_FILES glob."
    exit 1
fi
PROTEIN_FILES_LIST="${PHYLO_DIR}/.protein_files_list.txt"
echo "${PROTEIN_FILES}" > "${PROTEIN_FILES_LIST}"

python3 << PYEOF
import csv, glob, os

# ---- load marker VOG accessions ----
marker_vogs = set()
with open("${MARKER_VOGS}") as f:
    for line in f:
        marker_vogs.add(line.split('\t')[0])
print(f"Using {len(marker_vogs)} terminase-large-subunit marker VOGs")

# ---- confirmed vOTU representatives ----
votu_reps = {}
with open("${VOTU_OUT}/vOTU_index.tsv") as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        votu_reps[row["representative"]] = row["votu_id"]

def reduce_to_base(s):
    for sep in ("__vs2", "__vib"):
        idx = s.find(sep)
        if idx != -1:
            return s[:idx]
    return s

rep_roots = {rep: reduce_to_base(rep) for rep in votu_reps}

# ---- scan per-batch annotation files (using the VERIFIED WORKING "scaffold"
#      column, not "fasta" -- see 05_prophage_flagging.sh comments) for
#      genes hitting a marker VOG, on a confirmed vOTU representative ----
terminase_gene_ids = {}  # votu_id -> gene_id (first hit found, one per vOTU)
batch_files = sorted(glob.glob("${DRAMV_OUT}/batch_*/annotate/annotations.tsv"))
print(f"Scanning {len(batch_files)} annotation files for marker gene hits...")

for fpath in batch_files:
    with open(fpath, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        if "scaffold" not in header or "vogdb_id" not in header or "gene_position" not in header:
            continue
        scaffold_idx = header.index("scaffold")
        vogdb_idx = header.index("vogdb_id")
        # gene identifier: reconstruct as "<scaffold>_<gene_position>" to match .faa headers
        gene_pos_idx = header.index("gene_position")
        for row in reader:
            if len(row) <= max(scaffold_idx, vogdb_idx, gene_pos_idx):
                continue
            vogdb_val = row[vogdb_idx]
            if not vogdb_val:
                continue
            if not any(v in vogdb_val for v in marker_vogs):
                continue
            scaffold_val = row[scaffold_idx]
            for rep, root in rep_roots.items():
                if scaffold_val.startswith(root):
                    votu_id = votu_reps[rep]
                    if votu_id not in terminase_gene_ids:
                        gene_id = f"{scaffold_val}_{row[gene_pos_idx]}"
                        terminase_gene_ids[votu_id] = gene_id
                    break

print(f"vOTUs with a terminase large subunit hit: {len(terminase_gene_ids)} / {len(votu_reps)}")

# ---- extract protein sequences from .faa files for these gene ids ----
wanted_gene_ids = set(terminase_gene_ids.values())
gene_to_votu = {v: k for k, v in terminase_gene_ids.items()}

out_faa = "${PHYLO_DIR}/terminase_sequences.faa"
n_extracted = 0
with open("${PROTEIN_FILES_LIST}") as pf:
    protein_files = [line.strip() for line in pf if line.strip()]

with open(out_faa, "w") as fout:
    for faa_path in protein_files:
        if not faa_path.strip():
            continue
        with open(faa_path) as fin:
            write_this = False
            for line in fin:
                if line.startswith(">"):
                    header_id = line[1:].split()[0].strip()
                    write_this = header_id in wanted_gene_ids
                    if write_this:
                        votu_id = gene_to_votu.get(header_id, header_id)
                        fout.write(f">{votu_id}\\n")
                        n_extracted += 1
                elif write_this:
                    fout.write(line)

print(f"Extracted {n_extracted} terminase protein sequences -> {out_faa}")
if n_extracted < len(terminase_gene_ids) * 0.8:
    print(f"WARNING: extracted fewer sequences ({n_extracted}) than expected "
          f"({len(terminase_gene_ids)}) -- check .faa header format matches "
          f"gene_id reconstruction (scaffold_geneposition) assumption above.")
PYEOF

log "Step 12a finished. Next: 12b_align_trim_iqtree.sh"
