#!/usr/bin/env bash
###############################################################################
# 12a_extract_terminase_marker.sh — Phase 5: marker gene discovery + extraction
# for real phylogenetic tree construction
#
# Uses terminase large subunit (TerL) as the marker gene -- the standard
# single-copy phylogenetic marker for Caudoviricetes (93.8% of this
# catalogue), analogous to GTDB-Tk's bacterial/archaeal marker set. TerL is
# universally conserved across tailed phages and is the field-standard
# choice for phage phylogenomics (same principle as GTDB-Tk's approach:
# identify a conserved single-copy marker, then align/tree on that alone
# rather than whole genomes, which don't share synteny across viral taxa).
#
# Step 1: discover real terminase-large-subunit VOG accessions from VOGDB's
#         own functional descriptions (same verified technique as
#         05_prophage_flagging.sh's marker discovery -- do not hardcode
#         accessions).
# Step 2: find genes hitting those VOGs on each vOTU representative (using
#         the CONFIRMED WORKING "scaffold" column + prefix-match join
#         technique from 05_prophage_flagging.sh -- do not reuse the
#         "fasta" column, which was proven unreliable in that file).
# Step 3: extract the corresponding PROTEIN sequence for that gene from
#         DRAM-v's predicted-protein output.
#
# *** VERIFY BEFORE RUNNING ***
# This assumes DRAM-v's protein FASTA (typically genes.faa, from its
# internal Prodigal call) exists per batch. Confirm with:
#   find /disk252/qamar/metagenome_analysis/phase_5_virome/05_dramv -iname "*.faa" | head -5
# and correct PROTEIN_FAA_PATTERN below if the real location/naming differs.
#
# Output: RESULTS_DIR/terminase_marker_vogs.tsv
#         RESULTS_DIR/terminase_sequences.faa   (one TerL protein per vOTU
#                                                 that has one -- NOT all
#                                                 vOTUs will, especially
#                                                 non-Caudoviricetes ones)
###############################################################################
source "/disk252/qamar/metagenome_analysis/phase_5_virome/scripts/config.sh"

RESULTS_DIR="${PHASE5_DIR}/08_results"
PHYLO_DIR="${PHASE5_DIR}/12_phylogeny"
mkdir -p "${PHYLO_DIR}"

VOG_ANNOT_GZ="/disk252/qamar/databases/dramv_config/database_download/vog.annotations.tsv.gz"
if [[ ! -f "${VOG_ANNOT_GZ}" ]]; then
    log "ERROR: ${VOG_ANNOT_GZ} not found. Adjust path if your install differs."
    exit 1
fi

###############################################################################
# Step 1: discover real terminase-large-subunit VOGs (not hardcoded)
###############################################################################
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

###############################################################################
# Step 2/3: find genes hitting these VOGs + extract protein sequences
# *** VERIFY: confirm this glob finds real protein FASTA files before trusting output ***
###############################################################################
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
# BUG FIX: PROTEIN_FILES is multi-line (one path per line from find). Embedding
# it directly inline into a Python string literal inside the heredoc breaks
# (unterminated string, since Python literals can't span raw newlines like
# that). Write it to a real file on disk instead and have Python read that.
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