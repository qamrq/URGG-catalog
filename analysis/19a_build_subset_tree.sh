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
#  SCRIPT  : 19a_build_subset_tree.sh  -  Subset tree from GTDB-Tk marker alignments (IQ-TREE)
#  AUTHOR  : QRQ
#  INPUT   : genome list, GTDB-Tk align folder, output folder (arguments)
#  OUTPUT  : bac120_tree.treefile, ar53_tree.treefile
# ==========================================================

set -euo pipefail

GENOME_LIST="$1"
ALIGN_DIR="$2"
OUTDIR="$3"

mkdir -p "$OUTDIR"

BAC_MSA_GZ="$ALIGN_DIR/gtdbtk.bac120.user_msa.fasta.gz"
AR_MSA_GZ="$ALIGN_DIR/gtdbtk.ar53.user_msa.fasta.gz"

for f in "$GENOME_LIST" "$BAC_MSA_GZ" "$AR_MSA_GZ"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

n_target=$(wc -l < "$GENOME_LIST")
echo "Target genome list: $n_target genomes"

for domain in bac120 ar53; do
    gz="$ALIGN_DIR/gtdbtk.${domain}.user_msa.fasta.gz"
    plain="$OUTDIR/gtdbtk.${domain}.user_msa.fasta"
    if [[ ! -f "$plain" ]]; then
        echo ">> Decompressing $domain alignment..."
        gunzip -c "$gz" > "$plain"
    fi
done

for domain in bac120 ar53; do
    echo ""
    echo "======================================================"
    echo "Domain: $domain"
    echo "======================================================"

    plain_msa="$OUTDIR/gtdbtk.${domain}.user_msa.fasta"
    subset_msa="$OUTDIR/${domain}_subset.fasta"

    # extract headers present in this domain's alignment
    grep "^>" "$plain_msa" | sed 's/^>//' | cut -d' ' -f1 > "$OUTDIR/${domain}_all_headers.txt"

    # intersect with target genome list
    comm -12 <(sort "$GENOME_LIST") <(sort "$OUTDIR/${domain}_all_headers.txt") \
        > "$OUTDIR/${domain}_matched_genomes.txt"

    n_matched=$(wc -l < "$OUTDIR/${domain}_matched_genomes.txt")
    echo "  Genomes from target list found in $domain alignment: $n_matched"

    if [[ "$n_matched" -eq 0 ]]; then
        echo "  No genomes for this domain -- skipping tree build."
        continue
    fi

    if command -v seqkit &>/dev/null; then
        seqkit grep -f "$OUTDIR/${domain}_matched_genomes.txt" "$plain_msa" > "$subset_msa"
    else
        echo "  ERROR: seqkit not found on PATH" >&2
        exit 1
    fi

    n_extracted=$(grep -c "^>" "$subset_msa" || true)
    echo "  Sequences extracted: $n_extracted"

    if [[ "$n_extracted" -lt 4 ]]; then
        echo "  Fewer than 4 sequences -- not enough to build a meaningful tree, skipping."
        continue
    fi

    trimmed="$OUTDIR/${domain}_trimmed.fasta"
    echo "  Running trimAl..."
    trimal -in "$subset_msa" -out "$trimmed" -automated1

    echo "  Running IQ-TREE (-m MFP -bb 1000 -redo)..."
    iqtree -s "$trimmed" -m MFP -bb 1000 -redo -nt 48 -pre "$OUTDIR/${domain}_tree"

    echo "  Tree written: $OUTDIR/${domain}_tree.treefile"
done

echo ""
echo "======================================================"
echo "Done."
echo "======================================================"
echo "Bacterial tree: $OUTDIR/bac120_tree.treefile"
echo "Archaeal tree:  $OUTDIR/ar53_tree.treefile"
echo ""
echo "NEXT: use tree_ring_annotation.csv (already built) as the ring-annotation"
echo "input for the ggtree figure script, joined by genome ID against each"
echo "tree's tip labels."
