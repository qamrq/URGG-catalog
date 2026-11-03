#!/usr/bin/env bash
###############################################################################
# 04b_drep_99_strain.sh
#
# Tier B: Strain-level dereplication (99% ANI) on ALL HQ genomes
# (completeness >90%, contamination <=5%, ALL sample sites)
#
# Output used for:
#   - SGVFinder2 SV reference database (all sites, all species)
#   - Strain-level codiversification analysis
#   - RUSITEC target selection
#
# Why ALL sites (not rumen-only):
#   SV analysis is being performed across all species metagenomes spanning
#   multiple GI tract sites. Restricting the reference to rumen-only would
#   cause reads from fecal/intestinal samples to fail to map or multi-map,
#   producing false negatives and false SV calls for non-rumen sites.
#
# Why HQ-only (>90% completeness):
#   Low-completeness MAGs have large missing genomic regions that SGVFinder2
#   would misinterpret as structural variants (deletions). This is a
#   technical quality requirement, not a biological filter.
#
# Run in parallel with 04a_drep_95_species.sh
# Thread allocation: 24 threads (48 go to Tier A)
###############################################################################
set -euo pipefail

# ── EDIT IF NEEDED ────────────────────────────────────────────────────────────
WORKDIR="/disk252/qamar/metagenome_analysis/drep_analysis"
THREADS=24
# ─────────────────────────────────────────────────────────────────────────────

METADATA="$WORKDIR/metadata.tsv"
GENOMEINFO="$WORKDIR/genomeInfo.csv"
GENOME_LIST="$WORKDIR/genome_list_all_HQ.txt"
OUTPUT="$WORKDIR/drep_99_strain"
LOGFILE="$WORKDIR/drep_99_strain_run.log"

# ── Build the HQ genome list (all sites, >90% completeness) if not present ───
# metadata.tsv columns:
#   1=genome  2=bioproject  3=host  4=species  5=sample_site
#   6=quality_tier  7=completeness  8=contamination  9=N50
if [[ ! -f "$GENOME_LIST" ]]; then
    echo ">> Building all-site HQ genome list from metadata.tsv..."
    awk -F'\t' 'NR>1 && $6=="HQ" {print $1}' "$METADATA" \
        | sed "s|^|$WORKDIR/genomes/|" > "$GENOME_LIST"
    n=$(wc -l < "$GENOME_LIST")
    echo "   Written $n HQ genomes (all sites) to $GENOME_LIST"
    echo ""

    # show site distribution of what's going in
    echo "   Site distribution of HQ genomes:"
    awk -F'\t' 'NR>1 && $6=="HQ" {print $5}' "$METADATA" \
        | sort | uniq -c | sort -rn \
        | awk '{printf "     %-25s %d genomes\n", $2, $1}'
    echo ""
fi

echo "========================================================"
echo "Tier B: Strain-level dereplication (99% ANI)"
echo "Input genomes: $(wc -l < "$GENOME_LIST") (all sites, HQ only)"
echo "Threads: $THREADS"
echo "Output: $OUTPUT"
echo "Started: $(date)"
echo "========================================================"

# validate inputs
for f in "$GENOME_LIST" "$GENOMEINFO"; do
    [[ -f "$f" ]] || { echo "ERROR: required file not found: $f" >&2; exit 1; }
done

/usr/bin/time -v dRep dereplicate "$OUTPUT" \
    -p $THREADS \
    --genomeInfo "$GENOMEINFO" \
    --ignoreGenomeQuality \
    -pa 0.95 -sa 0.99 -nc 0.30 \
    -cm larger \
    --S_algorithm skani \
    --multiround_primary_clustering \
    --primary_chunksize 2000 \
    -g "$GENOME_LIST" \
    2>&1 | tee "$LOGFILE"

echo ""
echo "========================================================"
echo "Tier B complete: $(date)"
echo "Strain representatives: $OUTPUT/dereplicated_genomes/"
echo "Rep count: $(ls $OUTPUT/dereplicated_genomes/ | wc -l)"
echo ""
echo "Key output tables (preserve these):"
echo "  $OUTPUT/data_tables/Cdb.csv   cluster assignments (strain)"
echo "  $OUTPUT/data_tables/Wdb.csv   winning representatives"
echo "  $OUTPUT/data_tables/Sdb.csv   genome scores"
echo "  $OUTPUT/data_tables/Ndb.csv   pairwise ANI cache"
echo ""
echo "Next step: build SGVFinder2 reference database from these representatives:"
echo "  mkdir -p $WORKDIR/sv_reference_db"
echo "  cp $OUTPUT/dereplicated_genomes/*.fa $WORKDIR/sv_reference_db/"
echo "========================================================"