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
#  SCRIPT  : 03_prep_basalt_for_drep.sh  -  Collect QC-passing bins, rename and build dRep inputs
#  AUTHOR  : QRQ
#  INPUT   : $URGG_BINS quality reports and bins
#  OUTPUT  : drep_analysis/genomes/
#            drep_analysis/metadata.tsv, genomeInfo.csv, genome lists
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_BINS="${URGG_BINS:-/path/to/basalt_binning}"

set -euo pipefail

WORKDIR="${URGG_PROJECT}/drep_analysis"

PROJECTS=(
"${URGG_BINS}/CNP0001390_MAGs/Yak_MAGs_QS50_n28125/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Sheep/0.MAGs/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Cattle_MAGs_QS50_n10630/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Antelope/0.MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB12797_DONE/Bestbinset/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
)

# Quality thresholds (from CheckM2 values in TSV)
CATALOG_MIN_COMP=70
HQ_MIN_COMP=90
MAX_CONT=5
 
mkdir -p "$WORKDIR/genomes"
 
if [[ ! -f "$WORKDIR/genome_list.txt" ]];      then : > "$WORKDIR/genome_list.txt"; fi
if [[ ! -f "$WORKDIR/genome_list_rumen.txt" ]]; then : > "$WORKDIR/genome_list_rumen.txt"; fi
if [[ ! -f "$WORKDIR/metadata.tsv" ]]; then
    echo -e "genome\tbioproject\thost\tspecies\tsample_site\tquality_tier\tcompleteness\tcontamination\tN50" \
        > "$WORKDIR/metadata.tsv"
fi
if [[ ! -f "$WORKDIR/genomeInfo.csv" ]]; then
    echo "genome,completeness,contamination" > "$WORKDIR/genomeInfo.csv"
fi
 
VALID_SITES="Rumen Cecum Colon Duodenum Jejunum Ileum Rumen_solid Rectum Reticulum Omasum Abomasum Small_intestine Large_intestine Feces Multi_site Unknown"
 
is_valid_site() {
    local s="$1"
    for v in $VALID_SITES; do
        [[ "$s" == "$v" ]] && return 0
    done
    return 1
}
 
total=0; kept=0; hq=0; rumen_count=0; dup_skip=0
 
for qtsv in "${PROJECTS[@]}"; do
 
    if [[ ! -f "$qtsv" ]]; then
        echo "ERROR: TSV not found: $qtsv" >&2; exit 1
    fi
 
    bindir=$(dirname "$qtsv")
    echo ">> Processing: $qtsv"
 
    delim=$'\t'
    header_raw=$(head -1 "$qtsv")
    n_commas=$( (grep -o ',' <<< "$header_raw" || true) | wc -l )
    n_tabs=$( (grep -o $'\t' <<< "$header_raw" || true) | wc -l )
    if [[ "$n_commas" -gt "$n_tabs" ]]; then
        delim=","
    fi
    echo "   detected delimiter: $( [[ "$delim" == "," ]] && echo 'comma' || echo 'tab' )"
 
    # Read header to get column indices (robust to order)
    header="$header_raw"
    IFS="$delim" read -ra cols <<< "$header"
 
    col_binid=-1; col_comp=-1; col_cont=-1; col_n50=-1
    col_site=-1; col_bp=-1; col_host=-1; col_species=-1
 
    for i in "${!cols[@]}"; do
        cname="${cols[$i]//$'\r'/}"
        case "$cname" in
            Bin_ID)                    col_binid=$i  ;;
            Completeness|Completenes)  col_comp=$i   ;;
            Contamination)             col_cont=$i   ;;
            N50)                       col_n50=$i    ;;
            sample_site|Sample_site)   col_site=$i   ;;
            bioproject)                col_bp=$i     ;;
            host)                      col_host=$i   ;;
            species)                   col_species=$i ;;
        esac
    done
 
    # Required columns (species is optional — older files may lack it)
    for cname in col_binid col_comp col_cont col_n50 col_site col_bp col_host; do
        if [[ "${!cname}" == "-1" ]]; then
            echo "ERROR: Missing column '${cname#col_}' in $qtsv" >&2; exit 1
        fi
    done
    if [[ "$col_species" == "-1" ]]; then
        echo "  NOTE: no 'species' column in this file — will record species as 'unknown'" >&2
    fi
 
    while IFS="$delim" read -ra fields; do
        for k in "${!fields[@]}"; do fields[$k]="${fields[$k]//$'\r'/}"; done
 
        [[ "${fields[$col_binid]:-}" == "Bin_ID" ]] && continue
        [[ -z "${fields[$col_binid]:-}" ]] && continue
 
        bin_id="${fields[$col_binid]}"
        comp="${fields[$col_comp]}"
        cont="${fields[$col_cont]}"
        n50="${fields[$col_n50]}"
        site="${fields[$col_site]}"
        bp="${fields[$col_bp]}"
        host="${fields[$col_host]}"
        species="unknown"
        [[ "$col_species" != "-1" ]] && species="${fields[$col_species]:-unknown}"
        total=$((total+1))
 
        # validate site (case-sensitive exact match against full vocabulary)
        if ! is_valid_site "$site"; then
            echo "  WARN: unrecognised site '$site' for bin $bin_id — using 'Unknown'" >&2
            site="Unknown"
        fi
 
        # quality gate
        pass=$(awk -v c="$comp" -v x="$cont" -v mc=$CATALOG_MIN_COMP -v xc=$MAX_CONT \
            'BEGIN{ print (c+0>=mc && x+0<=xc) ? 1 : 0 }')
        [[ "$pass" == "0" ]] && continue
 
        # quality tier
        tier=$(awk -v c="$comp" -v hq=$HQ_MIN_COMP \
            'BEGIN{ print (c+0>hq) ? "HQ" : "MQ" }')
 
        # locate bin file in same folder as TSV
        src=""
        for ext in "" ".fa" ".fasta" ".fna"; do
            cand="$bindir/${bin_id}${ext}"
            [[ -f "$cand" ]] && { src="$cand"; break; }
        done
        if [[ -z "$src" ]]; then
            echo "  WARN: bin file not found for '$bin_id' in $bindir" >&2
            continue
        fi
 
        # build origin-encoding filename
        safe_id=$(echo "$bin_id" | sed 's/[^A-Za-z0-9._-]/_/g')
        newname="${bp}__${host}__${site}__${safe_id}.fa"
 
        # skip duplicates
        if [[ -f "$WORKDIR/genomes/$newname" ]]; then
            echo "  WARN: duplicate filename '$newname' — skipping" >&2
            dup_skip=$((dup_skip+1))
            continue
        fi
 
        cp "$src" "$WORKDIR/genomes/$newname"
 
        echo "$WORKDIR/genomes/$newname" >> "$WORKDIR/genome_list.txt"
        echo -e "${newname}\t${bp}\t${host}\t${species}\t${site}\t${tier}\t${comp}\t${cont}\t${n50}" \
            >> "$WORKDIR/metadata.tsv"
        echo "${newname},${comp},${cont}" >> "$WORKDIR/genomeInfo.csv"
 
        if [[ "$site" == "Rumen" ]]; then
            echo "$WORKDIR/genomes/$newname" >> "$WORKDIR/genome_list_rumen.txt"
            rumen_count=$((rumen_count+1))
        fi
 
        kept=$((kept+1))
        [[ "$tier" == "HQ" ]] && hq=$((hq+1))
 
    done < "$qtsv"
 
done
 
echo ""
echo "─────────────────────────────────────────────"
echo "Scanned:          $total bins"
echo "Passed (>=${CATALOG_MIN_COMP}% comp, <=${MAX_CONT}% cont): $kept"
echo "  of which HQ (>${HQ_MIN_COMP}% comp): $hq"
echo "  of which Rumen:                $rumen_count"
echo "  duplicates skipped:            $dup_skip"
echo ""
echo "Output files:"
echo "  All sites:   $WORKDIR/genome_list.txt        ($(wc -l < "$WORKDIR/genome_list.txt") genomes)"
echo "  Rumen only:  $WORKDIR/genome_list_rumen.txt  ($(wc -l < "$WORKDIR/genome_list_rumen.txt") genomes)"
echo "  Metadata:    $WORKDIR/metadata.tsv   (columns: genome,bioproject,host,species,sample_site,quality_tier,completeness,contamination,N50)"
echo "  genomeInfo:  $WORKDIR/genomeInfo.csv"
echo ""
echo "To build strain catalog genome list (Rumen HQ only):"
echo "  awk -F'\t' 'NR==1||(\$6==\"HQ\"&&\$5==\"Rumen\"){print \$1}' $WORKDIR/metadata.tsv \\"
echo "    | tail -n+2 | sed \"s|^|$WORKDIR/genomes/|\" > $WORKDIR/genome_list_rumen_HQ.txt"
