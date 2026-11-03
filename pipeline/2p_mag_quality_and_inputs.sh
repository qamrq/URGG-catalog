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
#  SCRIPT  : 2p_mag_quality_and_inputs.sh  -  MAG quality (CheckM2), host-name harmonisation and dereplication inputs
#  AUTHOR  : QRQ
#  INPUT   : $URGG_BINS/<bioproject>/<bin folders>
#  OUTPUT  : checkm2_summary.tsv per folder
#            drep_analysis/genomes/, metadata.tsv, genomeInfo.csv, genome lists
#  USAGE   : bash 2p_mag_quality_and_inputs.sh [checkm2 names inputs]   (no argument: all steps in this order)
# ==========================================================

export URGG_BINS="${URGG_BINS:-/path/to/basalt_binning}"
export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"

set -euo pipefail

step_checkm2() (
set -euo pipefail

THREADS=32
EXTENSION="fa"

MAG_DIRS=(
    "${URGG_BINS}/PRJCA008889_DONE/Gazelle/BestBinset_outlier_refined_filtrated/"
)

# confirm checkm2 is available
if ! command -v checkm2 &>/dev/null; then
    echo "ERROR: checkm2 not found. Add it to PATH first." >&2
    exit 1
fi

echo "========================================================"
echo "Batch CheckM2 — ${#MAG_DIRS[@]} folders"
echo "Threads: $THREADS | Extension: .$EXTENSION"
echo "Started: $(date)"
echo "========================================================"
echo ""

PASS=0; FAIL=0; SKIP=0

for magdir in "${MAG_DIRS[@]}"; do

    echo "── $(date '+%H:%M:%S') | $(basename $(dirname $magdir))/$(basename $magdir)"

    # check directory exists
    if [[ ! -d "$magdir" ]]; then
        echo "  ERROR: directory not found — skipping" >&2
        FAIL=$((FAIL+1)); continue
    fi

    # count MAG files
    n_mags=$(find "$magdir" -maxdepth 1 -name "*.${EXTENSION}" | wc -l)
    if [[ "$n_mags" -eq 0 ]]; then
        echo "  WARN: no .${EXTENSION} files found — skipping" >&2
        FAIL=$((FAIL+1)); continue
    fi
    echo "  MAGs found: $n_mags"

    checkm2_out="$magdir/checkm2_output"
    summary_tsv="$magdir/checkm2_summary.tsv"

    # skip if already successfully completed
    if [[ -f "$checkm2_out/quality_report.tsv" && -f "$summary_tsv" ]]; then
        n_done=$(wc -l < "$summary_tsv")
        echo "  Already done ($((n_done-1)) MAGs in summary) — skipping"
        echo "  (delete checkm2_output/ to force rerun)"
        SKIP=$((SKIP+1)); continue
    fi

    # run CheckM2
    echo "  Running CheckM2..."
    checkm2 predict \
        --input "$magdir" \
        --output-directory "$checkm2_out" \
        --threads "$THREADS" \
        --extension "$EXTENSION" \
        --force \
        2>&1 | grep -E "INFO|WARNING|ERROR|complete|Finished" || true

    if [[ ! -f "$checkm2_out/quality_report.tsv" ]]; then
        echo "  ERROR: CheckM2 failed — quality_report.tsv not produced" >&2
        FAIL=$((FAIL+1)); continue
    fi

    # format to clean summary TSV
    python3 - << PYEOF
import csv, sys

infile  = "$checkm2_out/quality_report.tsv"
outfile = "$summary_tsv"

with open(infile, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)

# verify required columns
required = {"Name", "Genome_Size", "Completeness", "Contamination", "Contig_N50"}
missing  = required - set(reader.fieldnames)
if missing:
    sys.exit(f"ERROR: CheckM2 report missing columns: {missing}")

with open(outfile, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Bin_ID", "Genome_size", "Completeness", "Contamination", "N50"])
    for row in rows:
        writer.writerow([
            row["Name"],
            row["Genome_Size"],
            row["Completeness"],
            row["Contamination"],
            row["Contig_N50"],
        ])

# quality distribution
hq = mq = fail_q = 0
for row in rows:
    c = float(row["Completeness"])
    x = float(row["Contamination"])
    if   c > 90  and x <= 5: hq     += 1
    elif c >= 70 and x <= 5: mq     += 1
    else:                    fail_q += 1

print(f"  Total: {len(rows)}  |  HQ (>90/<5): {hq}  |  MQ (>=70/<5): {mq}  |  Fail: {fail_q}")
print(f"  Summary saved: {outfile}")
PYEOF

    PASS=$((PASS+1))
    echo ""

done

echo "========================================================"
echo "Finished: $(date)"
echo "Completed: $PASS  |  Skipped (already done): $SKIP  |  Failed: $FAIL"
echo ""
echo "Output files per folder:"
echo "  checkm2_output/quality_report.tsv   (full CheckM2 output)"
echo "  checkm2_summary.tsv                 (Bin_ID|Genome_size|Completeness|Contamination|N50)"
echo ""
echo "Next steps for each folder:"
echo "  1. Add columns: sample_site, bioproject, host to checkm2_summary.tsv"
echo "  2. Rename to:   Best_binset_quality_report_final.tsv"
echo "  3. Add path to: prep_basalt_for_drep.sh PROJECTS array"
echo ""
echo "For multi-site projects (e.g. PRJNA656389), run:"
echo "  python3 add_sample_site_to_quality_report.py \\"
echo "      --quality  <folder>/checkm2_summary.tsv \\"
echo "      --metadata <project_metadata.csv> \\"
echo "      --output   <folder>/Best_binset_quality_report_final.tsv"
)

step_names() (
set -euo pipefail

MODE="${NAME_MODE:---dry-run}"

PROJECT_FILES=(
"${URGG_BINS}/PRJNA214227/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA486341/Goat/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA486341/Camel_drom/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA987743_MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJNA983076_DONE/Giraffa_camelopardalis/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA983076_DONE/Camelus_bactrianus/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA972320/Goat/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA657455/Yak/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"${URGG_BINS}/PRJNA657455/Water_deer/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"${URGG_BINS}/PRJNA657455/Sheep/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"${URGG_BINS}/PRJNA656389_MAGs/MAGs_4960/checkm2_summary_with_sample_site.tsv"
"${URGG_BINS}/PRJNA603266/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA60251/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA526070_DONE/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA512570/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA485488/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA448333/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJNA362629/BestBinset_outlier_refined_filtrated/quality_report_with_sample_site.tsv"
"${URGG_BINS}/PRJNA301235/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA298982/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJNA260105/BestBinset/BestBinSet_quality_report.tsv"
"${URGG_BINS}/PRJNA1069655/BestBinset/Best_binset_quality_report_with_site.tsv"
"${URGG_BINS}/PRJNA419351_DONE/PRJNA419351_DONE_Disk252/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
"${URGG_BINS}/PRJNA1034944_DONE/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJEB8939/BestBinset_outlier_refined_filtrated_retrieved/BestBinset_outlier_refined_filtrated_retrieved_quality_report.tsv"
"${URGG_BINS}/PRJEB39057_MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB38078/Giraffe/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJEB38078/Gaur/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJEB38078/Camelus_dromedarius/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJEB38078/Camelus_bact/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJEB34458_MAGs/Sheep/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB34458_MAGs/Reindeer/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB34458_MAGs/Red_deer/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB34458_MAGs/Cattle/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB21624_PRJEB31266_MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB21624_MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB14312/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJEB10338/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJCA008889_DONE/Sheep/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJCA008889_DONE/Goat/BestBinset/Best_binset_quality_report.tsv"
"${URGG_BINS}/PRJCA008889_DONE/Antelope/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"${URGG_BINS}/PRJCA008889_DONE/Gazelle/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
"${URGG_BINS}/Gp0348364_Gp0348367_JGI_dairy_cattle_USA/Gp0348364_Gp0348367_JGI_dairy_cattle_USA.csv"
"${URGG_BINS}/Gp0306217_JGI_tropical_cattle_AUS/Gp0306217_JGI_tropical_cattle_AUS_Best_binset_quality_report.csv"
"${URGG_BINS}/Gp0294691_Gp0294694_JGI_beef_cattle_CAN/Gp0294691_Gp0294694_JGI_beef_cattle_CAN.csv"
"${URGG_BINS}/Gp0290826_JGI_sheep_NZ/Gp0290826_JGI_sheep_NZ.csv"
"${URGG_BINS}/Gp0138594_JGI_beef_cattle_BRA/Gp0138594_JGI_beef_cattle_BRA.csv"
"${URGG_BINS}/CNP0001390_MAGs/Yak_MAGs_QS50_n28125/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Sheep/0.MAGs/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Cattle_MAGs_QS50_n10630/checkm2_summary.tsv"
"${URGG_BINS}/CNP0001390_MAGs/Antelope/0.MAGs/checkm2_summary.tsv"
"${URGG_BINS}/PRJEB12797_DONE/Bestbinset/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
)

if [[ "$MODE" != "--dry-run" && "$MODE" != "--apply" ]]; then
    echo "Usage: $0 [--dry-run|--apply]" >&2
    exit 1
fi

echo "Mode: $MODE"
echo ""

for f in "${PROJECT_FILES[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "SKIP (not found): $f" >&2
        continue
    fi

    # detect delimiter
    case "$f" in
        *.csv) delim="," ;;
        *.tsv) delim=$'\t' ;;
        *)
            n_commas=$(head -1 "$f" | grep -o "," | wc -l)
            n_tabs=$(head -1 "$f" | grep -o $'\t' | wc -l)
            if [[ "$n_commas" -gt "$n_tabs" ]]; then delim=","; else delim=$'\t'; fi
            ;;
    esac

    # find column indices for species, sample_site, host, bioproject
    header=$(head -1 "$f")
    IFS="$delim" read -ra cols <<< "$header"
    col_species=-1; col_site=-1; col_host=-1; col_bp=-1
    for i in "${!cols[@]}"; do
        cname="${cols[$i]//$'\r'/}"
        case "$cname" in
            species) col_species=$i ;;
            sample_site|Sample_site) col_site=$i ;;
            host) col_host=$i ;;
            bioproject) col_bp=$i ;;
        esac
    done

    n_changes=0
    tmpfile=$(mktemp)

    while IFS="$delim" read -ra fields; do
        for k in "${!fields[@]}"; do
            orig="${fields[$k]}"
            val="${orig//$'\r'/}"
            val="${val#"${val%%[![:space:]]*}"}"
            val="${val%"${val##*[![:space:]]}"}"

            if [[ $k -eq $col_species && "$val" == *" "* ]]; then
                val="${val// /_}"
            fi

            if [[ "$val" != "$orig" ]]; then
                n_changes=$((n_changes+1))
            fi
            fields[$k]="$val"
        done
        # rejoin with original delimiter
        line=""
        for k in "${!fields[@]}"; do
            if [[ $k -eq 0 ]]; then line="${fields[$k]}"; else line="${line}${delim}${fields[$k]}"; fi
        done
        echo "$line" >> "$tmpfile"
    done < "$f"

    echo "── $(basename $f)"
    echo "   delim='$delim'  species_col=$col_species  site_col=$col_site  host_col=$col_host  bp_col=$col_bp"
    echo "   cells changed (whitespace/species-underscore fixes): $n_changes"

    if [[ "$MODE" == "--apply" && "$n_changes" -gt 0 ]]; then
        timestamp=$(date +%Y%m%d_%H%M%S)
        base="${f%.*}"
        ext="${f##*.}"
        backup_path="${base}_backup_before_standardize_${timestamp}.${ext}"
        cp "$f" "$backup_path"
        mv "$tmpfile" "$f"
        echo "   APPLIED. Backup saved: $backup_path"
    else
        rm -f "$tmpfile"
        [[ "$n_changes" -gt 0 ]] && echo "   (dry-run — no changes written; rerun with --apply)"
    fi
    echo ""

done

echo "Done."
[[ "$MODE" == "--dry-run" ]] && echo "This was a DRY RUN. Rerun with --apply to write changes."
)

step_inputs() (
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
)

for s in ${@:-checkm2 names inputs}; do
  echo "==> 2p_mag_quality_and_inputs.sh: $s"
  "step_$s"
done
