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
#  SCRIPT  : 02_1_standardize_species_names.sh  -  Harmonise host species names in quality reports
#  AUTHOR  : QRQ
#  INPUT   : $URGG_BINS/<bioproject>/.../quality_report.tsv (list in script)
#  OUTPUT  : same files, host names as Genus_species
# ==========================================================

export URGG_BINS="${URGG_BINS:-/path/to/basalt_binning}"

set -euo pipefail

MODE="${1:---dry-run}"

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
