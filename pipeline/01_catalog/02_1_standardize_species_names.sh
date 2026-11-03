#!/usr/bin/env bash
###############################################################################
# standardize_species_names.sh
#
# Fixes inconsistent species naming across quality-report files:
#   "Bos taurus"  -> "Bos_taurus"
#   "Ovis aries"  -> "Ovis_aries"
#   (any "Genus species" with a literal space -> "Genus_species")
#
# Also trims trailing/leading whitespace from sample_site, host, bioproject,
# and species columns generally (catches the "Rumen " issue you found).
#
# IMPORTANT: this script EDITS FILES IN PLACE (with a .bak backup created
# alongside each file first). Review the dry-run output before confirming.
#
# USAGE:
#   bash standardize_species_names.sh --dry-run     # preview changes only
#   bash standardize_species_names.sh --apply       # actually modify files
###############################################################################
set -euo pipefail

MODE="${1:---dry-run}"

# ── EDIT THIS: same file list as your scan script ─────────────────────────────
PROJECT_FILES=(
"/disk252/qamar/basalt_binning/PRJNA214227/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA486341/Goat/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA486341/Camel_drom/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA987743_MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJNA983076_DONE/Giraffa_camelopardalis/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA983076_DONE/Camelus_bactrianus/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA972320/Goat/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA657455/Yak/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA657455/Water_deer/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA657455/Sheep/BestBinset/Best_binset_quality_report_with_sample_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA656389_MAGs/MAGs_4960/checkm2_summary_with_sample_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA603266/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA60251/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA526070_DONE/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA512570/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA485488/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA448333/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA362629/BestBinset_outlier_refined_filtrated/quality_report_with_sample_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA301235/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA298982/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA260105/BestBinset/BestBinSet_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJNA1069655/BestBinset/Best_binset_quality_report_with_site.tsv"
"/disk252/qamar/basalt_binning/PRJNA419351_DONE/PRJNA419351_DONE_Disk252/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJNA1034944_DONE/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB8939/BestBinset_outlier_refined_filtrated_retrieved/BestBinset_outlier_refined_filtrated_retrieved_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB39057_MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB38078/Giraffe/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB38078/Gaur/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB38078/Camelus_dromedarius/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB38078/Camelus_bact/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB34458_MAGs/Sheep/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB34458_MAGs/Reindeer/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB34458_MAGs/Red_deer/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB34458_MAGs/Cattle/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB21624_PRJEB31266_MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB21624_MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB14312/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJEB10338/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJCA008889_DONE/Sheep/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJCA008889_DONE/Goat/BestBinset/Best_binset_quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJCA008889_DONE/Antelope/BestBinset_outlier_refined_filtrated/quality_report.tsv"
"/disk252/qamar/basalt_binning/PRJCA008889_DONE/Gazelle/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/Gp0348364_Gp0348367_JGI_dairy_cattle_USA/Gp0348364_Gp0348367_JGI_dairy_cattle_USA.csv"
"/disk252/qamar/basalt_binning/Gp0306217_JGI_tropical_cattle_AUS/Gp0306217_JGI_tropical_cattle_AUS_Best_binset_quality_report.csv"
"/disk252/qamar/basalt_binning/Gp0294691_Gp0294694_JGI_beef_cattle_CAN/Gp0294691_Gp0294694_JGI_beef_cattle_CAN.csv"
"/disk252/qamar/basalt_binning/Gp0290826_JGI_sheep_NZ/Gp0290826_JGI_sheep_NZ.csv"
"/disk252/qamar/basalt_binning/Gp0138594_JGI_beef_cattle_BRA/Gp0138594_JGI_beef_cattle_BRA.csv"
"/disk252/qamar/basalt_binning/CNP0001390_MAGs/Yak_MAGs_QS50_n28125/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/CNP0001390_MAGs/Sheep/0.MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/CNP0001390_MAGs/Cattle_MAGs_QS50_n10630/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/CNP0001390_MAGs/Antelope/0.MAGs/checkm2_summary.tsv"
"/disk252/qamar/basalt_binning/PRJEB12797_DONE/Bestbinset/BestBinset_outlier_refined_filtrated/checkm2_summary.tsv"
)
# ─────────────────────────────────────────────────────────────────────────────

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
            val="${orig//$'\r'/}"          # strip CR
            val="${val#"${val%%[![:space:]]*}"}"   # trim leading space
            val="${val%"${val##*[![:space:]]}"}"   # trim trailing space

            # only convert "Genus species" -> "Genus_species" in the species column
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