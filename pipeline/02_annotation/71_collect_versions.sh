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
#  SCRIPT  : 71_collect_versions.sh  -  Report software and database versions
#  AUTHOR  : QRQ
#  INPUT   : conda environments, databases
#  OUTPUT  : version report (stdout)
# ==========================================================

export URGG_PROJECT="${URGG_PROJECT:-/path/to/metagenome_analysis}"
export URGG_HOME="${URGG_HOME:-/path/to/server_home}"

set -uo pipefail
SCR=${SCR:-$(cd "$(dirname "$0")" && pwd)}
source "${SCR}/config.sh" 2>/dev/null || echo "NOTE: config.sh not sourced; database paths from it will be blank"
hr () { printf '\n==================== %s ====================\n' "$1"; }
CONDA_BASE=$(conda info --base 2>/dev/null)
[[ -n "$CONDA_BASE" ]] && source "$CONDA_BASE/etc/profile.d/conda.sh"

hr "Tools (searched in every conda environment)"
TOOLS="run_dbcan diamond hmmsearch iqtree2 iqtree trimal FastTree fasttree mafft mmseqs amrfinder emapper.py gtdbtk checkm2 checkm prodigal dRep virsorter vibrant checkv genomad mob_recon coverm"
for env in $(conda env list 2>/dev/null | awk '!/^#/ && NF {print $1}'); do
  for t in $TOOLS; do
    p=$(conda run -n "$env" which "$t" 2>/dev/null) || continue
    [[ -z "$p" ]] && continue
    case "$t" in
      FastTree|fasttree) v=$(conda run -n "$env" "$t" -expert 2>&1 | head -1) ;;
      trimal)            v=$(conda run -n "$env" trimal --version 2>&1 | head -2 | tr '\n' ' ') ;;
      hmmsearch)         v=$(conda run -n "$env" hmmsearch -h 2>&1 | grep -m1 HMMER) ;;
      run_dbcan)         v=$( (conda run -n "$env" run_dbcan version 2>&1 || conda run -n "$env" run_dbcan --version 2>&1; conda run -n "$env" pip show dbcan 2>/dev/null | grep -i '^version') | grep -v -i error | head -2 | tr '\n' ' ') ;;
      emapper.py)        v=$(conda run -n "$env" emapper.py --version 2>&1 | tail -2 | tr '\n' ' ') ;;
      mob_recon)         v=$(conda run -n "$env" mob_recon --version 2>&1 | head -1) ;;
      *)                 v=$(conda run -n "$env" "$t" --version 2>&1 | head -2 | tr '\n' ' ') ;;
    esac
    printf '%-12s %-16s %s\n' "$env" "$t" "$v"
  done
done

hr "Databases"
echo "dbCAN database directory (CAZyDB release is in the file names):"
for d in ${DBCAN_DB:-} ${URGG_PROJECT}/databases/dbcan* ${URGG_HOME}/*/dbcan_db* ; do
  [[ -d "$d" ]] && { echo "  $d"; ls -l "$d" | grep -i -E 'cazy|dbcan|hmm|sub|fam' | awk '{print "    ", $6, $7, $8, $9}'; }
done
echo "VFDB database used by 04b (VFDB_DMND): ${VFDB_DMND:-not set}"
[[ -n "${VFDB_DMND:-}" ]] && { ls -l "${VFDB_DMND}"* 2>/dev/null; echo "  -> set A if the name contains 'setA' (core), set B if 'setB' (full)"; }
echo "VFDB sequence count (set A core ~4,000 proteins; set B full ~30,000):"
for db in ${VFDB_DMND:-} $(find ${URGG_HOME} -maxdepth 6 \( -iname 'VFDB*.dmnd' -o -iname 'VFDB_set*_pro.fas*' \) 2>/dev/null | head -6); do
  [[ -e "$db" || -e "${db}.dmnd" ]] || continue
  case "$db" in
    *.fas|*.fasta|*.faa) printf '  %s  %s sequences\n' "$db" "$(grep -c '^>' "$db")" ;;
    *.gz)                printf '  %s  %s sequences\n' "$db" "$(zcat "$db" | grep -c '^>')" ;;
    *) for env in $(conda env list 2>/dev/null | awk '!/^#/ && NF {print $1}'); do
         n=$(conda run -n "$env" diamond dbinfo -d "$db" 2>/dev/null | grep -i -E 'sequences|Database' | tr '\n' ' ') && [[ -n "$n" ]] && { echo "  $db  $n"; break; }
       done ;;
  esac
done
VFH=${URGG_PROJECT}/annotation/diamond_vfdb/diamond_vfdb_filtered.tsv
[[ -f "$VFH" ]] && echo "  distinct VFDB subject IDs hit in $VFH: $(cut -f2 "$VFH" | sort -u | wc -l) (more than ~4,000 is impossible for set A)"
echo "AMRFinderPlus database (AMRFINDER_DB): ${AMRFINDER_DB:-not set}"
[[ -n "${AMRFINDER_DB:-}" ]] && cat "${AMRFINDER_DB}/version.txt" 2>/dev/null
echo "GTDB-Tk reference data (GTDBTK_DATA_PATH): ${GTDBTK_DATA_PATH:-not set}"
[[ -n "${GTDBTK_DATA_PATH:-}" ]] && cat "${GTDBTK_DATA_PATH}/metadata/metadata.txt" 2>/dev/null | head -5

hr "R packages"
Rscript -e 'pk <- c("castor","phytools","ape","glmmTMB","MCMCglmm","lme4","data.table","ggplot2","ggforce","ggalluvial","patchwork");
  cat(R.version.string, "\n"); for (p in pk) cat(sprintf("%-11s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")))'

hr "Tree for Pagel: selection rule (select_tree_genomes.R header and filters)"
f=$(ls ${SCR}/*select_tree_genomes*.R 2>/dev/null | head -1)
if [[ -n "$f" ]]; then echo "$f"; grep -n -E '^#|filter|sample|completeness|contamination|n_per|top|head\(|slice|set.seed' "$f" | head -60
else echo "select_tree_genomes.R not found under ${SCR}"; fi
TREE=${URGG_PROJECT}/annotation/iqtree_output/tree_build_with_lineage
echo -e "\nIQ-TREE log of the 372-genome tree (version, model, bootstrap):"
ls ${TREE}/*.log ${TREE}/*.iqtree 2>/dev/null | head -3
grep -h -m3 -E 'IQ-TREE.*(version|built)|Best-fit model|Command:|-B |bootstrap' ${TREE}/*.log ${TREE}/*.iqtree 2>/dev/null | head -8
echo -e "\ntrimAl / alignment commands recorded near the tree:"
grep -h -r -m5 -E 'trimal|mafft|-gt |-automated1' ${TREE} 2>/dev/null | head -5
