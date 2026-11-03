# Copy this file to config/local_paths.R (git-ignored) and set the two paths for
# your machine. Environment variables URGG_PROJECT / URGG_FIGURES override it.
#
# URGG_PROJECT: folder that holds the project layout used on the analysis server
#   annotation/results_phase3/...      (pipeline outputs: s49_, s63_, s66_, s80_ ...)
#   annotation/gtdbtk_output_v232/...  (GTDB-Tk summaries)
#   website_release/web_data/...       (web tables, e.g. web_votus.tsv)
# On a PC, copy only the files a script needs into the same relative sub-folders.

# Analysis server
URGG_PROJECT <- "/disk252/qamar/metagenome_analysis"
URGG_FIGURES <- "/disk252/qamar/metagenome_analysis/figures_out"

# Windows PC (example)
# URGG_PROJECT <- "D:/PD-ZJU/MAG-SAG/urgg_project"
# URGG_FIGURES <- "D:/PD-ZJU/MAG-SAG/Figures"
