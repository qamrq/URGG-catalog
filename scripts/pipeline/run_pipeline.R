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
#  SCRIPT  : run_pipeline.R  -  Run construction steps 01-15 in order
#  AUTHOR  : QRQ
#  INPUT   : steps 01-15 in scripts/pipeline/ and config.R
#  OUTPUT  : the outputs of each step (see the banner of each script)
# ==========================================================
#
# Usage, from the repository root:
#   Rscript scripts/pipeline/run_pipeline.R            # all steps
#   Rscript scripts/pipeline/run_pipeline.R 06 09      # steps 06 to 09
#   Rscript scripts/pipeline/run_pipeline.R 03         # one step
# Step 02 is run in its Part B form (mapping). Build the HoCoRT index once with
# build_hocort_index() in 02_hocort_host_removal.R before the first run.

PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", "scripts/pipeline")
Sys.setenv(PIPELINE_DIR = PIPELINE_DIR)
source(file.path(PIPELINE_DIR, "config.R"))

steps <- data.frame(
  id   = sprintf("%02d", 1:15),
  file = c("01_fastp_qc.R", "02_hocort_host_removal.R", "03_megahit_assembly.R",
           "04_basalt_binning.R", "05_checkm2_quality.R", "06_prepare_drep_inputs.R",
           "07_drep_stage1_per_host.R", "08_drep_stage2_species.R", "09_drep_strain_99.R",
           "10_gtdbtk_classify.R", "11_gene_prediction_prodigal.R",
           "12_gene_catalog_mmseqs2.R", "13_eggnog_annotation.R", "14_dbcan_cazyme.R",
           "15_diamond_nr_unassigned.R"),
  fun  = c("step_fastp_qc", "step_hocort_map", "step_megahit", "step_basalt",
           "step_checkm2", "step_prepare_drep", "step_drep_stage1", "step_drep_stage2",
           "step_drep_strain", "step_gtdbtk", "step_prodigal", "step_gene_catalog",
           "step_eggnog", "step_dbcan", "step_diamond_nr"),
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
from <- if (length(args) >= 1) args[1] else "01"
to   <- if (length(args) >= 2) args[2] else if (length(args) == 1) args[1] else "15"

selected <- steps[steps$id >= from & steps$id <= to, ]
if (nrow(selected) == 0) stop("No steps between ", from, " and ", to)

for (i in seq_len(nrow(selected))) {
  s <- selected[i, ]
  log_msg("===== Step ", s$id, ": ", s$file, " =====")
  source(file.path(PIPELINE_DIR, s$file))
  do.call(s$fun, list())
}
log_msg("Done: steps ", from, " to ", to)
