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
#  SCRIPT  : run_pipeline.R  -  Construction pipeline, steps 01-15 in order
#  AUTHOR  : QRQ
#  INPUT   : steps 01-15 in scripts/pipeline/
#  OUTPUT  : the outputs of each step (see the header of each script)
# ==========================================================

source("scripts/pipeline/config.R")

steps <- data.frame(
  file = c("01_fastp_qc.R", "02_hocort_host_removal.R", "03_megahit_assembly.R",
           "04_basalt_binning.R", "05_checkm2_quality.R", "06_prepare_drep_inputs.R",
           "07_drep_stage1_per_host.R", "08_drep_stage2_species.R", "09_drep_strain_99.R",
           "10_gtdbtk_classify.R", "11_gene_prediction_prodigal.R", "12_gene_catalog_mmseqs2.R",
           "13_eggnog_annotation.R", "14_dbcan_cazyme.R", "15_diamond_nr_unassigned.R"),
  run  = c("step_fastp_qc", "step_hocort_map", "step_megahit",
           "step_basalt", "step_checkm2", "step_prepare_drep",
           "step_drep_stage1", "step_drep_stage2", "step_drep_strain",
           "step_gtdbtk", "step_prodigal", "step_gene_catalog",
           "step_eggnog", "step_dbcan", "step_diamond_nr"),
  stringsAsFactors = FALSE
)

for (f in steps$file) source(file.path(PIPELINE_DIR, f))

for (i in seq_len(nrow(steps))) {
  log_msg("===== Step ", sprintf("%02d", i), ": ", steps$file[i], " =====")
  get(steps$run[i])()
}
log_msg("Pipeline finished")
