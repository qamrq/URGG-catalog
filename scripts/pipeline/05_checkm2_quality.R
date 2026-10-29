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
#  SCRIPT  : 05_checkm2_quality.R  -  Completeness and contamination (CheckM2 v1.1.0)
#  AUTHOR  : QRQ
#  INPUT   : BIN_DIR/*.fa (BASALT bins)
#  OUTPUT  : QC_OUT/quality_report.tsv
# ==========================================================
#
# Threshold applied downstream (step 06): >= 70% completeness, <= 5% contamination.

BIN_DIR <- file.path(BASALT_DIR, "bins")

step_checkm2 <- function(threads = 30) {
  ensure_dir(QC_OUT)
  run_tool("checkm2", c("predict", "--threads", threads,
                        "--input", BIN_DIR,
                        "--output-directory", QC_OUT,
                        "-x", "fa"))
  log_msg("CheckM2 complete: ", file.path(QC_OUT, "quality_report.tsv"))
}
