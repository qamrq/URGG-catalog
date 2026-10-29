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

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

BIN_DIR <- file.path(BASALT_DIR, "bins")   # edit if your bins are elsewhere

step_checkm2 <- function(threads = 30) {
  ensure_dir(QC_OUT)
  # One-time database download (run once, then comment out):
  # run_tool("checkm2", c("database", "--download"))
  run_tool("checkm2", c("predict", "--threads", threads,
                        "--input", BIN_DIR,
                        "--output-directory", QC_OUT,
                        "-x", "fa"))
  log_msg("CheckM2 complete: ", file.path(QC_OUT, "quality_report.tsv"))
}

if (sys.nframe() == 0L) step_checkm2()
