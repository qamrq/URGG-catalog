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
#  SCRIPT  : 01_fastp_qc.R  -  Read quality filtering (fastp v1.0.1, default parameters)
#  AUTHOR  : QRQ
#  INPUT   : RAW_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#  OUTPUT  : QC_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#             QC_DIR/fastp_reports/<sample>.html, .json
# ==========================================================

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_fastp_qc <- function(threads = THREADS_FASTP) {
  ensure_dir(QC_DIR, "fastp_reports")
  r1s <- list.files(RAW_DIR, pattern = "_1\\.fastq\\.gz$", full.names = TRUE)
  if (length(r1s) == 0) stop("No *_1.fastq.gz files in ", RAW_DIR)

  for (r1 in r1s) {
    sample <- sub("_1\\.fastq\\.gz$", "", basename(r1))
    r2   <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", r1)
    out1 <- file.path(QC_DIR, paste0(sample, "_1.fastq.gz"))
    out2 <- file.path(QC_DIR, paste0(sample, "_2.fastq.gz"))
    if (file.exists(out2) && file.size(out2) > 0) next   # already done

    run_tool("fastp", c("-i", r1, "-I", r2, "-o", out1, "-O", out2,
                        "-h", file.path(QC_DIR, "fastp_reports", paste0(sample, ".html")),
                        "-j", file.path(QC_DIR, "fastp_reports", paste0(sample, ".json")),
                        "-w", threads))
    log_msg("QC complete: ", sample)
  }
  log_msg("fastp QC complete: ", QC_DIR)
}

THREADS_FASTP <- 16

if (sys.nframe() == 0L) step_fastp_qc()
