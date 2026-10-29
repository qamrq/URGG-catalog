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
#  SCRIPT  : 03_megahit_assembly.R  -  Per-sample metagenome assembly (MEGAHIT v1.2.9)
#  AUTHOR  : QRQ
#  INPUT   : HOCORT_DIR/<sample>_1_clean.fastq.gz, <sample>_2_clean.fastq.gz
#  OUTPUT  : ASM_DIR/<sample>_megahit_output/final.contigs.fa
# ==========================================================
#
# Settings (manuscript): meta-sensitive preset, contigs >= 1,000 bp,
# k-list 27,47,67,87,107,127 (from megahit_assembly.txt).

step_megahit <- function(threads = THREADS, mem = 0.8) {
  ensure_dir(ASM_DIR)
  r1s <- list.files(HOCORT_DIR, pattern = "_1_clean\\.fastq\\.gz$", full.names = TRUE)
  if (length(r1s) == 0) stop("No *_1_clean.fastq.gz files in ", HOCORT_DIR)

  for (r1 in r1s) {
    sample <- sub("_1_clean\\.fastq\\.gz$", "", basename(r1))
    r2  <- sub("_1_clean\\.fastq\\.gz$", "_2_clean.fastq.gz", r1)
    out <- file.path(ASM_DIR, paste0(sample, "_megahit_output"))
    if (file.exists(file.path(out, "final.contigs.fa"))) next   # already done

    run_tool("megahit", c("-1", r1, "-2", r2,
                          "--k-list", "27,47,67,87,107,127",
                          "--min-contig-len", "1000",
                          "--presets", "meta-sensitive",
                          "--memory", mem,
                          "--num-cpu-threads", threads,
                          "--continue",
                          "--out-dir", out))
    log_msg("Assembly completed: ", sample)
  }
  log_msg("All MEGAHIT assemblies complete: ", ASM_DIR)
}
