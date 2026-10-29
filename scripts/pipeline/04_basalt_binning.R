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
#  SCRIPT  : 04_basalt_binning.R  -  Binning and refinement (BASALT autobinning)
#  AUTHOR  : QRQ
#  INPUT   : ASM_DIR/*_megahit_output/final.contigs.fa (step 03)
#             HOCORT_DIR/*_1_clean.fastq.gz, *_2_clean.fastq.gz (step 02)
#  OUTPUT  : BASALT_DIR/ (bins, quality reports, best bin set)
# ==========================================================
#
# Settings (manuscript): minimum completeness 70%, maximum contamination 10%.
# Run as: BASALT ... --module autobinning --sensitive quick -t 128 -m 500
#
# Working folder layout (BASALT needs these names):
#   BASALT_DIR/<sample>_final_contigs.fa   assembly per sample (symlink to step 03)
#   BASALT_DIR/<sample>_R1.fastq.gz        clean reads, renamed from _1_clean (symlink)
#   BASALT_DIR/<sample>_R2.fastq.gz        clean reads, renamed from _2_clean (symlink)
#   BASALT_DIR/assemblies.txt, shortreads.txt
# Symlinks are used instead of the mv renaming in the original runs, so the
# HoCoRT and MEGAHIT outputs stay in place.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

# Link a file into the BASALT folder under a new name (replaces an existing link)
link_into <- function(target, link) {
  lnk <- Sys.readlink(link)                       # "" if not a link, NA on error
  is_link <- !is.na(lnk) && nzchar(lnk)
  if (is_link || file.exists(link)) file.remove(link)
  file.symlink(target, link)
}

step_basalt <- function(threads = THREADS, mem_gb = 500) {
  ensure_dir(BASALT_DIR)

  # 1) Link each MEGAHIT assembly as <sample>_final_contigs.fa
  asm_dirs <- list.files(ASM_DIR, pattern = "_megahit_output$", full.names = TRUE)
  for (d in asm_dirs) {
    sample <- sub("_megahit_output$", "", basename(d))
    link_into(file.path(d, "final.contigs.fa"),
              file.path(BASALT_DIR, paste0(sample, "_final_contigs.fa")))
  }

  # 2) Link clean reads as <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
  r1s <- list.files(HOCORT_DIR, pattern = "_1_clean\\.fastq\\.gz$", full.names = TRUE)
  for (r1 in r1s) {
    sample <- sub("_1_clean\\.fastq\\.gz$", "", basename(r1))
    link_into(r1, file.path(BASALT_DIR, paste0(sample, "_R1.fastq.gz")))
    link_into(sub("_1_clean\\.fastq\\.gz$", "_2_clean.fastq.gz", r1),
              file.path(BASALT_DIR, paste0(sample, "_R2.fastq.gz")))
  }

  # 3) Assembly list (comma-separated) and short-read pairs "R1,R2/R1,R2/..."
  asm <- sort(list.files(BASALT_DIR, pattern = "_final_contigs\\.fa$"))
  writeLines(paste(asm, collapse = ","), file.path(BASALT_DIR, "assemblies.txt"))

  r1_links <- sort(list.files(BASALT_DIR, pattern = "_R1\\.fastq\\.gz$"))
  pairs <- paste0(r1_links, ",", sub("_R1\\.fastq\\.gz$", "_R2.fastq.gz", r1_links))
  writeLines(paste(pairs, collapse = "/"), file.path(BASALT_DIR, "shortreads.txt"))

  # 4) Binning (run from the BASALT folder so the relative names resolve)
  old <- setwd(BASALT_DIR)
  on.exit(setwd(old), add = TRUE)
  run_tool("BASALT", c("-a", paste(asm, collapse = ","),
                       "-s", paste(pairs, collapse = "/"),
                       "--module", "autobinning", "--sensitive", "quick",
                       "-t", threads, "-m", mem_gb,
                       "--min-cpn", 70, "--max-ctn", 10))
  log_msg("BASALT autobinning complete: ", BASALT_DIR)
}

if (sys.nframe() == 0L) step_basalt()
