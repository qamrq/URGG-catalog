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
#  SCRIPT  : 02_hocort_host_removal.R  -  Host and diet read removal (HoCoRT)
#  AUTHOR  : QRQ
#  INPUT   : QC_DIR/<sample>_1.fastq.gz, <sample>_2.fastq.gz
#             REF_DIR/*.fasta.gz (host and plant references)
#  OUTPUT  : INDEX_DIR/<name> (Bowtie2 indexes)
#             HOCORT_DIR/<sample>_1_clean.fastq.gz, <sample>_2_clean.fastq.gz
# ==========================================================
#
# Methods (manuscript): reads were mapped to the matching host genome and dietary
# plant references with HoCoRT. Unmapped read pairs are the clean reads used for
# assembly (step 03). Part A builds one Bowtie2 index per reference set;
# Part B maps each sample to its index and keeps the unmapped reads.

# Part A: build one Bowtie2 index for a reference set.
#   name  : index prefix, e.g. "bharal_sheep_human_plant"
#   fasta : reference FASTA (.fasta.gz)
build_hocort_index <- function(name, fasta) {
  ensure_dir(INDEX_DIR)
  run_tool("hocort", c("index", "bowtie2", "--input", fasta,
                       "--output", file.path(INDEX_DIR, name)))
}

# Part B: map every QC'd sample to one index and keep the unmapped pairs.
#   index : index prefix to use for these samples (host-specific)
step_hocort_map <- function(index = file.path(INDEX_DIR, "panhomo_matv1")) {
  ensure_dir(HOCORT_DIR, "log_files")
  r1s <- list.files(QC_DIR, pattern = "_1\\.fastq\\.gz$", full.names = TRUE)
  if (length(r1s) == 0) stop("No *_1.fastq.gz files in ", QC_DIR)

  for (r1 in r1s) {
    sample <- sub("_1\\.fastq\\.gz$", "", basename(r1))
    r2   <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", r1)
    out1 <- file.path(HOCORT_DIR, paste0(sample, "_1_clean.fastq.gz"))
    out2 <- file.path(HOCORT_DIR, paste0(sample, "_2_clean.fastq.gz"))
    if (file.exists(out2) && file.size(out2) > 0) next   # already done

    run_tool("hocort", c("map", "bowtie2", "-x", index,
                         "-i", r1, r2, "--filter", "true",
                         "-o", out1, out2,
                         "preset='end-to-end'",
                         paste0("--log-file=", file.path(HOCORT_DIR, "log_files",
                                                         paste0(sample, "_map.log")))))
    log_msg("HoCoRT complete: ", sample)
  }
  log_msg("HoCoRT host removal complete: ", HOCORT_DIR)
}
