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
#  SCRIPT  : 11_gene_prediction_prodigal.R  -  Gene prediction per genome (Prodigal v2.6.3)
#  AUTHOR  : QRQ
#  INPUT   : GENOME_DIR/*.fa
#  OUTPUT  : GENES/<genome>.faa, GENES/<genome>.gff
# ==========================================================
#
# Prodigal in single-genome mode, with meta mode as a fallback for genomes that
# single mode rejects. Runs in parallel (one process per genome).
# Resumable: genomes that already have a .faa file are skipped.
# Note: the proposal named MetaGeneMark; Prodigal is used because MetaGeneMark
# needs a licence key and is awkward to run in batches.

predict_one <- function(fa, genes_dir) {
  base <- sub("\\.fa$", "", basename(fa))
  faa  <- file.path(genes_dir, paste0(base, ".faa"))
  gff  <- file.path(genes_dir, paste0(base, ".gff"))
  if (file.exists(faa) && file.size(faa) > 0) return(invisible(0L))

  status <- run_tool("prodigal", c("-i", fa, "-a", faa, "-o", gff, "-f", "gff", "-p", "single", "-q"),
                     must = FALSE)
  if (!identical(as.integer(status), 0L)) {
    status <- run_tool("prodigal", c("-i", fa, "-a", faa, "-o", gff, "-f", "gff", "-p", "meta", "-q"),
                       must = FALSE)
  }
  invisible(as.integer(status))
}

step_prodigal <- function(jobs = PARALLEL_JOBS) {
  genes_dir <- ensure_dir(GENES)
  fas <- list.files(GENOME_DIR, pattern = "\\.fa$", full.names = TRUE)
  log_msg("Predicting genes on ", length(fas), " genomes with ", jobs, " jobs")

  # Windows cannot fork, so run serially there (the pipeline runs on Linux)
  cores <- if (.Platform$OS.type == "windows") 1L else jobs
  status <- parallel::mclapply(fas, predict_one, genes_dir = genes_dir, mc.cores = cores)

  failed <- sum(unlist(status) != 0L)
  log_msg("Prodigal complete: ", length(list.files(genes_dir, pattern = "\\.faa$")),
          " protein files in ", genes_dir, "; failed: ", failed)
}
