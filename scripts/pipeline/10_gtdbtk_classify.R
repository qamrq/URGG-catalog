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
#  SCRIPT  : 10_gtdbtk_classify.R  -  Taxonomy classification (GTDB-Tk 2.7.2, release 232)
#  AUTHOR  : QRQ
#  INPUT   : GENOME_DIR/*.fa (species representatives, step 08)
#  OUTPUT  : ANNOT/gtdbtk_output/gtdbtk.bac120.summary.tsv, gtdbtk.ar53.summary.tsv
#             ANNOT/gtdbtk_output/gtdbtk.all.summary.tsv
# ==========================================================
#
# GTDB-Tk is the heaviest single step (about 320 GB of RAM). Use --batchfile
# splitting if memory is short.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_gtdbtk <- function(threads = THREADS) {
  out <- ensure_dir(ANNOT, "gtdbtk_output")
  Sys.setenv(GTDBTK_DATA_PATH = GTDBTK_DB)
  n <- length(list.files(GENOME_DIR, pattern = "\\.fa$"))
  log_msg("GTDB-Tk classify_wf on ", n, " genomes, reference ", GTDBTK_DB)

  run_tool("gtdbtk", c("classify_wf",
                       "--genome_dir", GENOME_DIR,
                       "--out_dir", out,
                       "--extension", "fa",
                       "--cpus", threads,
                       "--pplacer_cpus", 8,
                       "--skip_ani_screen"),
           env = ENV_GTDBTK,
           log = file.path(out, "gtdbtk_run.log"))

  # Merge bacterial and archaeal summaries into one table
  bac <- file.path(out, "gtdbtk.bac120.summary.tsv")
  ar  <- file.path(out, "gtdbtk.ar53.summary.tsv")
  parts <- Filter(file.exists, c(bac, ar))
  merged <- do.call(rbind, lapply(parts, read.delim, stringsAsFactors = FALSE,
                                  colClasses = "character"))
  write.table(merged, file.path(out, "gtdbtk.all.summary.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  log_msg("Combined taxonomy: ", file.path(out, "gtdbtk.all.summary.tsv"),
          " (", nrow(merged), " genomes)")
}

if (sys.nframe() == 0L) step_gtdbtk()
