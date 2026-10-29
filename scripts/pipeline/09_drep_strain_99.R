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
#  SCRIPT  : 09_drep_strain_99.R  -  dRep strain set at 99% ANI (near-complete genomes)
#  AUTHOR  : QRQ
#  INPUT   : DREP_DIR/metadata.tsv (quality_tier == "HQ"), DREP_DIR/genomeInfo.csv
#  OUTPUT  : DREP_DIR/drep_99_strain/dereplicated_genomes/ (23,999 strain representatives)
#             DREP_DIR/drep_99_strain/data_tables/Cdb.csv, Wdb.csv, Sdb.csv, Ndb.csv
# ==========================================================
#
# Strain-level (99% ANI) clustering of the near-complete genomes (completeness
# > 90%, contamination <= 5%). Used for the viral analyses.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_drep_strain <- function(workdir = DREP_DIR, threads = 24) {
  meta <- read.delim(file.path(workdir, "metadata.tsv"), stringsAsFactors = FALSE,
                     colClasses = "character")
  hq <- meta[meta$quality_tier == "HQ", ]
  genome_list <- file.path(workdir, "genome_list_HQ.txt")
  writeLines(file.path(workdir, "genomes", hq$genome), genome_list)
  log_msg("Near-complete (HQ) genomes: ", nrow(hq))

  out <- file.path(workdir, "drep_99_strain")
  if (dir.exists(out)) stop("Output directory exists: ", out, "\nDelete it to start fresh.")

  run_tool("dRep", c("dereplicate", out,
                     "-p", threads,
                     "--genomeInfo", file.path(workdir, "genomeInfo.csv"),
                     "--ignoreGenomeQuality",
                     "-pa", "0.95", "-sa", "0.99", "-nc", "0.30",
                     "-cm", "larger",
                     "--S_algorithm", "skani",
                     "--multiround_primary_clustering",
                     "--primary_chunksize", "2000",
                     "-g", genome_list),
           log = file.path(workdir, "drep_99_strain_run.log"))

  log_msg("Strain representatives: ", file.path(out, "dereplicated_genomes"),
          " (", count_files(file.path(out, "dereplicated_genomes")), ")")
}

if (sys.nframe() == 0L) step_drep_strain()
