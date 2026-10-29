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
#  SCRIPT  : 08_drep_stage2_species.R  -  dRep stage 2: cross-host species catalog at 95% ANI
#  AUTHOR  : QRQ
#  INPUT   : DREP_DIR/stage1_per_host/<subgroup>_drep/dereplicated_genomes/*.fa (step 07)
#             DREP_DIR/genomeInfo.csv
#  OUTPUT  : DREP_DIR/drep_95_species_final/dereplicated_genomes/ (25,451 species representatives)
#             DREP_DIR/drep_95_species_final/data_tables/Cdb.csv, Wdb.csv
# ==========================================================
#
# Takes all Stage-1 representatives, pools them, and runs one final 95% dRep to
# collapse redundancy across hosts and BioProjects.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_drep_stage2 <- function(workdir = DREP_DIR, threads = THREADS) {
  stage1      <- file.path(workdir, "stage1_per_host")
  genome_info <- file.path(workdir, "genomeInfo.csv")
  pooled_list <- file.path(workdir, "stage2_pooled_reps.txt")
  out         <- file.path(workdir, "drep_95_species_final")
  log_file    <- file.path(workdir, "drep_95_stage2_run.log")

  # Step 1: collect all Stage-1 representatives
  index <- read.delim(file.path(stage1, "subgroup_index.tsv"), stringsAsFactors = FALSE)
  reps <- character()
  n_sub <- 0
  for (sg in index$subgroup) {
    rdir <- file.path(stage1, paste0(sg, "_drep"), "dereplicated_genomes")
    if (dir.exists(rdir)) {
      reps <- c(reps, list.files(rdir, pattern = "\\.fa$", full.names = TRUE))
      n_sub <- n_sub + 1
    } else {
      warning("No representatives found for subgroup '", sg, "' at ", rdir)
    }
  }
  writeLines(reps, pooled_list)
  log_msg("Pooled ", length(reps), " representatives from ", n_sub, " subgroups")
  if (length(reps) == 0) stop("No Stage-1 representatives found. Did step 07 complete?")

  # Step 2: refuse to overwrite an earlier result
  if (dir.exists(out)) stop("Output directory exists: ", out,
                            "\nDelete it to start fresh.")

  log_msg("STAGE 2: final cross-host dereplication (95% ANI), input ", length(reps),
          " genomes, threads ", threads)
  run_tool("dRep", c("dereplicate", out,
                     "-p", threads,
                     "--genomeInfo", genome_info,
                     "--ignoreGenomeQuality",
                     "-pa", "0.95", "-sa", "0.95", "-nc", "0.30",
                     "-cm", "larger",
                     "--S_algorithm", "skani",
                     "--multiround_primary_clustering",
                     "--primary_chunksize", "2000",
                     "-g", pooled_list),
           log = log_file)

  log_msg("STAGE 2 complete. Final species catalog: ",
          file.path(out, "dereplicated_genomes"),
          " (", count_files(file.path(out, "dereplicated_genomes")), " representatives)")
  log_msg("Cluster assignments: ", file.path(out, "data_tables", "Cdb.csv"))
  log_msg("To map an original genome to its final species cluster, chain the Stage-1 Cdb",
          " into the Stage-2 Cdb: Stage 2 only saw Stage-1 winners.")
}

if (sys.nframe() == 0L) step_drep_stage2()
