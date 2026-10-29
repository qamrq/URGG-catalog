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
#  SCRIPT  : 07_drep_stage1_per_host.R  -  dRep stage 1: per-host dereplication at 95% ANI
#  AUTHOR  : QRQ
#  INPUT   : DREP_DIR/metadata.tsv, DREP_DIR/genomeInfo.csv, DREP_DIR/genomes/
#  OUTPUT  : DREP_DIR/stage1_per_host/<subgroup>_drep/dereplicated_genomes/
#             DREP_DIR/stage1_per_host/subgroup_index.tsv
# ==========================================================
#
# Problem: one 95% dRep run on all 117,145 genomes ran out of memory (about
# 508 GB). Genomes are split by host, and any host above MAX_PER_GROUP genomes
# is split again by BioProject. Each subgroup is dereplicated separately.
# Subgroups that already have output are skipped, so a rerun only redoes the
# failed ones.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

MAX_PER_GROUP <- 15000

# Build the dRep argument list for one run (manuscript settings)
drep_args <- function(out_dir, genome_info, genome_list, threads) {
  c("dereplicate", out_dir,
    "-p", threads,
    "--genomeInfo", genome_info,
    "--ignoreGenomeQuality",
    "-pa", "0.95", "-sa", "0.95", "-nc", "0.30",
    "-cm", "larger",
    "--S_algorithm", "skani",
    "--multiround_primary_clustering",
    "--primary_chunksize", "2000",
    "-g", genome_list)
}

# Write one genome list per subgroup and return an index (subgroup, list)
make_subgroups <- function(meta, genome_dir, lists_dir, max_per_group = MAX_PER_GROUP) {
  ensure_dir(lists_dir)
  index <- list()
  for (host in unique(meta$host)) {
    h <- meta[meta$host == host, ]
    if (nrow(h) <= max_per_group) {
      subs <- list(list(name = paste0("host__", safe_name(host)), rows = h))
      log_msg("   ", host, " (", nrow(h), ") -> single subgroup")
    } else {
      log_msg("   ", host, " (", nrow(h), ") exceeds ", max_per_group, ", splitting by BioProject:")
      subs <- lapply(unique(h$bioproject), function(bp) {
        list(name = paste0("host__", safe_name(host), "__bp__", safe_name(bp)),
             rows = h[h$bioproject == bp, ])
      })
    }
    for (s in subs) {
      list_path <- file.path(lists_dir, paste0(s$name, ".txt"))
      writeLines(file.path(genome_dir, s$rows$genome), list_path)
      index[[length(index) + 1]] <- data.frame(subgroup = s$name, list = list_path,
                                               n = nrow(s$rows), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, index)
}

step_drep_stage1 <- function(workdir = DREP_DIR, threads = THREADS) {
  metadata    <- file.path(workdir, "metadata.tsv")
  genome_info <- file.path(workdir, "genomeInfo.csv")
  genome_dir  <- file.path(workdir, "genomes")
  stage1      <- ensure_dir(workdir, "stage1_per_host")
  lists_dir   <- ensure_dir(stage1, "genome_lists")

  meta <- read.delim(metadata, stringsAsFactors = FALSE, colClasses = "character")
  log_msg("STAGE 1: per-host dereplication (95% ANI), max per group = ", MAX_PER_GROUP)

  index <- make_subgroups(meta, genome_dir, lists_dir)
  write.table(index[, c("subgroup", "list")], file.path(stage1, "subgroup_index.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg(">> Created ", nrow(index), " subgroups. Index: ", file.path(stage1, "subgroup_index.tsv"))

  done <- 0; skipped <- 0; failed <- 0
  for (i in seq_len(nrow(index))) {
    sub <- index$subgroup[i]; list_path <- index$list[i]; n <- index$n[i]
    out <- file.path(stage1, paste0(sub, "_drep"))
    reps_dir <- file.path(out, "dereplicated_genomes")
    log_msg("── [", i, "/", nrow(index), "] ", sub, " (", n, " genomes) ──")

    if (count_files(reps_dir) > 0) {
      log_msg("   already done (", count_files(reps_dir), " reps), skipping")
      skipped <- skipped + 1
      next
    }
    if (dir.exists(out)) unlink(out, recursive = TRUE)

    # single-genome subgroups: dRep cannot dereplicate one genome, so copy it
    if (n <= 1) {
      ensure_dir(reps_dir)
      if (n == 1) file.copy(readLines(list_path), reps_dir)
      log_msg("   ", n, " genome(s) copied directly, no dRep needed")
      done <- done + 1
      next
    }

    status <- run_tool("dRep", drep_args(out, genome_info, list_path, threads),
                       log = file.path(stage1, paste0(sub, ".log")), must = FALSE)
    if (identical(as.integer(status), 0L)) {
      log_msg("   done: ", count_files(reps_dir), " reps")
      done <- done + 1
    } else {
      log_msg("   FAILED, see ", file.path(stage1, paste0(sub, ".log")))
      failed <- failed + 1
    }
  }

  log_msg("STAGE 1 complete. Done: ", done, "  Skipped: ", skipped, "  Failed: ", failed)
  if (failed > 0) log_msg("Some subgroups failed. Rerun to retry only the failed ones.")
  else log_msg("All subgroups done. Proceed to step 08.")
  invisible(index)
}

if (sys.nframe() == 0L) step_drep_stage1()
