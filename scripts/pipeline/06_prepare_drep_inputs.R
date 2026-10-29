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
#  SCRIPT  : 06_prepare_drep_inputs.R  -  Quality gate and dRep inputs
#  AUTHOR  : QRQ
#  INPUT   : ROOT/source_data/list_bioproject (BioProject and host per line)
#             BASALT_DIR/<BioProject>/BestBinset/Best_binset_quality_report.tsv and bin folder
#  OUTPUT  : WORKDIR/genomes/ (renamed genomes)
#             WORKDIR/genome_list.txt, metadata.tsv, genomeInfo.csv
# ==========================================================
#
# Processes BASALT Best_binset_quality_report.tsv files for each BioProject in
# list_bioproject. Applies the quality gate, tags each genome HQ or MQ, copies
# passing genomes into one flat pool with origin-encoding names, and writes the
# three files dRep needs (genome_list.txt, metadata.tsv, genomeInfo.csv).
# Genomes are tagged by quality_tier, not split into HQ and MQ folders.

WORKDIR        <- DREP_DIR      # shared dRep folder used by steps 06-09
PROJECT_LIST   <- file.path(ROOT, "source_data", "list_bioproject")
CATALOG_MIN_COMP <- 70    # floor to enter the pool at all
HQ_MIN_COMP      <- 90    # > 90 => HQ / strain-catalog eligible
MAX_CONT         <- 5
BIN_EXTENSIONS   <- c("", ".fa", ".fasta", ".fna")

step_prepare_drep <- function(project_list = PROJECT_LIST, workdir = WORKDIR) {
  projects <- read.delim(project_list, comment.char = "#", stringsAsFactors = FALSE,
                         colClasses = "character")
  projects <- projects[projects$bioproject != "bioproject" & nzchar(projects$bioproject), ]
  log_msg("Projects in ", project_list, ": ", nrow(projects))

  genome_dir <- ensure_dir(workdir, "genomes")
  rows <- list()
  total <- 0L

  for (i in seq_len(nrow(projects))) {
    bp   <- projects$bioproject[i]
    host <- projects$host[i]
    qtsv <- file.path(BASALT_DIR, bp, "BestBinset", "Best_binset_quality_report.tsv")
    bindir <- file.path(BASALT_DIR, bp, "BestBinset")
    log_msg(">> ", bp, " (", host, ")  table=", qtsv)
    if (!file.exists(qtsv)) { warning("Quality table not found: ", qtsv); next }

    q <- read.delim(qtsv, stringsAsFactors = FALSE, colClasses = c(Bin_ID = "character"))
    total <- total + nrow(q)

    for (k in seq_len(nrow(q))) {
      bin_id <- q$Bin_ID[k]
      comp   <- as.numeric(q$Completeness[k])
      cont   <- as.numeric(q$Contamination[k])
      n50    <- q$N50[k]
      if (is.na(comp) || is.na(cont)) next
      if (comp < CATALOG_MIN_COMP || cont > MAX_CONT) next

      tier <- if (comp > HQ_MIN_COMP) "HQ" else "MQ"

      # BASALT bin_id is usually the file name, sometimes without the extension
      cands <- file.path(bindir, paste0(bin_id, BIN_EXTENSIONS))
      src <- cands[file.exists(cands)][1]
      if (is.na(src)) {
        warning("Bin file not found for '", bin_id, "' under ", bindir)
        next
      }

      newname <- paste0(bp, "__", safe_name(host), "__", safe_name(bin_id), ".fa")
      file.copy(src, file.path(genome_dir, newname), overwrite = TRUE)

      rows[[length(rows) + 1]] <- data.frame(
        genome = newname, bioproject = bp, host = host, quality_tier = tier,
        completeness = comp, contamination = cont, N50 = n50,
        stringsAsFactors = FALSE)
    }
  }

  meta <- if (length(rows)) do.call(rbind, rows) else
    data.frame(genome = character(), bioproject = character(), host = character(),
               quality_tier = character(), completeness = numeric(),
               contamination = numeric(), N50 = numeric(), stringsAsFactors = FALSE)

  write.table(meta, file.path(workdir, "metadata.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  writeLines(file.path(genome_dir, meta$genome), file.path(workdir, "genome_list.txt"))
  write.csv(data.frame(genome = meta$genome, completeness = meta$completeness,
                       contamination = meta$contamination),
            file.path(workdir, "genomeInfo.csv"), row.names = FALSE, quote = FALSE)

  kept <- nrow(meta)
  hq   <- sum(meta$quality_tier == "HQ")
  log_msg("-------------------------------------------")
  log_msg("Scanned:      ", total, " bins")
  log_msg("Passed (>=", CATALOG_MIN_COMP, " comp, <=", MAX_CONT, " cont): ", kept)
  log_msg("  of which HQ (>", HQ_MIN_COMP, " comp): ", hq)
  log_msg("Pool genomes: ", genome_dir)
  invisible(meta)
}
