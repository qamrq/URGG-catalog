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
#  SCRIPT  : make_mock_63.R  -  Mock thiamine table for layout tests
#  AUTHOR  : QRQ
#  INPUT   : representative_host_range.csv (columns: representative, host_species_list)
#  OUTPUT  : mock_63_genome_thiamine_guilds.tsv (layout test only, not real data)
# ==========================================================
#
# Layout-test stand-in for 63_genome_thiamine_guilds.tsv. NOT real data. It is
# used only to test the figure layout and must never be used for results.

HOST_RANGE_CSV <- Sys.getenv("HOST_RANGE_CSV",
  "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig2/supp/representative_host_range.csv")
OUT_TSV <- Sys.getenv("MOCK_OUT",
  "C:/qrq_metagenome_pipeline/test/mock_63_genome_thiamine_guilds.tsv")

make_mock_63 <- function(host_range_csv = HOST_RANGE_CSV, out_tsv = OUT_TSV, seed = 1) {
  set.seed(seed)
  sh <- read.csv(host_range_csv, stringsAsFactors = FALSE)[, c("representative", "host_species_list")]
  names(sh)[names(sh) == "representative"] <- "genome"
  sh <- sh[!duplicated(sh$genome), ]
  sh <- head(sh, 25451)
  n <- nrow(sh)

  i <- seq_len(n)
  tiers <- ifelse(i <= 454, "Core / broad (>3 host species)",
           ifelse(i <= 454 + 4811, "Shared (2-3 host species)", "Host-specific (1 host species)"))

  guilds <- c("1_de_novo", "2_partial_de_novo", "3_precursor_salvage",
              "4_thiamine_salvage", "5_none_detected")
  share <- list("Core / broad (>3 host species)"     = c(.12, .19, .14, .33, .22),
                "Shared (2-3 host species)"          = c(.044, .16, .12, .39, .28),
                "Host-specific (1 host species)"     = c(.025, .10, .10, .44, .33))
  guild <- vapply(tiers, function(t) sample(guilds, 1, prob = share[[t]] / sum(share[[t]])),
                  character(1), USE.NAMES = FALSE)

  out <- data.frame(
    genome = sh$genome,
    guild = guild,
    sharing_category = tiers,
    hmp_denovo = as.integer(guild %in% guilds[1:2]),
    transport = as.integer(runif(n) < 0.5),
    tpk_thiN = as.integer(runif(n) < 0.6),
    tmp_kinase_thiL = as.integer(runif(n) < 0.3),
    stringsAsFactors = FALSE)

  write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  print(table(out$guild))
  invisible(out)
}

if (sys.nframe() == 0L) make_mock_63()
