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
#  SCRIPT  : 12a_cazy_join_tiers.R  -  Join host-sharing tiers to the CAZy table
#  AUTHOR  : QRQ
#  INPUT   : CAZy long table, representative_host_range.csv (arguments 1-2)
#  OUTPUT  : genome_tier_cazy_wide.csv (argument 3 = output folder)
# ==========================================================

suppressMessages({
  library(dplyr); library(readr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 12a_cazy_join_tiers.R <cazy_long.csv> ",
       "<representative_host_range.csv> <output_dir>")
}
cazy_path <- args[1]
tier_path <- args[2]
outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cazy <- read_csv(cazy_path, show_col_types = FALSE)
tiers <- read_csv(tier_path, show_col_types = FALSE)

if (!"representative" %in% colnames(tiers)) {
  stop("representative_host_range.csv missing expected 'representative' column. ",
       "Actual columns: ", paste(colnames(tiers), collapse = ", "))
}
if (!"genome" %in% colnames(cazy)) {
  stop("CAZy long table missing expected 'genome' column. ",
       "Actual columns: ", paste(colnames(cazy), collapse = ", "))
}

cazy_genomes <- unique(cazy$genome)
tier_genomes <- unique(tiers$representative)
overlap <- length(intersect(cazy_genomes, tier_genomes))

message(sprintf("CAZy table genomes: %d", length(cazy_genomes)))
message(sprintf("Tier table genomes (representatives): %d", length(tier_genomes)))
message(sprintf("Overlap on direct ID match: %d", overlap))

if (overlap == 0) {
  stop("ZERO overlap between CAZy genome IDs and representative_host_range.csv's ",
       "'representative' column. These do not use the same naming convention ",
       "(e.g. one may have file extensions or path prefixes the other lacks). ",
       "Reconcile the ID formats first.\n",
       "Example CAZy genome ID: ", cazy_genomes[1], "\n",
       "Example tier representative ID: ", tier_genomes[1])
} else if (overlap < 0.5 * min(length(cazy_genomes), length(tier_genomes))) {
  warning(sprintf(
    "Only %.1f%% of the smaller table's genomes matched; check the ID formats.",
    100 * overlap / min(length(cazy_genomes), length(tier_genomes))))
}

joined_long <- cazy %>%
  inner_join(tiers, by = c("genome" = "representative"))

n_joined <- n_distinct(joined_long$genome)
message(sprintf("Genomes successfully joined (long form): %d", n_joined))

write_csv(joined_long, file.path(outdir, "genome_tier_cazy_joined.csv"))

cazy_wide <- cazy %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = cazy_family, values_from = present, values_fill = 0L)

joined_wide <- cazy_wide %>%
  inner_join(tiers, by = c("genome" = "representative"))

write_csv(joined_wide, file.path(outdir, "genome_tier_cazy_wide.csv"))

message("\nOutputs written to: ", normalizePath(outdir))
message("  genome_tier_cazy_joined.csv  (", nrow(joined_long), " rows, long form)")
message("  genome_tier_cazy_wide.csv    (", nrow(joined_wide), " genomes, wide form)")
message("\nNEXT: spot-check this joined table reproduces known published stats")
message("(e.g. CBM98/GT117/GH190 odds ratios by sharing_category) before using")
message("it as the basis for tree genome selection or ring annotation.")
