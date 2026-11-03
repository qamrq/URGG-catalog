#!/usr/bin/env Rscript
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
#  SCRIPT  : 07_rebuild_ring_annotation_lineage76.R  -  Ring annotation table for the RUC tree
#  AUTHOR  : QRQ
#  INPUT   : lineage genome IDs, rerun GTDB-Tk summary, CAZy and cooperative presence tables
#  OUTPUT  : annotation/iqtree_output/tree_build_lineage_only_fresh/tree_ring_annotation_lineage76.csv
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages(library(tidyverse))

genome_list_path <- paste0(URGG_PROJECT, "/annotation/results_phase3/s32_taxonomy_novelty/lineage_genome_ids.txt")

gtdbtk_rerun_summary <- paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232_lineage_rerun/classify/gtdbtk.bac120.summary.tsv")

cazy_wide_path <- paste0(URGG_PROJECT, "/annotation/iqtree_output/genome_tier_cazy_wide.csv")

cooperative_path <- paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_category_presence_per_genome.csv")

out_path <- paste0(URGG_PROJECT, "/annotation/iqtree_output/tree_build_lineage_only_fresh/tree_ring_annotation_lineage76.csv")

genome_ids <- read_lines(genome_list_path) %>% discard(~ .x == "")
stopifnot(length(genome_ids) == 76)

# --- Ring 1: Genus (from classification string, g__ field) ---
gtdbtk <- read_tsv(gtdbtk_rerun_summary, show_col_types = FALSE) %>%
  select(genome = user_genome, classification) %>%
  filter(genome %in% genome_ids) %>%
  mutate(
    genus = str_extract(classification, "g__[^;]*") %>% str_remove("^g__"),
    genus = if_else(genus == "" | is.na(genus), "Unclassified", genus)
  ) %>%
  select(genome, genus)

cat(sprintf("Genus: matched %d/76 genomes\n", nrow(gtdbtk)))
cat("Genus distribution:\n")
print(count(gtdbtk, genus, sort = TRUE))

# --- Ring 2: CAZyme headline families ---
cazy_headline_families <- c("GH190", "GT121", "GT26")

cazy <- read_csv(cazy_wide_path, show_col_types = FALSE) %>%
  filter(genome %in% genome_ids) %>%
  select(genome, any_of(cazy_headline_families))

missing_cazy_cols <- setdiff(cazy_headline_families, colnames(cazy))
if (length(missing_cazy_cols) > 0) {
  warning(sprintf(
    "CAZy columns not found in source table: %s — check exact column naming in genome_tier_cazy_wide.csv",
    paste(missing_cazy_cols, collapse = ", ")
  ))
}

cat(sprintf("\nCAZyme: matched %d/76 genomes\n", nrow(cazy)))

# --- Ring 3: Cooperative-trait categories ---
cooperative_categories <- c(
  "antibiotic_degradation", "biofilm", "Extracellular",
  "quorum_sensing", "secretion_system", "siderophores"
)

cooperative_long <- read_csv(cooperative_path, show_col_types = FALSE)

cat("\nCategory values found in source table:\n")
print(count(cooperative_long, cooperative_categories))

cooperative <- cooperative_long %>%
  filter(genome %in% genome_ids) %>%
  select(genome, cooperative_categories, present) %>%
  pivot_wider(names_from = cooperative_categories, values_from = present) %>%
  select(genome, any_of(cooperative_categories))

missing_coop_cols <- setdiff(cooperative_categories, colnames(cooperative))
if (length(missing_coop_cols) > 0) {
  warning(sprintf(
    "Cooperative-trait columns not found in source table: %s — check exact column naming",
    paste(missing_coop_cols, collapse = ", ")
  ))
}

cat(sprintf("Cooperative traits: matched %d/76 genomes\n", nrow(cooperative)))

ring_table <- tibble(genome = genome_ids) %>%
  left_join(gtdbtk, by = "genome") %>%
  left_join(cazy, by = "genome") %>%
  left_join(cooperative, by = "genome") %>%
  mutate(across(all_of(cazy_headline_families), ~ replace_na(.x, 0))) %>%
  mutate(across(all_of(cooperative_categories), ~ replace_na(.x, 0)))

to_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x != 0))
  x_upper <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x_upper %in% c("TRUE", "T", "1", "YES") ~ 1L,
    x_upper %in% c("FALSE", "F", "0", "NO") ~ 0L,
    TRUE ~ NA_integer_
  )
}

ring_table <- ring_table %>%
  mutate(across(all_of(c(cazy_headline_families, cooperative_categories)), to_binary))

n_bad_binary <- ring_table %>%
  filter(if_any(all_of(c(cazy_headline_families, cooperative_categories)), is.na)) %>%
  nrow()
if (n_bad_binary > 0) {
  warning(sprintf(
    "%d genome(s) have an unparseable binary value after normalization — inspect before writing iTOL files.",
    n_bad_binary
  ))
}

n_incomplete <- ring_table %>%
  filter(if_any(-genome, is.na)) %>%
  nrow()

cat(sprintf("\n%d/76 genomes have at least one NA across the merged ring table\n", n_incomplete))
if (n_incomplete > 0) {
  cat("Genomes with missing data (inspect before plotting):\n")
  ring_table %>%
    filter(if_any(-genome, is.na)) %>%
    pull(genome) %>%
    print()
}

write_csv(ring_table, out_path)
cat("\nWrote:", out_path, "\n")
