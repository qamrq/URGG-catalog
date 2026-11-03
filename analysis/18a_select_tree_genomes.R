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
#  SCRIPT  : 18a_select_tree_genomes.R  -  Select tier-balanced genomes for the subset tree
#  AUTHOR  : QRQ
#  INPUT   : genome_tier_cazy_wide.csv, output folder (arguments)
#  OUTPUT  : tree genome list
# ==========================================================

suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 18a_select_tree_genomes.R <genome_tier_cazy_wide.csv> ",
       "<output_dir> [fam_substrate_mapping.tsv] [target_n]")
}
wide_path <- args[1]
outdir <- args[2]
subst_map_path <- if (length(args) >= 3) args[3] else NA
target_n <- if (length(args) >= 4) as.integer(args[4]) else 1200L

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

wide <- read_csv(wide_path, show_col_types = FALSE)
message("Loaded ", nrow(wide), " genomes, ", ncol(wide), " columns.")

if (!"sharing_category" %in% colnames(wide)) {
  stop("Column 'sharing_category' not found. Actual columns: ",
       paste(head(colnames(wide), 30), collapse = ", "), " ...")
}
core_label <- grep("core", unique(wide$sharing_category), ignore.case = TRUE, value = TRUE)
if (length(core_label) != 1) {
  stop("Could not uniquely detect core/broad label. Candidates: ",
       paste(core_label, collapse = ", "))
}
message("Core/broad label: '", core_label, "'")

tier_counts <- wide %>% count(sharing_category)
message("\nTier counts:")
print(tier_counts)

phylum_col <- NULL
candidates <- grep("phylum", colnames(wide), ignore.case = TRUE, value = TRUE)
if (length(candidates) >= 1) {
  phylum_col <- candidates[1]
  message("\nUsing taxonomy column: '", phylum_col, "'")
} else {
  lineage_candidates <- grep("classification|lineage|taxonomy", colnames(wide),
                             ignore.case = TRUE, value = TRUE)
  if (length(lineage_candidates) >= 1) {
    lin_col <- lineage_candidates[1]
    message("\nNo direct phylum column found. Attempting to parse phylum from '",
            lin_col, "' (expects GTDB-style 'd__...;p__...;c__...' format).")
    wide <- wide %>%
      mutate(phylum_parsed = str_extract(.data[[lin_col]], "p__[^;]+"))
    if (all(is.na(wide$phylum_parsed))) {
      stop("Could not parse a phylum field from '", lin_col, "'. ",
           "Inspect this column's actual format before proceeding -- do not ",
           "guess further. Sample values: ",
           paste(head(unique(wide[[lin_col]]), 3), collapse = " | "))
    }
    phylum_col <- "phylum_parsed"
  } else {
    stop("No usable taxonomy column found (checked for 'phylum', 'classification', ",
         "'lineage', 'taxonomy'). Actual columns: ",
         paste(colnames(wide), collapse = ", "),
         "\nSTOP and identify the correct taxonomy column before proceeding.")
  }
}

n_phyla <- n_distinct(wide[[phylum_col]], na.rm = TRUE)
message("Distinct phyla detected: ", n_phyla)
print(wide %>% count(.data[[phylum_col]], sort = TRUE) %>% head(20))

core_genomes <- wide %>% filter(sharing_category == core_label)
n_core <- nrow(core_genomes)
message("\n>> Including all ", n_core, " core/broad genomes.")

remaining_budget <- max(target_n - n_core, 0)
non_core <- wide %>% filter(sharing_category != core_label)

if (remaining_budget > 0 && nrow(non_core) > 0) {

  tiers_present <- unique(non_core$sharing_category)
  n_tiers <- length(tiers_present)
  per_tier_budget <- floor(remaining_budget / n_tiers)

  sampled_list <- list()

  for (tier in tiers_present) {

    # Correct: define tier_genomes inside the loop
    tier_genomes <- non_core %>% filter(sharing_category == tier)

    # Identify phyla present in this tier
    phyla_in_tier <- unique(tier_genomes[[phylum_col]])
    phyla_in_tier <- phyla_in_tier[!is.na(phyla_in_tier)]
    n_phyla_tier <- length(phyla_in_tier)

    # Budget per phylum
    per_phylum_budget <- ceiling(per_tier_budget / max(n_phyla_tier, 1))

    # Correct sampling using group_modify()
    tier_sample <- tier_genomes %>%
      filter(!is.na(.data[[phylum_col]])) %>%
      group_by(.data[[phylum_col]]) %>%
      group_modify(~ {
        n_group <- nrow(.x)
        slice_sample(.x, n = min(per_phylum_budget, n_group))
      }) %>%
      ungroup()

    sampled_list[[tier]] <- tier_sample

    message(sprintf(
      "  Tier '%s': sampled %d genomes across %d phyla",
      tier, nrow(tier_sample), n_phyla_tier
    ))
  }

  non_core_sample <- bind_rows(sampled_list)

} else {
  non_core_sample <- non_core[0, ]
}

selected <- bind_rows(core_genomes, non_core_sample)

message("\n>> Total selected genomes: ", nrow(selected),
        " (target was ~", target_n, ")")

if (!"genome" %in% colnames(selected)) {
  stop("Expected 'genome' column not found in wide table.")
}
writeLines(selected$genome, file.path(outdir, "tree_genome_list.txt"))

ring_cols <- c("genome", "sharing_category")
for (fam in c("CBM98", "GT117", "GH190")) {
  if (fam %in% colnames(selected)) {
    ring_cols <- c(ring_cols, fam)
  } else {
    warning("Family column '", fam, "' not found in wide table -- skipping this ring.")
  }
}
ring_table <- selected %>% select(all_of(ring_cols))

if (!is.na(subst_map_path) && file.exists(subst_map_path)) {
  message("\n>> Mapping families to substrates using: ", subst_map_path)
  subst_map <- read_tsv(subst_map_path, show_col_types = FALSE)
  message("   Substrate mapping columns: ", paste(colnames(subst_map), collapse = ", "))

  subst_long <- subst_map %>%
    separate_rows(Substrate_high_level, sep = ";") %>%
    mutate(Substrate_high_level = str_trim(tolower(Substrate_high_level)))

  get_families_for_substrate <- function(target, subst_long_tbl, subst_map_tbl) {
    fams <- subst_long_tbl %>%
      filter(Substrate_high_level == tolower(target)) %>%
      pull(Family) %>% unique()

    if (length(fams) == 0) {
      message("   NOTE: '", target, "' not found as a standalone Substrate_high_level ",
              "token -- falling back to word-boundary search in Substrate_curated.")
      pattern <- paste0("\\b", target, "\\b")
      fams <- subst_map_tbl %>%
        filter(str_detect(tolower(Substrate_curated), pattern)) %>%
        pull(Family) %>% unique()
    }
    fams
  }

  chitin_families <- get_families_for_substrate("chitin", subst_long, subst_map)
  pectin_families <- get_families_for_substrate("pectin", subst_long, subst_map)

  message("\n   Chitin-associated families (n=", length(chitin_families), "): ",
          paste(head(chitin_families, 15), collapse = ", "),
          if (length(chitin_families) > 15) " ..." else "")
  message("   Pectin-associated families (n=", length(pectin_families), "): ",
          paste(head(pectin_families, 15), collapse = ", "),
          if (length(pectin_families) > 15) " ..." else "")
  message("\n   VERIFY the family lists above look biologically sensible before",
          " trusting the resulting substrate rings.")

  available_cols <- colnames(selected)
  chitin_cols_present <- intersect(chitin_families, available_cols)
  pectin_cols_present <- intersect(pectin_families, available_cols)

  message("\n   Of the mapped families, present as columns in the wide table:")
  message("     Chitin: ", length(chitin_cols_present), " of ", length(chitin_families))
  message("     Pectin: ", length(pectin_cols_present), " of ", length(pectin_families))

  if (length(chitin_cols_present) > 0) {
    selected <- selected %>%
      mutate(chitin_capacity = as.integer(
        rowSums(across(all_of(chitin_cols_present))) > 0))
    ring_table <- selected %>% select(all_of(ring_cols), chitin_capacity)
  } else {
    warning("No chitin-associated families found as columns in the wide table -- ",
            "chitin ring not added.")
  }

  if (length(pectin_cols_present) > 0) {
    selected <- selected %>%
      mutate(pectin_capacity = as.integer(
        rowSums(across(all_of(pectin_cols_present))) > 0))
    ring_table <- selected %>%
      select(all_of(ring_cols),
             any_of(c("chitin_capacity", "pectin_capacity")))
  } else {
    warning("No pectin-associated families found as columns in the wide table -- ",
            "pectin ring not added.")
  }

} else {
  message("\n>> No fam-substrate-mapping.tsv provided (or path not found) -- ",
          "chitin/pectin substrate rings SKIPPED. Provide the mapping file ",
          "path as the 3rd argument to include them.")
}

write_csv(ring_table, file.path(outdir, "tree_ring_annotation.csv"))

message("\n============================================================")
message("OUTPUTS")
message("============================================================")
message("  tree_genome_list.txt     (", nrow(selected), " genome IDs, for ",
        "subsetting the GTDB-Tk alignment)")
message("  tree_ring_annotation.csv (", nrow(ring_table), " rows x ",
        ncol(ring_table), " cols: genome, tier, + available family rings)")
message("\nNEXT STEPS:")
message("  1. Subset gtdbtk.bac120.user_msa.fasta.gz / .ar53. to tree_genome_list.txt")
message("  2. trimAl")
message("  3. IQ-TREE (-m MFP -bb 1000 -redo)")
message("  4. ggtree figure using tree_ring_annotation.csv")
