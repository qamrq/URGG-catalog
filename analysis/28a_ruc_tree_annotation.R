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
#  SCRIPT  : 28a_ruc_tree_annotation.R  -  Ring annotation and iTOL datasets for the RUC tree
#  AUTHOR  : QRQ
#  INPUT   : lineage genome IDs, rerun GTDB-Tk summary, CAZy and cooperative presence tables
#  OUTPUT  : annotation/iqtree_output/tree_build_lineage_only_fresh/ (ring table, itol_datasets/)
#  USAGE   : Rscript 28a_ruc_tree_annotation.R [rings itol]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("rings")) local({
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
})

if (run_part("itol")) local({
suppressPackageStartupMessages(library(tidyverse))

ring_table_path <- paste0(URGG_PROJECT, "/annotation/iqtree_output/tree_build_lineage_only_fresh/tree_ring_annotation_lineage76.csv")

out_dir <- paste0(URGG_PROJECT, "/annotation/iqtree_output/tree_build_lineage_only_fresh/itol_datasets")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ring_table <- read_csv(ring_table_path, show_col_types = FALSE)
stopifnot(nrow(ring_table) == 76)

binary_check_cols <- c("GH190", "GT121", "GT26",
                        "antibiotic_degradation", "biofilm", "Extracellular",
                        "quorum_sensing", "secretion_system", "siderophores")
for (col in binary_check_cols) {
  vals <- unique(ring_table[[col]])
  if (!all(vals %in% c(0, 1))) {
    stop(sprintf(
      "Column '%s' contains non-0/1 values (%s) — rerun rebuild_ring_annotation_lineage76.R with the binary-encoding fix before generating iTOL files.",
      col, paste(vals, collapse = ", ")
    ))
  }
}
cat("Binary column type check passed — all indicator columns are clean 0/1.\n\n")

ring_table <- ring_table %>%
  mutate(host = str_split_fixed(genome, "__", 4)[, 2])

cat("Host species distribution:\n")
print(count(ring_table, host, sort = TRUE))

hosts <- sort(unique(ring_table$host))
host_palette <- setNames(
  RColorBrewer::brewer.pal(max(3, length(hosts)), "Set3")[seq_along(hosts)],
  hosts
)

branch_color_lines <- c(
  "TREE_COLORS",
  "SEPARATOR TAB",
  "DATA"
)
branch_color_lines <- c(
  branch_color_lines,
  paste(ring_table$genome, "branch", host_palette[ring_table$host], "normal", "1", sep = "\t")
)
writeLines(branch_color_lines, file.path(out_dir, "01_host_species_branch_colors.txt"))

write_csv(
  tibble(host = names(host_palette), color = host_palette),
  file.path(out_dir, "00_host_color_legend_REFERENCE_ONLY.csv")
)

genera <- sort(unique(ring_table$genus))
genus_palette <- setNames(
  RColorBrewer::brewer.pal(max(3, length(genera)), "Paired")[seq_along(genera)],
  genera
)

colorstrip_header <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tGenus",
  "COLOR\t#000000",
  "COLOR_BRANCHES\t0",
  "LEGEND_TITLE\tGenus",
  paste0("LEGEND_SHAPES\t", paste(rep(1, length(genera)), collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(genus_palette, collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(genera, collapse = "\t")),
  "DATA"
)
colorstrip_data <- paste(ring_table$genome, genus_palette[ring_table$genus], ring_table$genus, sep = "\t")
writeLines(c(colorstrip_header, colorstrip_data), file.path(out_dir, "02_genus_colorstrip.txt"))

cazy_families <- c("GH190", "GT121", "GT26")
cazy_colors <- c(GH190 = "#4DAF4A", GT121 = "#984EA3", GT26 = "#FF7F00")

cazy_header <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tCAZyme headline families",
  "COLOR\t#000000",
  paste0("FIELD_SHAPES\t", paste(rep(2, length(cazy_families)), collapse = "\t")),
  paste0("FIELD_LABELS\t", paste(cazy_families, collapse = "\t")),
  paste0("FIELD_COLORS\t", paste(cazy_colors[cazy_families], collapse = "\t")),
  "LEGEND_TITLE\tCAZyme family",
  paste0("LEGEND_SHAPES\t", paste(rep(2, length(cazy_families)), collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(cazy_colors[cazy_families], collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(cazy_families, collapse = "\t")),
  "DATA"
)
cazy_data <- ring_table %>%
  select(genome, all_of(cazy_families)) %>%
  unite("vals", all_of(cazy_families), sep = "\t") %>%
  mutate(line = paste(genome, vals, sep = "\t")) %>%
  pull(line)
writeLines(c(cazy_header, cazy_data), file.path(out_dir, "03_cazyme_binary.txt"))

coop_categories <- c(
  "antibiotic_degradation", "biofilm", "Extracellular",
  "quorum_sensing", "secretion_system", "siderophores"
)
coop_colors <- setNames(
  RColorBrewer::brewer.pal(6, "Dark2"),
  coop_categories
)

coop_header <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tCooperative-trait categories",
  "COLOR\t#000000",
  paste0("FIELD_SHAPES\t", paste(rep(2, length(coop_categories)), collapse = "\t")),
  paste0("FIELD_LABELS\t", paste(coop_categories, collapse = "\t")),
  paste0("FIELD_COLORS\t", paste(coop_colors[coop_categories], collapse = "\t")),
  "LEGEND_TITLE\tCooperative trait",
  paste0("LEGEND_SHAPES\t", paste(rep(2, length(coop_categories)), collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(coop_colors[coop_categories], collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(coop_categories, collapse = "\t")),
  "DATA"
)
coop_data <- ring_table %>%
  select(genome, all_of(coop_categories)) %>%
  unite("vals", all_of(coop_categories), sep = "\t") %>%
  mutate(line = paste(genome, vals, sep = "\t")) %>%
  pull(line)
writeLines(c(coop_header, coop_data), file.path(out_dir, "04_cooperative_binary.txt"))

cat("\nAll iTOL dataset files written to:", out_dir, "\n")
cat("Files:\n")
print(list.files(out_dir))
})
