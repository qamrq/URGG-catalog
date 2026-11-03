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
#  SCRIPT  : 08_generate_itol_datasets_lineage76.R  -  iTOL datasets for the RUC tree (Fig. 5a)
#  AUTHOR  : QRQ
#  INPUT   : annotation/iqtree_output/tree_build_lineage_only_fresh/tree_ring_annotation_lineage76.csv
#  OUTPUT  : annotation/iqtree_output/tree_build_lineage_only_fresh/itol_datasets/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

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
