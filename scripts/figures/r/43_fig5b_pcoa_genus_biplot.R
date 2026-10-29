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
#  SCRIPT  : 43_fig5b_pcoa_genus_biplot.R  -  Fig. 5b PCoA with genus vectors (envfit biplot)
#  AUTHOR  : QRQ
#  INPUT   : Abundance matrix (MATRIX_F), sample manifest (MANIFEST_F), GTDB summary (GTDB_SUMMARY)
#  OUTPUT  : Fig. 5b biplot (PDF, PNG) in OUT_DIR
# ==========================================================

CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/figures_script"
source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(vegan); library(cowplot)
})

MATRIX_F   <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv"
MANIFEST_F <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_4_SV_analysis/working/samples_manifest.tsv"
GTDB_SUMMARY <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv"
OUT_DIR    <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity"

species_palette <- c(
  "Ovis aries" = "#2F5496", "Pantholops hodgsonii" = "#4776B3", "Bos taurus" = "#5C8FD1",
  "Bos grunniens" = "#7DA8E6", "Bos frontalis" = "#9FBEEF", "Capra hircus" = "#3C8F8A",
  "Bubalus bubalis" = "#56B4B0", "Procapra picticaudata" = "#8FBF8A",
  "Camelus bactrianus" = "#BF8F00", "Camelus dromedarius" = "#D9AE4E",
  "Alces alces" = "#3C7A3C", "Alces americanus" = "#6AA06A",
  "Hydropotes inermis" = "#8FBF6A", "Capreolus pygargus" = "#B5D98A",
  "Giraffa camelopardalis" = "#8B5A9E"
)

message("Reading abundance matrix ...")
abund <- read_tsv(MATRIX_F, show_col_types = FALSE) %>% filter(Genome != "unmapped")
mat <- as.matrix(abund[, -1]); rownames(mat) <- abund$Genome
mat <- apply(mat, 2, as.numeric); rownames(mat) <- abund$Genome
sample_mat <- t(mat)

manifest <- read_tsv(MANIFEST_F, show_col_types = FALSE) %>% distinct(sample_id, species)
species_lookup <- setNames(gsub("_", " ", manifest$species), manifest$sample_id)

meta <- tibble(sample = rownames(sample_mat), species = species_lookup[rownames(sample_mat)]) %>%
  filter(!is.na(species))
sample_mat <- sample_mat[meta$sample, , drop = FALSE]

message("Computing Bray-Curtis PCoA ...")
bc_dist <- vegdist(sample_mat, method = "bray")
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
var_explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

message("Aggregating to genus level and fitting envfit vectors ...")
gtdb <- read_tsv(GTDB_SUMMARY, show_col_types = FALSE)
extract_genus <- function(classification) {
  m <- str_extract(classification, "g__[^;]*")
  ifelse(is.na(m) | m == "g__", NA_character_, sub("^g__", "", m))
}
genome_genus <- gtdb %>% transmute(genome = user_genome, genus = extract_genus(classification)) %>%
  filter(!is.na(genus))

genus_mat_long <- as.data.frame(sample_mat) %>%
  mutate(sample = rownames(sample_mat)) %>%
  pivot_longer(-sample, names_to = "genome", values_to = "relabund") %>%
  inner_join(genome_genus, by = "genome") %>%
  group_by(sample, genus) %>%
  summarise(relabund = sum(relabund), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = relabund, values_fill = 0)

genus_mat <- as.matrix(genus_mat_long[, -1])
rownames(genus_mat) <- genus_mat_long$sample
genus_mat <- genus_mat[rownames(sample_mat), , drop = FALSE]

top_genera <- names(sort(colMeans(genus_mat), decreasing = TRUE))[1:15]

set.seed(42)
fit <- envfit(pcoa$points, genus_mat[, top_genera], permutations = 999)

vectors <- as.data.frame(scores(fit, "vectors")) %>%
  mutate(genus = rownames(.), r2 = fit$vectors$r, p = fit$vectors$pvals) %>%
  filter(p < 0.05)

cat(sprintf("Genera with significant vectors (p<0.05) out of top 15 tested: %d\n", nrow(vectors)))
print(as.data.frame(vectors %>% select(genus, r2, p) %>% arrange(desc(r2))))

arrow_mult <- 0.9 * max(abs(pcoa$points)) / max(abs(as.matrix(vectors[, 1:2])))
vectors <- vectors %>% mutate(Dim1 = Dim1 * arrow_mult, Dim2 = Dim2 * arrow_mult)

pcoa_df <- tibble(sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2]) %>%
  left_join(meta, by = "sample")

p_biplot <- ggplot(pcoa_df, aes(PC1, PC2, color = species)) +
  geom_point(size = 1.1, alpha = 0.5) +
  stat_ellipse(aes(group = species), linewidth = 0.3, alpha = 0.4, level = 0.68) +
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
              inherit.aes = FALSE, color = "black", linewidth = 0.4,
              arrow = arrow(length = unit(0.15, "cm"))) +
  geom_text(data = vectors, aes(Dim1, Dim2, label = genus),
            inherit.aes = FALSE, size = 2.6, fontface = "italic", color = "black",
            nudge_x = 0.02, nudge_y = 0.02) +
  scale_color_manual(values = species_palette, name = NULL, guide = "none") +
  labs(title = "Genus-vector biplot on the host-species Bray-Curtis ordination",
       x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
       y = sprintf("PCoA2 (%.1f%%)", var_explained[2]),
       subtitle = "Arrows: significant genus abundance vectors (envfit, p<0.05, 999 permutations)") +
  theme_catalogue_cartesian(base_size = 9)

export_catalogue_plot(p_biplot, OUT_DIR, "Fig5b_pcoa_genus_biplot", width = 8, height = 7)

cat("\nWritten to:", OUT_DIR, "\n")
cat("This can either REPLACE Fig 5b in the main text (if the vectors add clear\n")
cat("value without overcrowding) or serve as a supplementary extension of it --\n")
cat("recommend checking the rendered image before deciding which.\n")
