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
#  SCRIPT  : 36_fig5_diversity_composite.R  -  Fig. 5 diversity composite: Shannon boxplot and PCoA
#  AUTHOR  : QRQ
#  INPUT   : Final abundance matrix (MATRIX_F), sample manifest (MANIFEST_F)
#  OUTPUT  : Fig. 5a, 5b composite (PDF, PNG)
# ==========================================================

CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/figures_script"
source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(forcats)
  library(vegan); library(cowplot)
})

MATRIX_F   <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv"
MANIFEST_F <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_4_SV_analysis/working/samples_manifest.tsv"
OUT_DIR    <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity"

# =============================================================================
# Hand-curated, family-grouped palette (dark -> light shading within each
# host family, never auto-generated) -- extends the palette pattern already
# established for this manuscript's other figures.
# =============================================================================
species_palette <- c(
  "Ovis aries"             = "#2F5496",
  "Pantholops hodgsonii"   = "#4776B3",
  "Bos taurus"              = "#5C8FD1",
  "Bos grunniens"           = "#7DA8E6",
  "Bos frontalis"           = "#9FBEEF",
  "Capra hircus"            = "#3C8F8A",
  "Bubalus bubalis"         = "#56B4B0",
  "Procapra picticaudata"   = "#8FBF8A",
  "Camelus bactrianus"      = "#BF8F00",
  "Camelus dromedarius"     = "#D9AE4E",
  "Alces alces"             = "#3C7A3C",
  "Alces americanus"        = "#6AA06A",
  "Hydropotes inermis"      = "#8FBF6A",
  "Capreolus pygargus"      = "#B5D98A",
  "Giraffa camelopardalis"  = "#8B5A9E"
)

message("Reading abundance matrix ...")
abund <- read_tsv(MATRIX_F, show_col_types = FALSE) %>% filter(Genome != "unmapped")
mat <- as.matrix(abund[, -1]); rownames(mat) <- abund$Genome
mat <- apply(mat, 2, as.numeric); rownames(mat) <- abund$Genome
sample_mat <- t(mat)
cat(sprintf("Genomes: %d | Samples: %d\n", ncol(sample_mat), nrow(sample_mat)))

manifest <- read_tsv(MANIFEST_F, show_col_types = FALSE) %>% distinct(sample_id, species)
species_lookup <- setNames(gsub("_", " ", manifest$species), manifest$sample_id)

shannon <- diversity(sample_mat, index = "shannon")
richness <- specnumber(sample_mat)

div_df <- tibble(
  sample = rownames(sample_mat),
  shannon = shannon,
  richness = richness,
  species = species_lookup[rownames(sample_mat)]
) %>% filter(!is.na(species))

cat(sprintf("Samples with diversity computed: %d\n", nrow(div_df)))

species_order <- div_df %>%
  group_by(species) %>%
  summarise(med = median(shannon), n = n(), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(species)

div_df <- div_df %>% mutate(species = factor(species, levels = species_order))
n_labels <- div_df %>% group_by(species) %>% summarise(n = n(), .groups = "drop")

message("Building panel A (Shannon boxplot) ...")
pA <- ggplot(div_df, aes(species, shannon, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.25, color = "grey20") +
  geom_text(data = n_labels, aes(x = species, y = -0.3, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.3, color = "grey40", vjust = 1) +
  scale_fill_manual(values = species_palette, guide = "none") +
  labs(title = "a", x = NULL, y = "Shannon diversity index") +
  coord_cartesian(ylim = c(-0.6, NA), clip = "off") +
  theme_catalogue_cartesian(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        plot.margin = margin(5, 5, 20, 5))

message("Computing Bray-Curtis distance and PCoA (may take a few minutes for 1,369 samples) ...")
bc_dist <- vegdist(sample_mat, method = "bray")
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
var_explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

pcoa_df <- tibble(
  sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(div_df %>% select(sample, species), by = "sample") %>%
  filter(!is.na(species))

cat(sprintf("PCoA variance explained: PC1=%.1f%%, PC2=%.1f%%\n", var_explained[1], var_explained[2]))

pB <- ggplot(pcoa_df, aes(PC1, PC2, color = species)) +
  geom_point(size = 1.1, alpha = 0.6) +
  stat_ellipse(aes(group = species), linewidth = 0.3, alpha = 0.5, level = 0.68) +
  scale_color_manual(values = species_palette, name = NULL,
                      labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(title = "b",
       x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
       y = sprintf("PCoA2 (%.1f%%)", var_explained[2])) +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

top_title <- ggdraw() + draw_label(
  "Figure 5. Community diversity is significantly structured by host species",
  fontface = "bold", size = 12, x = .5, hjust = .5)

final_plot <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.15))
final_with_title <- plot_grid(top_title, final_plot, ncol = 1, rel_heights = c(.08, 1))

export_catalogue_plot(final_with_title, OUT_DIR, "Fig5_diversity_composite", width = 13, height = 6)

cat("\nWritten to:", OUT_DIR, "\n")
cat("Panel (a): Shannon diversity by host, ordered by descending median (matches Table 8)\n")
cat("Panel (b): Bray-Curtis PCoA, 68% confidence ellipses per host\n")
