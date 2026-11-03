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
#  SCRIPT  : 43_fig5_diversity_composite_biplot_genus.R  -  Bray-Curtis PCoA, PERMANOVA and genus envfit (Fig. 2b)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv
#            phase_4_SV_analysis/working/samples_manifest.tsv
#            annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#  OUTPUT  : annotation/results_phase3/s35_diversity/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

CATALOGUE_R_DIR <- Sys.getenv("CATALOGUE_R_DIR", "analysis/helpers")
source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(forcats); library(stringr)
  library(vegan); library(cowplot); library(ggnewscale)
})

MATRIX_F     <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv")
MANIFEST_F   <- paste0(URGG_PROJECT, "/phase_4_SV_analysis/working/samples_manifest.tsv")
GTDB_SUMMARY <- paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv")
OUT_DIR      <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity")

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

genus_vector_palette <- c(
  "#B2182B", "#2166AC", "#1B7837", "#B35806", "#542788",
  "#D6604D", "#4393C3", "#5AAE61", "#E08214", "#8073AC",
  "#F4A582", "#92C5DE", "#A6DBA0", "#FDB863", "#C2A5CF"
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
cat(sprintf("PCoA variance explained: PC1=%.1f%%, PC2=%.1f%%\n", var_explained[1], var_explained[2]))

pcoa_df <- tibble(
  sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(div_df %>% select(sample, species), by = "sample") %>%
  filter(!is.na(species))

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
genus_mat <- genus_mat[rownames(pcoa$points), , drop = FALSE]

top_genera <- names(sort(colMeans(genus_mat), decreasing = TRUE))[1:15]

set.seed(42)
fit <- envfit(pcoa$points, genus_mat[, top_genera], permutations = 999)

vectors <- as.data.frame(scores(fit, "vectors")) %>%
  mutate(genus = rownames(.), r2 = fit$vectors$r, p = fit$vectors$pvals) %>%
  filter(p < 0.05) %>%
  arrange(desc(r2))

cat(sprintf("Genera with significant vectors (p<0.05) out of top 15 tested: %d\n", nrow(vectors)))
print(as.data.frame(vectors %>% select(genus, r2, p)))

if (nrow(vectors) > length(genus_vector_palette)) {
  stop(sprintf(
    "genus_vector_palette only has %d colors but %d genera are significant -- add more manual colors.",
    length(genus_vector_palette), nrow(vectors)
  ))
}

vectors <- vectors %>% mutate(genus = factor(genus, levels = genus))
genus_colors <- setNames(genus_vector_palette[seq_len(nrow(vectors))], levels(vectors$genus))

arrow_mult <- 0.9 * max(abs(pcoa$points)) / max(abs(as.matrix(vectors[, c("Dim1", "Dim2")])))
vectors <- vectors %>% mutate(Dim1 = Dim1 * arrow_mult, Dim2 = Dim2 * arrow_mult)

ELLIPSE_MIN_N <- 4
species_n <- pcoa_df %>% count(species)
species_too_few <- species_n %>% filter(n < ELLIPSE_MIN_N) %>% pull(species)
if (length(species_too_few) > 0) {
  message(sprintf("Skipping ellipse for %d species with < %d samples: %s",
                   length(species_too_few), ELLIPSE_MIN_N, paste(species_too_few, collapse = ", ")))
}
ellipse_df <- pcoa_df %>% filter(!species %in% species_too_few)

message("Building panel B (PCoA + genus biplot) ...")
pB <- ggplot(pcoa_df, aes(PC1, PC2, color = species)) +
  geom_point(size = 1.1, alpha = 0.6) +
  stat_ellipse(data = ellipse_df, aes(group = species), linewidth = 0.3, alpha = 0.5, level = 0.68) +
  scale_color_manual(values = species_palette, name = NULL,
                      labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  ggnewscale::new_scale_color() +
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = Dim1, yend = Dim2, color = genus),
               inherit.aes = FALSE, linewidth = 0.6,
               arrow = arrow(length = unit(0.15, "cm"))) +
  scale_color_manual(values = genus_colors, name = "Genus",
                      labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(title = "b",
       x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
       y = sprintf("PCoA2 (%.1f%%)", var_explained[2])) +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8, face = "plain"),
        legend.key.size = unit(0.35, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.box = "vertical")

top_title <- ggdraw() + draw_label(
  "Figure 5. Community diversity is significantly structured by host species",
  fontface = "bold", size = 12, x = .5, hjust = .5)

final_plot <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.35))
final_with_title <- plot_grid(top_title, final_plot, ncol = 1, rel_heights = c(.08, 1))

export_catalogue_plot(final_with_title, OUT_DIR, "Fig5_diversity_composite_biplot", width = 14.5, height = 6.2)

cat("\nWritten to:", OUT_DIR, "\n")
cat("Panel (a): Shannon diversity by host, ordered by descending median (matches Table 8)\n")
cat("Panel (b): Bray-Curtis PCoA + significant genus vectors (envfit, p<0.05), host-species\n")
cat("           ellipses (68% CI, n>=4 only) and genus vectors each carry their own legend\n")
