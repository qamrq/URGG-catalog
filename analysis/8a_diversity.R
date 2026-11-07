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
#  SCRIPT  : 8a_diversity.R  -  Alpha and beta diversity, PERMANOVA and genus envfit (Fig. 2b)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv
#            phase_4_SV_analysis/working/samples_manifest.tsv
#            annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#  OUTPUT  : annotation/results_phase3/s35_diversity/
#  USAGE   : Rscript 8a_diversity.R [alpha_beta pcoa_envfit]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("alpha_beta")) local({
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })

library(vegan)

ABUND_F <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv")
MANIFEST_F <- paste0(URGG_PROJECT, "/phase_4_SV_analysis/working/samples_manifest.tsv")
OUT_DIR <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(ABUND_F, MANIFEST_F)) if (!file.exists(f)) stop("Missing: ", f)

abund <- read_tsv(ABUND_F, show_col_types = FALSE)
abund <- abund %>% filter(Genome != "unmapped")

mat <- as.matrix(abund[, -1])
rownames(mat) <- abund$Genome
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- abund$Genome

cat(sprintf("Genomes: %d | Samples: %d\n", nrow(mat), ncol(mat)))

manifest <- read_tsv(MANIFEST_F, show_col_types = FALSE) %>%
  distinct(sample_id, bioproject, species)

sample_mat <- t(mat)

shannon <- diversity(sample_mat, index = "shannon")
richness <- specnumber(sample_mat)

alpha_df <- tibble(sample_id = rownames(sample_mat),
                   shannon = shannon, richness = richness) %>%
  inner_join(manifest, by = "sample_id")

n_matched <- nrow(alpha_df)
n_samples_in_matrix <- ncol(mat)
cat(sprintf("Samples matched to host metadata: %d of %d in abundance matrix\n",
            n_matched, n_samples_in_matrix))
if (n_matched < 0.9 * n_samples_in_matrix) {
  warning("Fewer than 90% of abundance-matrix samples matched manifest metadata - ",
          "check that CoverM column renaming in 35b succeeded (see its ",
          "'COULD NOT MATCH' warnings).")
}

write_csv(alpha_df, file.path(OUT_DIR, "alpha_diversity_per_sample.csv"))

alpha_summary <- alpha_df %>%
  group_by(species) %>%
  summarise(n_samples = n(),
            median_shannon = round(median(shannon), 3),
            median_richness = round(median(richness), 1), .groups = "drop") %>%
  arrange(desc(median_shannon))
print(as.data.frame(alpha_summary))
write_csv(alpha_summary, file.path(OUT_DIR, "alpha_diversity_by_host_summary.csv"))

ALL_CATALOGUE_SPECIES <- c("Ovis_aries","Pantholops_hodgsonii","Bos_taurus","Bos_grunniens",
  "Capra_hircus","Bubalus_bubalis","Bos_indicus","Hydropotes_inermis","Alces_alces",
  "Procapra_picticaudata","Cervus_elaphus","Rangifer_tarandus","Alces_americanus",
  "Giraffa_camelopardalis","Camelus_bactrianus","Bos_frontalis","Camelus_dromedarius")
species_absent <- setdiff(ALL_CATALOGUE_SPECIES, unique(alpha_df$species))
cat(sprintf("\nCatalogue host species with NO read-mapping coverage in this analysis: %d of %d\n",
            length(species_absent), length(ALL_CATALOGUE_SPECIES)))
if (length(species_absent) > 0) cat("  ", paste(species_absent, collapse = ", "), "\n")

if (n_distinct(alpha_df$species) > 2) {
  kw_shannon <- kruskal.test(shannon ~ species, data = alpha_df)
  cat(sprintf("\nKruskal-Wallis (Shannon ~ host species): chi-sq = %.2f, p = %s\n",
              kw_shannon$statistic, format.pval(kw_shannon$p.value)))
}

p_alpha <- ggplot(alpha_df, aes(reorder(species, shannon, median), shannon)) +
  geom_boxplot(fill = "#2F5496", alpha = .6, outlier.size = .5) +
  coord_flip() +
  labs(title = "Alpha diversity (Shannon index) by host species",
       x = NULL, y = "Shannon index") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic"))
ggsave(file.path(OUT_DIR, "Fig_alpha_diversity_by_host.pdf"), p_alpha, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "Fig_alpha_diversity_by_host.png"), p_alpha, width = 7, height = 6, dpi = 300)

bc_dist <- vegdist(sample_mat, method = "bray")

pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
var_explained <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

pcoa_df <- tibble(sample_id = rownames(sample_mat),
                  PCo1 = pcoa$points[,1], PCo2 = pcoa$points[,2]) %>%
  inner_join(manifest, by = "sample_id")

write_csv(pcoa_df, file.path(OUT_DIR, "beta_diversity_pcoa_coordinates.csv"))

p_beta <- ggplot(pcoa_df, aes(PCo1, PCo2, color = species)) +
  geom_point(size = 2, alpha = .7) +
  labs(title = "Beta diversity (Bray-Curtis PCoA) by host species",
       x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
       y = sprintf("PCo2 (%.1f%%)", var_explained[2]),
       color = "Host species") +
  theme_minimal(base_size = 10) +
  theme(legend.text = element_text(face = "italic", size = 7))
ggsave(file.path(OUT_DIR, "Fig_beta_diversity_pcoa.pdf"), p_beta, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "Fig_beta_diversity_pcoa.png"), p_beta, width = 8, height = 6, dpi = 300)

meta_ordered <- manifest %>% filter(sample_id %in% rownames(sample_mat)) %>%
  arrange(match(sample_id, rownames(sample_mat)))
stopifnot(identical(meta_ordered$sample_id, rownames(sample_mat)))

permanova <- adonis2(bc_dist ~ species, data = meta_ordered, permutations = 999)
print(permanova)
write_csv(as.data.frame(permanova) %>% tibble::rownames_to_column("term"),
         file.path(OUT_DIR, "permanova_host_species.csv"))

cat("\nAll outputs written to:", OUT_DIR, "\n")
})

if (run_part("pcoa_envfit")) local({
CATALOGUE_R_DIR <- Sys.getenv("CATALOGUE_R_DIR", "analysis/helpers")
source(file.path(CATALOGUE_R_DIR, "catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "catalogue_export.R"))

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
})
