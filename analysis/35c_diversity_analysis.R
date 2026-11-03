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
#  SCRIPT  : 35c_diversity_analysis.R  -  Alpha and beta diversity by host
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv
#            phase_4_SV_analysis/working/samples_manifest.tsv
#  OUTPUT  : annotation/results_phase3/s35_diversity/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })

if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan", repos = "https://cloud.r-project.org")
library(vegan)

ABUND_F <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv")
MANIFEST_F <- paste0(URGG_PROJECT, "/phase_4_SV_analysis/working/samples_manifest.tsv")
OUT_DIR <- paste0(URGG_PROJECT, "/annotation/results_phase3/s35_diversity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(ABUND_F, MANIFEST_F)) if (!file.exists(f)) stop("Missing: ", f)

message("Reading abundance matrix ...")
abund <- read_tsv(ABUND_F, show_col_types = FALSE)
abund <- abund %>% filter(Genome != "unmapped")

mat <- as.matrix(abund[, -1])
rownames(mat) <- abund$Genome
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- abund$Genome

cat(sprintf("Genomes: %d | Samples: %d\n", nrow(mat), ncol(mat)))

message("Reading manifest for host/species metadata ...")
manifest <- read_tsv(MANIFEST_F, show_col_types = FALSE) %>%
  distinct(sample_id, bioproject, species)

message("Computing alpha diversity ...")
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
  warning("Fewer than 90% of abundance-matrix samples matched manifest metadata — ",
          "check that CoverM column renaming in 35b succeeded (see its ",
          "'COULD NOT MATCH' warnings).")
}

write_csv(alpha_df, file.path(OUT_DIR, "alpha_diversity_per_sample.csv"))

cat("\n=== Alpha diversity summary by host species ===\n")
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

message("Computing beta diversity ...")
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

message("Running PERMANOVA (999 permutations) ...")
meta_ordered <- manifest %>% filter(sample_id %in% rownames(sample_mat)) %>%
  arrange(match(sample_id, rownames(sample_mat)))
stopifnot(identical(meta_ordered$sample_id, rownames(sample_mat)))

permanova <- adonis2(bc_dist ~ species, data = meta_ordered, permutations = 999)
cat("\n=== PERMANOVA: host species effect on community composition ===\n")
print(permanova)
write_csv(as.data.frame(permanova) %>% tibble::rownames_to_column("term"),
         file.path(OUT_DIR, "permanova_host_species.csv"))

cat("\nAll outputs written to:", OUT_DIR, "\n")
