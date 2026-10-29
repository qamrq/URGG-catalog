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
#  SCRIPT  : 42_supp_genus_abundance_overview.R  -  Supplementary genus abundance overview
#  AUTHOR  : QRQ
#  INPUT   : Genus abundance table (GENUS_F)
#  OUTPUT  : Supplementary figure (PDF, PNG) in OUT_DIR
# ==========================================================

FIG_R_DIR <- "scripts/figures/r"
source(file.path(FIG_R_DIR, "00_install_packages.R"))
source(file.path(FIG_R_DIR, "01_catalogue_theme.R"))
source(file.path(FIG_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({ library(dplyr); library(readr); library(cowplot); library(forcats) })

GENUS_F <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity/table_genus_compositional_overview.csv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s35_diversity"

PREVALENCE_FLOOR <- 20  # same floor already applied for the manuscript's real numbers

genus_summary <- read_csv(GENUS_F, show_col_types = FALSE) %>%
  filter(n_samples_present >= PREVALENCE_FLOOR) %>%
  arrange(desc(median_relabund)) %>%
  slice_head(n = 20) %>%
  mutate(genus = fct_reorder(genus, median_relabund),
        rank_group = ifelse(row_number() <= 10, "Top 10", "Rank 11-20"))

# =============================================================================
# PANEL A -- top 20 by median abundance, top 10 highlighted
# =============================================================================
pA <- ggplot(genus_summary, aes(median_relabund, genus, fill = rank_group)) +
  geom_col(width = .7) +
  scale_fill_manual(values = c("Top 10" = "#D62728", "Rank 11-20" = "grey65"), name = NULL) +
  labs(title = "a", x = "Median relative abundance (%)", y = NULL,
       subtitle = sprintf("Top 20 genera (prevalence \u2265%d samples); top 10 highlighted", PREVALENCE_FLOOR)) +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))

# =============================================================================
# PANEL B -- same genera, colored by CV (stability vs. variability)
# =============================================================================
pB <- ggplot(genus_summary, aes(median_relabund, genus, color = cv_relabund, size = prevalence_pct)) +
  geom_point() +
  scale_color_gradient(low = "#2F5496", high = "#D62728", name = "CV") +
  scale_size_continuous(name = "Prevalence (%)", range = c(2, 6)) +
  labs(title = "b", x = "Median relative abundance (%)", y = NULL,
       subtitle = "Same genera, coloured by coefficient of variation") +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))

top_title <- ggdraw() + draw_label(
  "Figure S_. Genus-level relative abundance overview across the full abundance-profiled catalogue",
  fontface = "bold", size = 11, x = .5, hjust = .5)

final_plot <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.1))
final_with_title <- plot_grid(top_title, final_plot, ncol = 1, rel_heights = c(.08, 1))

export_catalogue_plot(final_with_title, OUT_DIR, "FigS_genus_abundance_overview", width = 12, height = 6)

cat("Written to:", OUT_DIR, "\n")
