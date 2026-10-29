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
#  SCRIPT  : 44_cazyme_consolidated_figure.R  -  Consolidated CAZyme figure
#  AUTHOR  : QRQ
#  INPUT   : table_cazyme_consolidated.csv (CONSOLIDATED_F)
#  OUTPUT  : CAZyme figure (PDF, PNG) in OUT_DIR
# ==========================================================

CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/figures_script"
source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({ library(dplyr); library(readr); library(cowplot); library(forcats) })

DATA_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s34_cazyme_core_specific"
CONSOLIDATED_F <- file.path(DATA_DIR, "table_cazyme_consolidated.csv")
OUT_DIR <- DATA_DIR

if (!file.exists(CONSOLIDATED_F)) {
  stop("Missing: ", CONSOLIDATED_F, "\nRun 40_cazyme_consolidated_table.R first.")
}

message("Reading consolidated table ...")
dat <- read_csv(CONSOLIDATED_F, show_col_types = FALSE)

substrate_palette <- c(
  "alpha-glucan" = "#2F5496",
  "peptidoglycan" = "#3C7A3C",
  "chitin" = "#BF8F00",
  "host glycan" = "#8B5A9E",
  "bacterial outer membrane (OM) polysaccharides" = "#D62728",
  "unresolved (no dbCAN_sub substrate call)" = "grey60"
)
dat <- dat %>% mutate(substrate = ifelse(is.na(substrate), "unresolved (no dbCAN_sub substrate call)", substrate))

missing_from_palette <- setdiff(unique(dat$substrate), names(substrate_palette))
if (length(missing_from_palette) > 0) {
  stop("Substrate(s) found with no assigned color: ", paste(missing_from_palette, collapse=", "),
       "\nAdd these to substrate_palette above before proceeding.")
}

family_order <- dat %>%
  arrange(substrate, desc(log2_OR_catalogue)) %>%
  pull(family)
dat <- dat %>% mutate(family = factor(family, levels = rev(family_order)))

message("Building panel A (catalogue-wide) ...")
pA <- ggplot(dat, aes(log2_OR_catalogue, family, color = substrate)) +
  geom_vline(xintercept = 0, linewidth = .4, color = "grey30") +
  geom_errorbarh(aes(xmin = log2(OR_catalogue_ci_low), xmax = log2(OR_catalogue_ci_high)),
                height = .25, linewidth = .5) +
  geom_point(size = 2.5) +
  geom_text(aes(label = ifelse(!is.na(p_adj_catalogue),
                               sprintf("p=%.1e", p_adj_catalogue), "")),
            hjust = -0.15, size = 2.3, color = "grey20", show.legend = FALSE) +
  scale_color_manual(values = substrate_palette, name = "Substrate") +
  labs(title = "a", x = expression(log[2]~"(OR, core/broad vs. host-specific)"), y = NULL,
       subtitle = "Catalogue-wide tier association (all 9 families)") +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "none")

message("Building panel B (lineage-specific) ...")
dat_lineage <- dat %>% filter(has_lineage_data)

pB <- ggplot(dat, aes(log2_OR_lineage, family, color = substrate)) +
  geom_vline(xintercept = 0, linewidth = .4, color = "grey30") +
  geom_errorbarh(data = dat_lineage,
                 aes(xmin = log2(OR_lineage_ci_low), xmax = log2(OR_lineage_ci_high)),
                height = .25, linewidth = .5) +
  geom_point(data = dat_lineage, size = 2.5) +
  geom_text(data = dat_lineage,
            aes(label = sprintf("p=%.1e", p_adj_lineage)),
            hjust = -0.15, size = 2.3, color = "grey20", show.legend = FALSE) +
  scale_color_manual(values = substrate_palette, name = "Substrate") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "b", x = expression(log[2]~"(OR, UMGS1840/HGM11507 lineage vs. rest of catalogue)"), y = NULL,
       subtitle = "Lineage-specific association (GT121/GT26/GH190 only -- tested for these 3 families)") +
  theme_catalogue_cartesian(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_blank())

top_row <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.15, 1))

title <- ggdraw() + draw_label(
  "Figure 4. CAZyme family enrichment, substrate, and lineage-specific association in one view",
  fontface = "bold", size = 12, x = .5, hjust = .5)

final_with_title <- plot_grid(title, top_row, ncol = 1, rel_heights = c(.10, 1))

export_catalogue_plot(final_with_title, OUT_DIR, "Fig4_cazyme_consolidated", width = 13, height = 6)

cat("\nWritten to:", OUT_DIR, "\n")
cat("Panel a: all 9 families' catalogue-wide tier enrichment\n")
cat("Panel b: the 3 families (GT121/GT26/GH190) with lineage-specific data --\n")
cat("         other 6 rows correctly show no point, not missing data\n")
cat("Both panels share the same family ordering (by substrate, then effect size)\n")
cat("and color scheme, so a reader can directly compare panel a vs b per family.\n")
