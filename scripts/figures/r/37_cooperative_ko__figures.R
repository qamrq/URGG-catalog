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
#  SCRIPT  : 37_cooperative_ko__figures.R  -  Cooperative-trait figure (Fig. 7)
#  AUTHOR  : QRQ
#  INPUT   : Proportion and category tables from 47 and 48 (PROPORTION_F, CATEGORY_F)
#  OUTPUT  : Fig. 7 (PDF) and table_cooperative_category_stats_R.csv in OUT_DIR
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork)
})

PROPORTION_F <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_proportion_per_genome.csv"
CATEGORY_F   <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_category_by_tier.csv"
OUT_DIR      <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s37_cooperative_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(PROPORTION_F, CATEGORY_F)) if (!file.exists(f)) stop("Required input missing: ", f)

TIER_LEVELS <- c("Host-specific (1 host species)", "Shared (2-3 host species)", "Core / broad (>3 host species)")
TIER_LABELS <- c("Host-specific", "Shared", "Core/broad")

# House palette: core/broad matches the table-header blue already used
# elsewhere in the manuscript (#2F5496); host-specific gets a distinct warm
# accent so the two enrichment directions read apart at a glance; shared and
# non-significant results are muted grey rather than a third saturated hue.
PAL <- c(
  "Core/broad"        = "#2F5496",
  "Shared"             = "#8FA9CE",
  "Host-specific"      = "#C0504D",
  "Core-enriched"      = "#2F5496",
  "Host-specific-enriched" = "#C0504D",
  "Not significant"    = "#B0B0B0"
)

theme_ms <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(color = "grey30", linewidth = 0.3),
      axis.ticks = element_line(color = "grey30", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

# =============================================================================
# PANEL A — aggregate proportion by tier
# =============================================================================
prop <- read_csv(PROPORTION_F, show_col_types = FALSE) %>%
  mutate(sharing_category = factor(sharing_category, levels = TIER_LEVELS, labels = TIER_LABELS))

panel_a <- ggplot(prop, aes(x = sharing_category, y = 100 * proportion_cooperative, fill = sharing_category)) +
  geom_violin(alpha = 0.5, color = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "grey20", fill = "white") +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = "Cooperative genes\n(% of genome, all 6 categories pooled)",
       title = "a") +
  theme_ms()

# =============================================================================
# STATS — per-category chi-square omnibus + Fisher OR/CI (core vs. host-specific)
# =============================================================================
cat_raw <- read_csv(CATEGORY_F, show_col_types = FALSE)

test_category <- function(df) {
  tiers <- df %>% select(sharing_category, n, n_genomes_total) %>%
    mutate(n_without = n_genomes_total - n)
  mat <- as.matrix(tiers %>% select(n, n_without))
  rownames(mat) <- tiers$sharing_category
  chi <- suppressWarnings(chisq.test(mat))
  chi_p <- if (any(chi$expected < 5)) {
    suppressWarnings(chisq.test(mat, simulate.p.value = TRUE, B = 10000))$p.value
  } else chi$p.value

  core <- tiers %>% filter(sharing_category == "Core / broad (>3 host species)")
  host <- tiers %>% filter(sharing_category == "Host-specific (1 host species)")
  ft <- fisher.test(matrix(c(core$n, core$n_without, host$n, host$n_without), nrow = 2, byrow = TRUE))

  tibble(
    OR = unname(ft$estimate), OR_ci_low = ft$conf.int[1], OR_ci_high = ft$conf.int[2],
    chi_p = chi_p, fisher_p = ft$p.value
  )
}

cat_stats <- cat_raw %>%
  group_by(cooperative_categories) %>%
  group_modify(~ test_category(.x)) %>%
  ungroup() %>%
  mutate(BH_p = p.adjust(fisher_p, method = "BH")) %>%
  mutate(
    direction = case_when(
      BH_p >= 0.05 ~ "Not significant",
      OR > 1 ~ "Core-enriched",
      OR < 1 ~ "Host-specific-enriched",
      TRUE ~ "Not significant"
    )
  )

write_csv(cat_stats, file.path(OUT_DIR, "table_cooperative_category_stats_R.csv"))
cat("\n=== Per-category stats (R-computed) ===\n")
print(as.data.frame(cat_stats))

# =============================================================================
# PANEL B — forest plot, styled after Hao et al. Fig. 2C/D
# =============================================================================
cat_stats_plot <- cat_stats %>%
  mutate(cooperative_categories = recode(cooperative_categories,
    "antibiotic_degradation" = "Antibiotic degradation",
    "siderophores" = "Siderophores",
    "biofilm" = "Biofilm formation",
    "secretion_system" = "Secretion systems",
    "quorum_sensing" = "Quorum sensing",
    "Extracellular" = "Extracellular proteins"
  )) %>%
  arrange(OR) %>%
  mutate(cooperative_categories = factor(cooperative_categories, levels = cooperative_categories))

panel_b <- ggplot(cat_stats_plot, aes(x = OR, y = cooperative_categories, color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = OR_ci_low, xmax = OR_ci_high), height = 0.15, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_color_manual(values = PAL, name = NULL) +
  labs(x = "Odds ratio (core/broad vs. host-specific, log scale)", y = NULL, title = "b") +
  theme_ms()

# =============================================================================
# PANEL C — prevalence % by category and tier
# =============================================================================
cat_long <- cat_raw %>%
  mutate(sharing_category = factor(sharing_category, levels = TIER_LEVELS, labels = TIER_LABELS)) %>%
  mutate(cooperative_categories = recode(cooperative_categories,
    "antibiotic_degradation" = "Antibiotic degradation",
    "siderophores" = "Siderophores",
    "biofilm" = "Biofilm formation",
    "secretion_system" = "Secretion systems",
    "quorum_sensing" = "Quorum sensing",
    "Extracellular" = "Extracellular proteins"
  ), .after = cooperative_categories) %>%
  mutate(cooperative_categories = factor(cooperative_categories, levels = levels(cat_stats_plot$cooperative_categories)))

panel_c <- ggplot(cat_long, aes(x = cooperative_categories, y = pct_genomes, fill = sharing_category)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = PAL, name = NULL) +
  labs(x = NULL, y = "% of genomes carrying \u22651 gene", title = "c") +
  coord_flip() +
  theme_ms()

# =============================================================================
# ASSEMBLE + EXPORT
# =============================================================================
fig7 <- (panel_a | panel_b) / panel_c + plot_layout(heights = c(1, 1.1))

ggsave(file.path(OUT_DIR, "Fig7_cooperative_traits.pdf"), fig7, width = 10, height = 9, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig7_cooperative_traits.png"), fig7, width = 10, height = 9, dpi = 600)

cat("\nFigure written to: ", OUT_DIR, "\n", sep = "")
cat("  Fig7_cooperative_traits.pdf / .png\n")
cat("  table_cooperative_category_stats_R.csv — stats behind panel b (use these OR/CI values in Table 11 if they differ from the earlier ad hoc check)\n")
