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
#  SCRIPT  : 31k_fig_quality_scatter.R  -  Completeness vs contamination scatter (supplement)
#  AUTHOR  : QRQ
#  INPUT   : final_reps: genome, completeness, contamination, quality_tier
#  OUTPUT  : Scatter figure (PDF, PNG) in OUT_DIR
# ==========================================================

###############################################################################
# (3b) MAG QUALITY SCATTER — completeness vs contamination, final representatives
#
# Drop-in for the existing pipeline. Expects, from section (3):
#   final_reps  : data.frame with genome, completeness, contamination, quality_tier
#   has_contam  : logical
#   OUT_DIR     : output directory
#   theme_pub   : the project ggplot theme
#
# Thresholds deliberately follow THIS PROJECT's QC (HQ >90% & <=5% contam;
# MQ 70-90% & <=5% contam), not generic MIMAG. The MIMAG 10% contamination
# band and the QS50 diagonal are omitted on purpose: every genome in this
# pool already cleared <=5% upstream, so both would sit outside the data.
###############################################################################

if (!has_contam) {

  message("No contamination column in the provenance chain — ",
          "skipping the quality scatter (section 3b).")

} else {

  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(patchwork)
  })

  ## --- guard: one row per representative -----------------------------------
  if (anyDuplicated(final_reps$genome)) {
    warning(sprintf(
      "%d representatives have >1 completeness/contamination record; taking the first per genome.",
      sum(duplicated(final_reps$genome))))
    final_reps <- final_reps %>% distinct(genome, .keep_all = TRUE)
  }

  qc_dat <- final_reps %>% filter(!is.na(completeness), !is.na(contamination))
  n_drop <- nrow(final_reps) - nrow(qc_dat)
  if (n_drop > 0)
    message(sprintf("Dropping %d representative(s) with missing completeness/contamination.", n_drop))

  ## --- palette: NAMED so tiers can't be mis-mapped -------------------------
  tier_cols <- c(
    "High-quality (HQ)"   = "#1B5E5A",
    "Medium-quality (MQ)" = "#C98B3A"
  )
  other_lab <- grep("^Other", unique(qc_dat$quality_tier), value = TRUE)
  if (length(other_lab)) tier_cols[other_lab] <- "#B0B0B0"

  tier_n <- qc_dat %>%
    count(quality_tier, name = "n") %>%
    mutate(lab = sprintf("%s  (%s; %.1f%%)", quality_tier,
                         format(n, big.mark = ","), 100 * n / sum(n)))
  lab_map <- setNames(tier_n$lab, tier_n$quality_tier)
  qc_dat  <- qc_dat %>% mutate(tier_lab = lab_map[quality_tier])
  leg_cols <- setNames(unname(tier_cols[names(lab_map)]), unname(lab_map))

  ## --- axis ranges follow the data, not fixed MIMAG limits -----------------
  x_lo  <- min(65, floor(min(qc_dat$completeness) / 5) * 5)
  y_hi  <- max(5.4, ceiling(max(qc_dat$contamination) * 1.08 * 2) / 2)

  ## --- point layer: rasterised if ggrastr is available ---------------------
  pt <- geom_point(alpha = 0.35, size = 0.45, stroke = 0, shape = 16)
  if (requireNamespace("ggrastr", quietly = TRUE))
    pt <- ggrastr::rasterise(pt, dpi = 450)   # keeps the PDF small & openable

  ## --- main panel ----------------------------------------------------------
  p_main <- ggplot(qc_dat, aes(completeness, contamination, colour = tier_lab)) +
    annotate("rect", xmin = 90, xmax = Inf, ymin = -Inf, ymax = 5,
             fill = "#1B5E5A", alpha = 0.055) +
    geom_vline(xintercept = c(70, 90), linetype = "dotted",
               colour = "grey20", linewidth = 0.35) +
    geom_hline(yintercept = 5, linetype = "dotted",
               colour = "grey20", linewidth = 0.35) +
    pt +
    annotate("text", x = 90.4, y = y_hi, label = "90%", hjust = 0, vjust = 1,
             size = 2.9, colour = "grey40") +
    annotate("text", x = 70.4, y = y_hi, label = "70%", hjust = 0, vjust = 1,
             size = 2.9, colour = "grey40") +
    annotate("text", x = x_lo + 0.4, y = 5.08, label = "5% contamination",
             hjust = 0, vjust = 0, size = 2.9, colour = "grey40") +
    scale_colour_manual(values = leg_cols, name = NULL) +
    coord_cartesian(xlim = c(x_lo, 100.4), ylim = c(-0.12, y_hi), expand = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    labs(x = "Completeness (%)", y = "Contamination (%)") +
    theme_pub +
    theme(legend.position = c(0.985, 0.985),
          legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 8.4))

  ## --- marginal histograms -------------------------------------------------
  p_top <- ggplot(qc_dat, aes(completeness)) +
    geom_histogram(bins = 50, fill = "#1B5E5A", alpha = 0.55, colour = NA) +
    geom_vline(xintercept = 90, linetype = "dotted", colour = "grey20",
               linewidth = 0.35) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = c(x_lo, 100.4), expand = FALSE) +
    theme_void()

  p_right <- ggplot(qc_dat, aes(y = contamination)) +
    geom_histogram(bins = 50, fill = "#C98B3A", alpha = 0.6, colour = NA) +
    geom_hline(yintercept = 5, linetype = "dotted", colour = "grey20",
               linewidth = 0.35) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(-0.12, y_hi), expand = FALSE) +
    theme_void()

  subtitle <- sprintf(
    "n = %s dereplicated representatives   \u00b7   median completeness %.1f%%   \u00b7   median contamination %.2f%%",
    format(nrow(qc_dat), big.mark = ","),
    median(qc_dat$completeness), median(qc_dat$contamination))

  pScatter <- (p_top + plot_spacer() + p_main + p_right) +
    plot_layout(ncol = 2, widths = c(4.4, 1), heights = c(1, 4.4)) +
    plot_annotation(
      title = "Quality of the final representative genomes",
      subtitle = subtitle,
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.6, colour = "grey35",
                                     margin = margin(b = 6))))

  ggsave(file.path(OUT_DIR, "Fig_S_completeness_contamination_scatter.pdf"),
         pScatter, width = 7.2, height = 6.6)
  ggsave(file.path(OUT_DIR, "Fig_S_completeness_contamination_scatter.png"),
         pScatter, width = 7.2, height = 6.6, dpi = 400)

  cat("  Fig_S_completeness_contamination_scatter\n")
}
