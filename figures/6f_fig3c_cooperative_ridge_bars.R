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
#  SCRIPT  : 6f_fig3c_cooperative_ridge_bars.R  -  Fig. 3c, cooperative genes per genome and by tier
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s80_figC_tables/figC_c_genome_props.csv
#            annotation/results_phase3/s80_figC_tables/figC_c_cooperative.csv
#  OUTPUT  : fig3/Fig3c_cooperative_ridge_bars.{png,pdf}
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "data")
URGG_FIGURES <- Sys.getenv("URGG_FIGURES", "figures_out")
urgg_in <- function(rel) {
  f <- file.path(URGG_PROJECT, rel)
  if (!file.exists(f)) stop("Input not found: ", f, " (set URGG_PROJECT)", call. = FALSE)
  f
}
urgg_out <- function(fig) { d <- file.path(URGG_FIGURES, fig); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
URGG_FONT <- local({
  f <- Sys.getenv("URGG_FONT", unset = "")
  if (!nzchar(f)) {
    f <- "sans"
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      sf <- systemfonts::system_fonts()
      for (fam in c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"))
        if (any(sf$family == fam & sf$italic)) { f <- fam; break }
    }
  }
  f
})
urgg_save <- function(p, dir, name, w_mm, h_mm) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 600, bg = "white",
                  device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm",
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  if (requireNamespace("svglite", quietly = TRUE))
    ggplot2::ggsave(file.path(dir, paste0(name, ".svg")), p, width = w_mm, height = h_mm, units = "mm", device = svglite::svglite)
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

IN_PROPS <- urgg_in("annotation/results_phase3/s80_figC_tables/figC_c_genome_props.csv")
IN_CATS  <- urgg_in("annotation/results_phase3/s80_figC_tables/figC_c_cooperative.csv")
OUT      <- urgg_out("fig3")

FONT <- URGG_FONT; FS <- 12
TIERS <- c("Core", "Shared", "Host-specific")
COL  <- c(Core = "#B06040", Shared = "#C89040", "Host-specific" = "#2898B8")
GREY_ALPHA <- 0.35
theme_fig <- function(base_size = FS) theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), plot.margin = margin(4, 8, 4, 4))

# ---- left: ridgeline ---------------------------------------------------------
genes <- read.csv(IN_PROPS, stringsAsFactors = FALSE)
genes <- genes[genes$n_ko > 0, ]
genes$pct <- 100 * genes$n_coop / genes$n_ko
stopifnot(setequal(unique(genes$tier), TIERS))
lo <- as.numeric(quantile(genes$pct, 0.001)); hi <- 6.0
dens <- lapply(TIERS, function(t) { v <- genes$pct[genes$tier == t]; d <- density(v, from = lo, to = hi, n = 400)
  data.frame(x = d$x, d = d$y, tier = t, n = length(v), med = median(v)) })
names(dens) <- TIERS
scale_max <- max(sapply(dens, function(z) max(z$d)))
spacing <- 1.25; base <- c(Core = 2 * spacing, Shared = spacing, "Host-specific" = 0)
ridge_df <- do.call(rbind, lapply(TIERS, function(t) data.frame(x = dens[[t]]$x, ymin = base[[t]], ymax = base[[t]] + dens[[t]]$d / scale_max, tier = t)))
ridge_df$tier <- factor(ridge_df$tier, levels = rev(TIERS))      # core drawn last (on top)
stats_df <- do.call(rbind, lapply(TIERS, function(t) data.frame(tier = t, n = dens[[t]]$n[1], med = dens[[t]]$med[1], base = base[[t]])))
message(sprintf("Median %% of KO-annotated genes: %s", paste(sprintf("%s %.2f", stats_df$tier, stats_df$med), collapse = ", ")))

p_ridge <- ggplot() +
  geom_ribbon(data = ridge_df, aes(x = x, ymin = ymin, ymax = ymax, fill = tier), colour = "white", linewidth = 0.3) +
  geom_segment(data = stats_df, aes(x = med, xend = med, y = base, yend = base + 0.92), colour = "#1A1A1A", linewidth = 0.6) +
  geom_text(data = stats_df, aes(x = med, y = base + 1.02, label = paste0("(n = ", format(n, big.mark = ",", trim = TRUE), ")"), colour = tier),
            family = FONT, fontface = "bold", size = FS / .pt, vjust = 0) +
  scale_fill_manual(values = COL, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = c(2, 4, 6), limits = c(0, 8.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 2 * spacing + 1.45), expand = c(0, 0)) +
  labs(x = "% of KO-annotated genes", y = NULL) + theme_fig() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

# ---- right: carriage bars + odds ratios (read from the table, not typed in) ---
ORDER <- c("Antibiotic degradation", "Siderophores", "Extracellular proteins", "Quorum sensing", "Secretion systems", "Biofilm formation")
cats <- read.csv(IN_CATS, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(setequal(cats$label, ORDER))
cats <- cats[match(ORDER, cats$label), ]
# grey = unadjusted difference not significant (adjusted p >= 0.05 or CI including 1), as in the paper's rule
cats$grey <- !(cats$p_u < 0.05 & (cats$lo_u > 1 | cats$hi_u < 1))
BARS <- data.frame(trait = factor(cats$label, levels = rev(ORDER)), Core = cats$core, Shared = cats$shared,
                   "Host-specific" = cats$host, or_u = cats$or_u, or_w = cats$or_w, grey = cats$grey, check.names = FALSE)
row_y <- setNames(rev(seq_len(nrow(BARS))) - 1, as.character(BARS$trait))
OFF <- c(Core = 0.28, Shared = 0.0, "Host-specific" = -0.28); H <- 0.27
bar_rows <- do.call(rbind, lapply(TIERS, function(t) data.frame(trait = BARS$trait, tier = t, value = BARS[[t]],
  y = as.numeric(row_y[as.character(BARS$trait)]) + OFF[[t]], grey = BARS$grey)))
bar_rows$tier <- factor(bar_rows$tier, levels = TIERS)
bar_rows$alpha <- ifelse(bar_rows$grey, GREY_ALPHA, 1)
bar_rows$label_inside <- bar_rows$value >= 10
bar_rows$lab_x <- ifelse(bar_rows$label_inside, bar_rows$value - 1, bar_rows$value + 1)
bar_rows$lab_hjust <- ifelse(bar_rows$label_inside, 1, 0)
bar_rows$lab_col <- ifelse(bar_rows$label_inside & bar_rows$alpha >= 1, "white", ifelse(bar_rows$label_inside, "#2B2B2B", "#7F7F7F"))
or_df <- data.frame(y = as.numeric(row_y[as.character(BARS$trait)]), lab = sprintf("OR %.2f | %.2f", BARS$or_u, BARS$or_w),
                    alpha = ifelse(BARS$grey, GREY_ALPHA, 1))

p_bars <- ggplot() +
  geom_rect(data = bar_rows, aes(xmin = 0, xmax = value, ymin = y - H / 2, ymax = y + H / 2, fill = tier, alpha = alpha), colour = NA) +
  geom_text(data = bar_rows, aes(x = lab_x, y = y, label = sprintf("%.1f", value), hjust = lab_hjust, colour = lab_col), family = FONT, size = (FS - 1) / .pt) +
  geom_text(data = or_df, aes(x = 106, y = y, label = lab, alpha = alpha), hjust = 0, family = FONT, fontface = "bold", size = FS / .pt, colour = "#2B2B2B") +
  scale_fill_manual(values = COL, name = NULL) + scale_alpha_identity() + scale_colour_identity() +
  scale_x_continuous(limits = c(0, 140), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = as.numeric(row_y), labels = names(row_y), limits = c(-0.6, nrow(BARS) - 0.4), expand = c(0, 0)) +
  labs(x = "Genomes carrying category (%)", y = NULL) + theme_fig() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = FS + 3, family = FONT), legend.position = "none")

fig_c <- p_ridge + p_bars + plot_layout(widths = c(0.22, 0.78)) +
  plot_annotation(caption = "OR: unadjusted | within lineages", theme = theme(plot.caption = element_text(hjust = 1, family = FONT, face = "bold")))
urgg_save(fig_c, OUT, "Fig3c_cooperative_ridge_bars", 360, 150)
