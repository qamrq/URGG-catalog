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
#  SCRIPT  : 10f_fig4e_hostpair_thiC_share.R  -  Fig. 4e, thiC share among shared vs all species per host pair
#  AUTHOR  : QRQ
#  INPUT   : see Fig4_common_setup.R
#  OUTPUT  : fig4/Fig4e_hostpair_thiC_share.{png,pdf}
# ==========================================================

.urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1]))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
  else getwd()
})
source(file.path(.urgg_script_dir, "_setup", "fig4_common_setup.R"))

# ---- c: shared vs all species, per pair ----
pc_d <- pp[n_shared >= 20 & host_a %in% keep & host_b %in% keep]
above <- pc_d[pct_denovo_among_shared > pct_denovo_among_union, .N]
wt <- if (nrow(pc_d) >= 5) suppressWarnings(wilcox.test(pc_d$pct_denovo_among_shared, pc_d$pct_denovo_among_union, paired = TRUE)) else NULL
abbr <- function(h) sub("^(\\w)\\w+ ", "\\1. ", h)
italic_legend <- function(...) {                                       # italic labels for species names, across ggplot2 versions
  if (packageVersion("ggplot2") >= "3.5.0") guide_legend(..., theme = theme(legend.text = element_text(face = "italic", size = 5.8)))
  else guide_legend(..., label.theme = element_text(face = "italic", size = 5.8))
}
# the five pairs sharing the most species get a colour (legend), the rest stay grey: no text on the plot
PAIR_COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
top_pairs <- pc_d[order(-n_shared)][1:min(5, .N)][, lab := sprintf("%s - %s (%s)", abbr(host_a), abbr(host_b), fmtn(n_shared))]
pc_d[, pair := "Other pairs"]
pc_d[top_pairs, on = .(host_a, host_b), pair := i.lab]
pair_levels <- c(top_pairs$lab, "Other pairs")
pc_d[, pair := factor(pair, levels = pair_levels)]
setorder(pc_d, -pair)                                                    # grey drawn first, coloured on top
rng <- range(c(pc_d$pct_denovo_among_union, pc_d$pct_denovo_among_shared), na.rm = TRUE)
lims <- c(floor(rng[1] / 5) * 5, ceiling(rng[2] / 5) * 5)
pBc <- ggplot(pc_d, aes(pct_denovo_among_union, pct_denovo_among_shared)) +
  geom_abline(slope = 1, intercept = 0, colour = NEUTRAL[["grey"]], linetype = "22", linewidth = 0.3) +
  geom_point(aes(size = n_shared, fill = pair), shape = 21, colour = "white", stroke = 0.35, alpha = 0.95) +
  scale_fill_manual(values = c(setNames(PAIR_COLS[seq_len(nrow(top_pairs))], top_pairs$lab), "Other pairs" = NEUTRAL[["light"]]),
                    breaks = pair_levels, name = "Host pair (shared species)") +
  scale_size_area(max_size = 5, name = "Shared species", breaks = c(50, 250, 1000), labels = fmtn) +
  guides(fill = italic_legend(order = 1, override.aes = list(size = 2.6)),
         size = guide_legend(order = 2, override.aes = list(fill = NEUTRAL[["light"]]))) +
  scale_x_continuous(limits = lims) + scale_y_continuous(limits = lims) +
  labs(x = "% thiC+ among all species of the two hosts", y = "% thiC+ among shared species",
       title = sprintf("%d of %d host pairs above the line%s", above, nrow(pc_d),
                       if (!is.null(wt)) sprintf(" (paired Wilcoxon p = %s)", format.pval(wt$p.value, digits = 2)) else "")) +
  theme_catalogue() + theme(plot.title = element_text(size = 6.4, face = "plain"), legend.position = "right", aspect.ratio = 1)
pBc <- pBc + theme(legend.text = element_text(size = 5.8), legend.spacing.y = unit(1.5, "mm"))

save_all(pBc, "Fig4e_hostpair_thiC_share", 90, 100)
