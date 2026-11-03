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
#  SCRIPT  : 14f_fig6e_viral_class_and_rank.R  -  Fig. 6e, viral class by tier and deepest assigned rank
#  AUTHOR  : QRQ
#  INPUT   : catalog_tables/web_votus.tsv
#  OUTPUT  : fig6/Fig6e_viral_class_and_rank.{png,pdf}
# ==========================================================

.urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1]))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
  else getwd()
})
source(file.path(.urgg_script_dir, "_setup", "urgg_paths.R"))

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

IN  <- urgg_in("catalog_tables/web_votus.tsv")
OUT <- urgg_out("fig6")
FONT <- URGG_FONT; FS <- 11
TIERS <- c("Host-specific", "Shared", "Core")
TIER_COL <- c("Host-specific" = "#2898B8", "Shared" = "#C89040", "Core" = "#B06040")
TEXT_COL <- c("Host-specific" = "white", "Shared" = "#1A1A1A", "Core" = "white")
MIN_INSIDE <- 6      # segments below this % of the class are labelled above the bar

v <- read.delim(IN, stringsAsFactors = FALSE, na.strings = character(0), check.names = FALSE)
vv <- v[v$phamer_call == "Viral", ]
vv$class[is.na(vv$class) | vv$class == ""] <- "No class assigned"
stopifnot(nrow(vv) == 6619, all(vv$sharing_tier %in% TIERS))
message(sprintf("Viral vOTUs: %d | %s", nrow(vv), paste(sprintf("%s %d", TIERS, as.integer(table(factor(vv$sharing_tier, TIERS)))), collapse = ", ")))

# ---- left: class by tier ------------------------------------------------------
tab <- as.data.frame.matrix(table(vv$class, factor(vv$sharing_tier, levels = TIERS)))
tab$total <- rowSums(tab); tab$class <- rownames(tab)
tab <- tab[order(-tab$total), ]
stopifnot(tab["Caudoviricetes", "total"] == 6284)
x_pos <- setNames(seq_len(nrow(tab)), tab$class)
long <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
  n <- unlist(tab[i, TIERS]); pct <- 100 * n / tab$total[i]
  data.frame(class = tab$class[i], x = i, tier = TIERS, n = n, pct = pct,
             ymin = c(0, head(cumsum(pct), -1)), ymax = cumsum(pct), total = tab$total[i])
}))
long$ymid <- (long$ymin + long$ymax) / 2
long$tier <- factor(long$tier, levels = TIERS)
inside <- long[long$pct >= MIN_INSIDE & long$n > 0, ]
inside$label <- sprintf("%s\n(%.1f%%)", format(inside$n, big.mark = ",", trim = TRUE), inside$pct)
small <- long[long$pct < MIN_INSIDE & long$n > 0, ]
callout <- do.call(rbind, lapply(split(small, small$x), function(z)
  data.frame(x = z$x[1], label = paste(sprintf("%s %s (%.1f%%)", z$tier, format(z$n, big.mark = ",", trim = TRUE), z$pct), collapse = "\n"))))
if (!is.null(callout)) callout$y <- 101.5 + 8 * ((seq_len(nrow(callout)) - 1) %% 2)   # alternate heights so neighbours do not overlap

p_left <- ggplot() +
  geom_rect(data = long, aes(xmin = x - 0.36, xmax = x + 0.36, ymin = ymin, ymax = ymax, fill = tier), colour = "white", linewidth = 0.3) +
  geom_text(data = inside, aes(x = x, y = ymid, label = label, colour = TEXT_COL[as.character(tier)]), family = FONT, size = (FS - 2) / .pt, lineheight = 0.9) +
  { if (!is.null(callout)) geom_text(data = callout, aes(x = x, y = y, label = label), family = FONT, size = (FS - 3) / .pt, lineheight = 0.9, vjust = 0, colour = "#1A1A1A") } +
  scale_fill_manual(values = TIER_COL, name = NULL) + scale_colour_identity() +
  scale_x_continuous(breaks = x_pos, labels = paste0(names(x_pos), "\nn = ", format(tab$total, big.mark = ",", trim = TRUE)), expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 124), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "vOTUs in the class (%)") +
  theme_classic(base_size = FS, base_family = FONT) +
  theme(axis.text.x = suppressWarnings(element_text(angle = 35, hjust = 1, vjust = 1, face = ifelse(names(x_pos) == "No class assigned", "plain", "italic"))), axis.line = element_line(linewidth = 0.3),
        legend.position = "top", plot.margin = margin(6, 10, 6, 6))

# ---- right: deepest assigned rank (root-only vOTUs omitted, as in the legend) --
RANKS <- c("Class", "Order", "Family", "Subfamily", "Genus", "Species")
rk <- as.data.frame(table(factor(vv$deepest_rank, levels = c(RANKS, "None (root only)"))))
names(rk) <- c("rank", "n"); rk$pct <- 100 * rk$n / nrow(vv)
message(sprintf("Root-only vOTUs omitted from the rank panel: %d", rk$n[rk$rank == "None (root only)"]))
rk <- rk[rk$rank %in% RANKS, ]; rk$rank <- factor(rk$rank, levels = rev(RANKS))
p_right <- ggplot(rk, aes(pct, rank)) +
  geom_col(width = 0.65, fill = "#5A6B7B") +
  geom_text(aes(label = sprintf("%s (%.1f%%)", format(n, big.mark = ",", trim = TRUE), pct)), hjust = -0.08, family = FONT, size = (FS - 2) / .pt) +
  scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "vOTUs called viral (%)", y = "Deepest assigned rank") +
  theme_classic(base_size = FS, base_family = FONT) + theme(axis.line = element_line(linewidth = 0.3), plot.margin = margin(6, 14, 6, 6))

fig <- p_left + p_right + plot_layout(widths = c(0.68, 0.32))
urgg_save(fig, OUT, "Fig6e_viral_class_and_rank", 300, 150)
