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
#  SCRIPT  : fig5_placement_tree.R  -  Target placement summary, tree layout (Fig. 5 supplement)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/placement_support_table.csv
#  OUTPUT  : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/fig_target_placement_tree.png
# ==========================================================

csv_path <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/placement_support_table.csv"
out_path <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/fig_target_placement_tree.png"

# Top to bottom in the panel
group_order <- c("RUC + GTDB UMGS1840", "GTDB SK-Y3", "GTDB RGIG8785", "Outgroup")
group_col <- c(
  "RUC + GTDB UMGS1840" = "#1F78B4",
  "GTDB SK-Y3"          = "#E31A1C",
  "GTDB RGIG8785"       = "#33A02C",
  "Outgroup"            = "#999999"
)

tab <- read.csv(csv_path, stringsAsFactors = FALSE, colClasses = c(
  run = "character", group = "character", n_tips = "integer",
  group_clade_UFBoot = "character"))

root_x <- 0.0
node_x <- 1.0
tip_x  <- 2.2
ylim   <- c(-0.6, length(group_order) - 0.4)

draw_tree <- function(run) {
  rows <- tab[tab$run == run, ]
  present <- group_order[group_order %in% rows$group]
  ys <- setNames(length(group_order) - 1 - (match(present, group_order) - 1), present)
  ingroup <- present[present != "Outgroup"]
  y_mid <- mean(ys[ingroup])

  # Root and ingroup/outgroup split
  segments(root_x, y_mid, node_x, y_mid, lwd = 1.4)
  segments(node_x, min(ys[ingroup]), node_x, max(ys[ingroup]), lwd = 1.4)

  for (g in present) {
    r <- rows[rows$group == g, ]
    y <- ys[[g]]
    if (g == "Outgroup") {
      segments(root_x, y_mid, root_x, y, lwd = 1.4)
      segments(root_x, y, tip_x, y, lwd = 1.4)
    } else {
      segments(node_x, y, tip_x, y, lwd = 1.4)
    }

    # Collapsed-clade triangle, width scaled to tip count
    height <- 0.32
    width  <- 0.35 + 0.22 * sqrt(r$n_tips) / sqrt(101) * 3
    polygon(c(tip_x, tip_x, tip_x + width),
            c(y - height / 2, y + height / 2, y),
            col = adjustcolor(group_col[[g]], alpha.f = 0.9),
            border = "black", lwd = 0.6)

    # Support value under the group label, clear of the branch lines
    if (is.na(r$group_clade_UFBoot)) {
      sup_text <- "UFBoot NA"
      sup_col  <- "#666666"
      sup_font <- 1
    } else {
      v <- as.integer(as.numeric(r$group_clade_UFBoot))
      sup_text <- sprintf("UFBoot %d", v)
      sup_col  <- if (v >= 95) "black" else "#D95F02"
      sup_font <- 2
    }
    text(tip_x + width + 0.06, y, sprintf("%s\n(n = %d tips)", g, r$n_tips),
         adj = c(0, 0.5), cex = 0.9)
    text(tip_x + width + 0.06, y - 0.27, sup_text, adj = c(0, 0.5),
         cex = 0.8, col = sup_col, font = sup_font)
  }
  title(main = run, adj = 0, cex.main = 1.1)
}

png(out_path, width = 11, height = 5.2, units = "in", res = 300, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5), oma = c(2.5, 0, 2, 0))
for (run in unique(tab$run)) {
  plot.new()
  plot.window(xlim = c(-0.2, 5.0), ylim = ylim)
  draw_tree(run)
}
par(op)
mtext("Placement of the target genome and reference groups across two outgroup choices",
      side = 3, line = 0.5, outer = TRUE, adj = 0, cex = 1.2, font = 2)
mtext(paste(
  "Collapsed clades; triangle width scales with tip count. UFBoot support under each group;",
  "values below 95 in orange. Relationships among groups are not shown in the source table."),
  side = 1, line = 0.5, outer = TRUE, adj = 0, cex = 0.7, col = "#444444")
invisible(dev.off())
