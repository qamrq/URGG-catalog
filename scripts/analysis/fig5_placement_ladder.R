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
#  SCRIPT  : fig5_placement_ladder.R  -  Target placement summary, ladder layout (Fig. 5 supplement)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/placement_support_table.csv
#  OUTPUT  : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/fig_target_placement_summary.png
# ==========================================================

csv_path <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/placement_support_table.csv"
out_path <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/fig_target_placement_summary.png"

# Group colours follow the existing Fig. 5 palette in the reference layout
group_col <- c(
  "Outgroup"            = "#595959",
  "GTDB SK-Y3"          = "#E8A02A",
  "RUC + GTDB UMGS1840" = "#8E4EA8",
  "GTDB RGIG8785"       = "#4A9ACB"
)
target_col <- "#E31A1C"

# Internal node UFBoot values, read from the reference layout (not in the CSV).
# Tip order top to bottom for each run, and node values A, B, C, D along the
# ladder: A = outgroup node, B = SK-Y3 node, C = RUC node, D = RGIG8785 + sheep node.
ladder <- list(
  "Run A: Bacilli outgroup" = list(
    order = c("Outgroup", "GTDB SK-Y3", "RUC + GTDB UMGS1840", "GTDB RGIG8785"),
    nodes = c(A = 100, B = 100, C = 47)),
  "Run B: Oscillospirales outgroup" = list(
    order = c("Outgroup", "GTDB RGIG8785", "RUC + GTDB UMGS1840", "GTDB SK-Y3"),
    nodes = c(A = 100, B = 79, C = 58))
)

# Group sizes and support values per run
tab <- read.csv(csv_path, stringsAsFactors = FALSE)

# Ladder geometry in image pixel-like units (y grows downward)
node_x <- c(A = 70, B = 157, C = 243, D = 330)
node_y <- c(A = 345, B = 393, C = 441, D = 490)
tip_x  <- 414
tip_y  <- c(150, 248, 345, 442, 540)   # four groups + unclassified sheep genome
tri_w  <- 69
tri_h  <- 50
label_x <- 496
lwd_spine <- 2.4 * 1.5

draw_run <- function(run, spec) {
  groups <- split(tab[tab$run == run, ], tab$group[tab$run == run])

  # Spine A -> B -> C -> D
  segments(node_x[c("A", "B", "C")], node_y[c("A", "B", "C")],
           node_x[c("B", "C", "D")], node_y[c("B", "C", "D")],
           col = "black", lwd = lwd_spine)

  # Tip branches: A to outgroup, B to SK-Y3 slot, C to RUC slot, D to RGIG8785
  # and the sheep genome
  segments(node_x[c("A", "B", "C", "D")], node_y[c("A", "B", "C", "D")],
           tip_x, tip_y[1:4], col = "black", lwd = lwd_spine)
  segments(node_x["D"], node_y["D"], tip_x, tip_y[5], col = "black", lwd = lwd_spine)

  # Collapsed group triangles and labels
  for (idx in seq_along(spec$order)) {
    g <- spec$order[idx]
    y <- tip_y[idx]
    r <- groups[[g]]
    polygon(c(tip_x, tip_x + tri_w, tip_x + tri_w),
            c(y, y - tri_h / 2, y + tri_h / 2),
            col = group_col[[g]], border = group_col[[g]], lwd = 0.8)
    label <- if (is.na(r$group_clade_UFBoot) || r$group_clade_UFBoot == "NA") {
      sprintf("%s\nn = %d", g, r$n_tips)
    } else {
      sprintf("%s\nn = %d; UFBoot %d", g, r$n_tips,
              as.integer(as.numeric(r$group_clade_UFBoot)))
    }
    text(label_x, y, label, adj = c(0, 0.5), cex = 1.05)
  }

  # Unclassified sheep genome: star and bold red label
  points(tip_x + 4, tip_y[5], pch = 8, cex = 2.6, col = target_col, lwd = 2.2)
  text(tip_x + 30, tip_y[5], "Unclassified sheep genome", adj = c(0, 0.5),
       col = target_col, font = 2, cex = 1.1)

  # Node support values: black at 95 or higher, red below 95
  for (key in names(spec$nodes)) {
    val <- spec$nodes[[key]]
    col <- if (val >= 95) "black" else "#D0312D"
    text(node_x[key] - 10, node_y[key] + 20, as.character(val),
         col = col, font = 2, cex = 1.1, adj = c(1, 0.5))
  }

  title(main = run, adj = 0, cex.main = 1.3, font.main = 2)
}

png(out_path, width = 11, height = 5.4, units = "in", res = 300, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5), oma = c(3, 0, 0, 0))
for (run in names(ladder)) {
  plot.new()
  plot.window(xlim = c(40, 1000), ylim = c(600, 100))
  draw_run(run, ladder[[run]])
}
par(op)
mtext(paste(
  "Red star: unclassified sheep genome (target).",
  "Node values: UFBoot support; black = 95 or higher, red = below 95 (unresolved).",
  "Triangles: collapsed groups.", sep = "\n"),
  side = 1, line = 1, adj = 0, cex = 0.9, outer = TRUE)
invisible(dev.off())
