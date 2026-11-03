#!/usr/bin/env Rscript
# Fig. 4f: host network of shared thiC-carrying species (edge width, shared species; colour, shared/all ratio)
# Part of Figure 4 (thiamine guilds). Sources Fig4_common_setup.R from the same folder; paths via URGG_PROJECT.
.urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1]))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
  else getwd()
})
source(file.path(.urgg_script_dir, "Fig4_common_setup.R"))


# ---- b: network of shared producers ----
fam_order <- c("Bovidae", "Cervidae", "Giraffidae", "Camelidae")
nodes_b <- data.table(host = keep)[, family := factor(host_family[host], levels = fam_order)]
nodes_b <- merge(nodes_b, hg[, .(host, n_genomes)], by = "host")
setorder(nodes_b, family, -n_genomes)
k <- nrow(nodes_b); ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k
nodes_b[, `:=`(x = cos(ang), y = sin(ang))]
ed <- pp[host_a %in% keep & host_b %in% keep & n_shared_denovo >= MIN_EDGE]
ed <- merge(ed, nodes_b[, .(host_a = host, xa = x, ya = y)], by = "host_a")
ed <- merge(ed, nodes_b[, .(host_b = host, xb = x, yb = y)], by = "host_b")
setorder(ed, n_shared_denovo)
pBb <- ggplot() +
  geom_segment(data = ed, aes(x = xa, y = ya, xend = xb, yend = yb, linewidth = n_shared_denovo, colour = denovo_shared_vs_union),
               lineend = "round", alpha = 0.9) +
  geom_point(data = nodes_b, aes(x, y, size = n_genomes), shape = 21, fill = species_palette[nodes_b$host], colour = "white", stroke = 0.7) +
  geom_text(data = nodes_b, aes(1.17 * x, 1.17 * y, label = host, hjust = ifelse(x > 0.15, 0, ifelse(x < -0.15, 1, 0.5)),
                                vjust = ifelse(abs(x) <= 0.15, ifelse(y > 0, 0, 1), 0.5)),
            size = 2.1, fontface = "italic", family = FONT, colour = NEUTRAL[["ink"]]) +
  scale_linewidth_continuous(range = c(0.25, 3.4), name = "Shared thiC+ species") +
  scale_size(range = c(1.6, 8), trans = "sqrt", name = "Genomes", breaks = c(150, 1000, 5000), labels = fmtn) +
  scale_colour_gradientn(colours = c(NEUTRAL[["light"]], SALVAGER[["light"]], SALVAGER[["mid"]], SALVAGER[["dark"]]),
                         values = c(0, 0.25, 0.6, 1), limits = c(0.75, 2.75), oob = scales::squish,
                         breaks = c(1, 1.5, 2, 2.5), name = "Producer share,\nshared / all species") +
  guides(size = guide_legend(override.aes = list(fill = NEUTRAL[["light"]], colour = NEUTRAL[["light"]]), order = 1),
         linewidth = guide_legend(order = 2, override.aes = list(colour = SALVAGER[["mid"]])), colour = guide_colourbar(order = 3, barwidth = unit(2.5, "mm"), barheight = unit(18, "mm"))) +
  coord_fixed(xlim = c(-1.6, 1.6), ylim = c(-1.22, 1.22), clip = "off") + theme_void(base_family = FONT) +
  theme(legend.position = "right", legend.box = "vertical", legend.spacing.y = unit(2, "mm"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6), legend.title = element_text(size = 6, face = "bold"),
        plot.margin = margin(4, 10, 4, 10))


save_all(pBb, "Fig4f_host_network", 183, 100)
