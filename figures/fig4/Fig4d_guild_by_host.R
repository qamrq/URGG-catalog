#!/usr/bin/env Rscript
# Fig. 4d: guild composition by host (hosts with >= MIN_HOST_N species; a species found in several hosts counts in each)
# Part of Figure 4 (thiamine guilds). Sources Fig4_common_setup.R from the same folder; paths via URGG_PROJECT.
.urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1]))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
  else getwd()
})
source(file.path(.urgg_script_dir, "Fig4_common_setup.R"))


ht <- if (file.exists(file.path(S66_DIR, "66_host_traits_all.csv"))) fread(file.path(S66_DIR, "66_host_traits_all.csv")) else NULL
hl <- melt(hg[host %in% keep], id.vars = c("host", "n_genomes"), variable.name = "guild", value.name = "pct")
hl[, guild := factor(as.character(guild), levels = rev(names(guild_palette)))]
ord <- hg[host %in% keep][order(`1_de_novo` + `2_partial_de_novo`), host]
hl[, host := factor(host, levels = ord)]
pBa <- ggplot(hl, aes(pct, host, fill = guild)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 8, fmtp(pct, 0), "")), position = position_stack(vjust = 0.5), size = 1.9,
            colour = ifelse(hl$guild %in% c("1_de_novo", "4_thiamine_salvage"), "white", NEUTRAL[["ink"]]), family = FONT) +
  scale_fill_manual(values = guild_palette, labels = guild_labels, breaks = names(guild_palette), name = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = "Genomes recovered from the host in each guild", y = NULL) + theme_catalogue() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank()) + guides(fill = guide_legend(nrow = 2))
pBa <- pBa + scale_y_discrete(labels = function(h) sprintf("%s (%s)", h, fmtn(hg[match(h, host), n_genomes]))) +
  theme(axis.text.y = element_text(face = "italic", size = 6.2))


save_all(pBa, "Fig4d_guild_by_host", 183, 100)
