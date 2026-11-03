#!/usr/bin/env Rscript
###############################################################################
# 66_thiamine_by_host.R
#
# Which ruminant hosts' gut genomes carry (and share) thiamine capacity?
# Uses the per-genome guilds from 63_ and the host list of each representative
# (representative_host_range.csv: host_species_list).
#
#   (a) Per host: % of genomes recovered from that host in each guild, and %
#       carrying thiC / complete-or-partial de novo / TPK / transporter,
#       with Wilson 95% CIs. Two versions: all genomes recovered from the host
#       (shared genomes count once in every host they occur in), and genomes
#       found ONLY in that host.
#   (b) Per host vs all other genomes: Fisher OR for de novo capacity (thiC),
#       BH across hosts.
#   (c) Host x host: species clusters shared by each pair of hosts, and how many
#       of those shared clusters can make thiamine de novo (thiC+), compared with
#       the de novo share among all clusters of the two hosts.
#   (d) The thiC+ clusters shared by the most hosts, with genus.
#
# Caveats printed in the report: host-level % reflects which lineages were
# recovered from each host (sampling depth, study number, diet), not a test of
# host effect; hosts with < MIN_HOST_N genomes are flagged.
#
# OUTPUTS (OUT_DIR)
#   66_host_guilds_all.csv, 66_host_guilds_hostonly.csv, 66_host_traits_all.csv,
#   66_host_traits_hostonly.csv, 66_host_thiC_fisher.csv,
#   66_hostpair_shared_producers.csv, 66_shared_producer_clusters.csv,
#   66_hostpair_heatmap.png (if ggplot2 is installed)
###############################################################################

suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
A  <- "/disk252/qamar/metagenome_analysis/annotation"
R3 <- file.path(A, "results_phase3")
GUILD_F    <- P("GUILD_F", file.path(R3, "s63_thiamine_guild_propagated/63_genome_thiamine_guilds.tsv"))
SHARING_F  <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
GTDB_BAC_F <- P("GTDB_BAC_F", file.path(A, "gtdbtk_output_v232/gtdbtk.bac120.summary.tsv"))
GTDB_AR_F  <- P("GTDB_AR_F", file.path(A, "gtdbtk_output_v232/gtdbtk.ar53.summary.tsv"))
OUT_DIR    <- P("OUT_DIR", file.path(R3, "s66_thiamine_by_host"))
MIN_HOST_N <- as.integer(P("MIN_HOST_N", "100"))
options(width = 220)
for (f in c(GUILD_F, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "66_run_report.txt"), split = TRUE)
hr  <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",")
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n; c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = round(100 * pmax(0, c - h), 1), hi = round(100 * pmin(1, c + h), 1))
}

# =============================================================================
# 0. Inputs
# =============================================================================
sh <- fread(SHARING_F)
setnames(sh, "representative", "genome")
if (!"host_species_list" %in% names(sh)) stop("representative_host_range.csv has no host_species_list column; columns are: ",
                                              paste(names(sh), collapse = ", "))
g <- fread(GUILD_F, select = c("genome", "guild", "hmp_denovo", "tpk_thiN", "transport", "tmp_kinase_thiL", "denovo_module_pct"))
g[, denovo_capable := hmp_denovo]                      # thiC: committed step; guilds 1 + 2
dat <- merge(sh[, .(genome, sharing_category, host_species_list)], g, by = "genome")
cat(sprintf("Genomes with guild calls and host lists: %s\n", fmt(nrow(dat))))
hs <- dat[, .(host = trimws(unlist(strsplit(host_species_list, "[;,|]")))), by = genome][host != ""]
hs <- unique(hs)
hs <- merge(hs, dat[, .(genome, sharing_category, guild, denovo_capable, tpk_thiN, transport, tmp_kinase_thiL)], by = "genome")
hs[, n_hosts := uniqueN(host), by = genome]
cat(sprintf("Genome-host memberships: %s across %d host species\n", fmt(nrow(hs)), uniqueN(hs$host)))
chk <- hs[, .(n_hosts = n_hosts[1]), by = .(genome, sharing_category)][, .(min = min(n_hosts), max = max(n_hosts)), by = sharing_category]
cat("Host counts per genome by tier (sanity check: host-specific = 1, shared = 2-3, core/broad > 3):\n"); print(chk)

tax <- rbindlist(lapply(c(GTDB_BAC_F, GTDB_AR_F), function(f)
  if (file.exists(f)) fread(f, sep = "\t", quote = "", select = c("user_genome", "classification")) else NULL))
if (nrow(tax)) {
  setnames(tax, c("genome", "classification"))
  tax[, genus := fifelse(grepl("g__", classification), sub("^.*g__([^;]*).*$", "\\1", classification), NA_character_)]
  tax[is.na(genus) | genus == "", genus := "(no genus)"]
}

# =============================================================================
# (a) Per-host composition
# =============================================================================
host_tables <- function(d, label) {
  n <- d[, .(n_genomes = uniqueN(genome)), by = host]
  gd <- dcast(d[, .N, by = .(host, guild)], host ~ guild, value.var = "N", fill = 0)
  gd <- merge(n, gd, by = "host")
  for (c in setdiff(names(gd), c("host", "n_genomes"))) gd[, (c) := round(100 * get(c) / n_genomes, 1)]
  tr <- d[, .(n_genomes = uniqueN(genome), n_denovo = sum(denovo_capable), n_TPK = sum(tpk_thiN),
              n_transport = sum(transport), n_thiL = sum(tmp_kinase_thiL)), by = host]
  for (v in c("denovo", "TPK", "transport", "thiL")) {
    ci <- wilson(tr[[paste0("n_", v)]], tr$n_genomes)
    tr[, (paste0("pct_", v)) := round(100 * get(paste0("n_", v)) / n_genomes, 1)]
    tr[, (paste0("pct_", v, "_ci")) := sprintf("%.1f-%.1f", ci$lo, ci$hi)]
  }
  tr[, flag := fifelse(n_genomes < MIN_HOST_N, sprintf("< %d genomes", MIN_HOST_N), "")]
  setorder(gd, -n_genomes); setorder(tr, -pct_denovo)
  cat(sprintf("\n[%s] guild %% by host:\n", label)); print(gd)
  cat(sprintf("\n[%s] capacity %% by host (sorted by de novo, thiC+):\n", label))
  print(tr[, .(host, n_genomes, pct_denovo, pct_denovo_ci, pct_TPK, pct_transport, pct_thiL, flag)])
  list(guilds = gd, traits = tr)
}
hr("(a) Thiamine guilds and capacity by host species")
cat("All genomes recovered from each host (a shared genome counts in every host it occurs in):\n")
all_h <- host_tables(hs, "all genomes")
fwrite(all_h$guilds, file.path(OUT_DIR, "66_host_guilds_all.csv")); fwrite(all_h$traits, file.path(OUT_DIR, "66_host_traits_all.csv"))
cat("\n\nGenomes found ONLY in that host (host-specific tier):\n")
only_h <- host_tables(hs[n_hosts == 1], "host-only genomes")
fwrite(only_h$guilds, file.path(OUT_DIR, "66_host_guilds_hostonly.csv")); fwrite(only_h$traits, file.path(OUT_DIR, "66_host_traits_hostonly.csv"))

# =============================================================================
# (b) Per host vs rest: de novo capacity
# =============================================================================
hr("(b) De novo capacity (thiC) per host vs all other genomes")
gl <- unique(hs[, .(genome, denovo_capable)])
fb <- rbindlist(lapply(sort(unique(hs$host)), function(h) {
  ing <- unique(hs[host == h, genome]); inh <- gl$genome %chin% ing
  a <- sum(gl$denovo_capable & inh); b <- sum(!gl$denovo_capable & inh)
  c <- sum(gl$denovo_capable & !inh); d <- sum(!gl$denovo_capable & !inh)
  ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
  data.table(host = h, n_genomes = a + b, pct_host = round(100 * a / (a + b), 1), pct_rest = round(100 * c / (c + d), 1),
             OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = ft$p.value)
}))
fb[, p_adj_BH := signif(p.adjust(p, "BH"), 3)]
fb[, call := fcase(p_adj_BH < 0.05 & OR_lo > 1, "higher", p_adj_BH < 0.05 & OR_hi < 1, "lower", default = "ns")]
print(fb[order(-OR), -"p"]); fwrite(fb, file.path(OUT_DIR, "66_host_thiC_fisher.csv"))
cat("Not taxonomy-adjusted: a host's value reflects which lineages were recovered from it.\n")

# =============================================================================
# (c) Host x host sharing of thiamine producers
# =============================================================================
hr("(c) Species clusters shared between host pairs, and how many make thiamine (thiC+)")
hosts <- sort(unique(hs$host))
pairs <- CJ(host_a = hosts, host_b = hosts)[host_a < host_b]
per_host <- split(hs[, .(genome, denovo_capable)], hs$host)
pp <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  a <- per_host[[pairs$host_a[i]]]; b <- per_host[[pairs$host_b[i]]]
  shared <- intersect(a$genome, b$genome); union_g <- union(a$genome, b$genome)
  dn <- unique(rbind(a, b))
  data.table(host_a = pairs$host_a[i], host_b = pairs$host_b[i], n_union = length(union_g), n_shared = length(shared),
             n_shared_denovo = dn[genome %chin% shared, sum(denovo_capable)],
             pct_denovo_among_shared = if (length(shared)) round(100 * dn[genome %chin% shared, mean(denovo_capable)], 1) else NA_real_,
             pct_denovo_among_union = round(100 * dn[, mean(denovo_capable)], 1))
}))
pp[, jaccard_all := round(n_shared / n_union, 4)]
pp[, denovo_shared_vs_union := round(pct_denovo_among_shared / pct_denovo_among_union, 2)]
setorder(pp, -n_shared_denovo)
cat("Host pairs ranked by number of shared thiC+ clusters (top 25):\n"); print(pp[1:min(25, .N)])
fwrite(pp, file.path(OUT_DIR, "66_hostpair_shared_producers.csv"))
cat(sprintf("\nPairs sharing >= 1 thiC+ cluster: %d of %d | pairs sharing clusters but no thiC+ cluster: %d\n",
            pp[n_shared_denovo > 0, .N], nrow(pp), pp[n_shared > 0 & n_shared_denovo == 0, .N]))
cat("Hosts that share no thiC+ cluster with any other host:\n")
sharing_hosts <- unique(c(pp[n_shared_denovo > 0, host_a], pp[n_shared_denovo > 0, host_b]))
print(setdiff(hosts, sharing_hosts))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  m <- rbind(pp[, .(x = host_a, y = host_b, n = n_shared_denovo)], pp[, .(x = host_b, y = host_a, n = n_shared_denovo)])
  ord <- all_h$traits[order(-n_genomes), host]
  m[, `:=`(x = factor(x, levels = ord), y = factor(y, levels = rev(ord)))]
  p <- ggplot(m, aes(x, y, fill = n)) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = 2.6) +
    scale_fill_gradient(low = "#f1f0ea", high = "#2a78d6", name = "Shared thiC+\nspecies clusters") +
    labs(x = NULL, y = NULL, title = "Thiamine-producing species shared between ruminant hosts") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(OUT_DIR, "66_hostpair_heatmap.png"), p, width = 8, height = 7, dpi = 300)
  cat("Heatmap written: 66_hostpair_heatmap.png\n")
}

# =============================================================================
# (d) The most widely shared thiamine producers
# =============================================================================
hr("(d) thiC+ species clusters shared by the most hosts")
sp <- hs[denovo_capable == TRUE, .(n_hosts = uniqueN(host), hosts = paste(sort(unique(host)), collapse = "; "),
                                   guild = guild[1]), by = genome][n_hosts >= 2]
if (nrow(tax)) sp <- merge(sp, tax[, .(genome, genus)], by = "genome", all.x = TRUE)
setorder(sp, -n_hosts)
cat(sprintf("thiC+ clusters found in >= 2 hosts: %s (complete de novo: %s)\n", fmt(nrow(sp)), fmt(sp[guild == "1_de_novo", .N])))
print(sp[1:min(20, .N)])
if (nrow(tax) && nrow(sp)) { cat("\nGenera contributing most shared thiC+ clusters:\n"); print(sp[, .N, by = genus][order(-N)][1:min(15, .N)]) }
fwrite(sp, file.path(OUT_DIR, "66_shared_producer_clusters.csv"))
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
