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
#  SCRIPT  : 33a_manuscript_numbers.R  -  Collect the numbers reported in the text
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/outputs of 13a-21a
#  OUTPUT  : 61_manuscript_numbers.txt
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
R3  <- paste0(URGG_PROJECT, "/annotation/results_phase3")
S49 <- P("S49", file.path(R3, "s49_trait_presence_taxonomy_glmm"))
S47 <- P("S47", file.path(R3, "s47_cooperative_taxonomy_glmm"))
S54 <- P("S54", file.path(R3, "s54_cazyme_propagated"))
SHARING_F <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
OUT_F <- P("OUT_F", "61_manuscript_numbers.txt")
NAMED <- c("GT19", "GT3", "CE11", "GT117", "CBM98", "CE2", "GT30", "GH189", "GH190", "GH123", "GH1", "GH38", "GH112", "CE5")
sink(OUT_F, split = TRUE)
options(width = 250, datatable.print.nrows = 500, datatable.print.class = FALSE)
hr <- function(t) cat("\n", strrep("=", 90), "\n", t, "\n", strrep("=", 90), "\n", sep = "")
sg <- function(d) { n <- names(d)[sapply(d, is.numeric)]; d[, (n) := lapply(.SD, signif, 3), .SDcols = n]; d[] }

# ---------------- A: all 49 tables ----------------
hr("A  49 presence models (raw p = uncorrected; BH = corrected within model level)")
for (f in sort(list.files(S49, pattern = "^table_s49_.*_presence_models\\.csv$", full.names = TRUE))) {
  d <- fread(f)
  lab <- sub("^table_s49_(.*)_presence_models\\.csv$", "\\1", basename(f))
  cat(sprintf("\n--- %s (%d traits x %d model levels) | %s ---\n", lab, uniqueN(d$trait), uniqueN(d$adjustment),
              format(file.mtime(f), "%Y-%m-%d %H:%M")))
  keep <- intersect(c("adjustment", "trait", "n_present", "pct_core", "pct_shared", "pct_host", "OR_fisher", "OR_fisher_low",
                      "OR_fisher_high", "p_fisher", "p_fisher_BH", "OR_glmm", "OR_glmm_low", "OR_glmm_high", "LRT_p", "LRT_p_BH",
                      "direction_glmm", "status"), names(d))
  x <- sg(d[, ..keep])
  if (uniqueN(d$trait) > 40) {
    x <- x[trait %in% NAMED | direction_glmm %in% c("Core-enriched", "Host-specific-enriched")]
    cat("(large table: showing named families and GLMM-significant traits only)\n")
  }
  print(x[order(adjustment, trait)])
}

# ---------------- B: cooperative pooled models ----------------
hr("B  Cooperative pooled proportion: all models (47)")
f <- file.path(S47, "table_s47_pooled_models.csv")
if (file.exists(f)) print(sg(fread(f))) else cat("Not found:", f, "\n")
f <- file.path(S47, "table_s47_category_presence_models.csv")
if (file.exists(f)) { cat("\nCategory presence models (47):\n"); print(sg(fread(f))) }

# ---------------- C: CAZy families named in text ----------------
hr("C  CAZy families named in the text (54 unadjusted + 49 adjusted)")
st <- fread(file.path(S54, "table_s54_cazyme_enrichment_stats_with_CI.csv"))
tc <- c("Core / broad (>3 host species)", "Shared (2-3 host species)", "Host-specific (1 host species)")
c1 <- st[cazy_family %in% NAMED, .(family = cazy_family, n_genomes, pct_core = get(tc[1]), pct_shared = get(tc[2]), pct_host = get(tc[3]),
                                    OR_unadj = OR_core_vs_hostspecific, lo = OR_ci_low, hi = OR_ci_high, p_unadj_BH = p_adj_BH, direction)]
g <- file.path(S49, "table_s49_cazy_prop_presence_models.csv")
if (file.exists(g)) {
  g <- fread(g)[adjustment == "adjusted for log genome size + MAG completeness"]
  g <- g[, intersect(c("trait", "OR_glmm", "OR_glmm_low", "OR_glmm_high", "LRT_p", "LRT_p_BH", "direction_glmm"), names(g)), with = FALSE]
  setnames(g, c("trait", "OR_glmm", "OR_glmm_low", "OR_glmm_high"), c("family", "OR_adj", "adj_lo", "adj_hi"), skip_absent = TRUE)
  c1 <- merge(c1, g, by = "family", all.x = TRUE)
}
print(sg(c1[order(-OR_unadj)]))

# ---------------- D: substrates ----------------
hr("D  Substrates, unadjusted (54) — all tested")
s <- fread(file.path(S54, "table_s54_substrate_tier_tests.csv"))
print(sg(s[, .(substrate, n_genomes, pct_core = get(tc[1]), pct_shared = get(tc[2]), pct_host = get(tc[3]),
               OR = OR_core_vs_hostspecific, lo = OR_ci_low, hi = OR_ci_high, p_raw = p_value, p_BH = p_adj_BH, direction)][order(-pct_core)]))

# ---------------- E: cross-host counts ----------------
hr("E  Cross-host clusters (Sankey)")
sh <- fread(SHARING_F)
cat(sprintf("Clusters with >= 2 host species: %d | of which with a genus assignment: %d\n",
            sh[n_host_species >= 2, .N], sh[n_host_species >= 2 & !is.na(genus) & genus != "", .N]))
print(sh[, .N, by = sharing_category])

# ---------------- F: versions ----------------
hr("F  Software versions")
for (p in c("glmmTMB", "TMB", "MCMCglmm", "ape", "phytools", "castor", "data.table", "ggplot2"))
  cat(sprintf("%-10s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")))
cat(R.version.string, "\n")
cat("External tools — run in the shell and paste:\n  iqtree2 --version | head -1\n  trimal --version\n  FastTree 2>&1 | head -1\n  run_dbcan version\n  emapper.py --version\n  mmseqs version\n")
sink()
cat("\nWritten: ", OUT_F, "\n", sep = "")
