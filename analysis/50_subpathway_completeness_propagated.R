#!/usr/bin/env Rscript
###############################################################################
# 50_subpathway_completeness_propagated.R
#
# Rebuilds the Section 3.4 metabolic sub-pathway completeness table (36h) with
# KO calls PROPAGATED from each non-redundant representative gene to every
# member gene of its MMseqs2 cluster.
#
# WHY: 36e built its genome x KO cache directly from the eggNOG annotations,
# which exist only for the 37.1M NR representative genes, and assigned each
# KO to the genome that owns the representative. The ~9.8M member genes that
# MMseqs2 collapsed into another genome's representative were never credited
# to their own genomes, so KO presence was undercounted — and not at random:
# genomes from densely sampled lineages lose the most. Propagation through the
# cluster map (the same route used for the cooperative-gene analysis, 36_/38_,
# and the thiamine guild table, 65c) credits every genome with the KOs of all
# of its genes.
#
# Marker sets and the two-stage test framework are copied VERBATIM from 36h,
# so any change in the results is due to KO counting alone.
#
# OUTPUTS (OUT_DIR)
#   50_marker_ko_presence_propagated.tsv       genome x marker-KO (long)
#   50_genome_subpathway_completeness_propagated.tsv   same layout as 36h's table
#   50_subpathway_hosttier_tests_propagated.tsv        same format as 36h's tests (regenerates Table 5)
#   50_old_vs_new_prevalence.tsv               36h (unpropagated) vs propagated, side by side
#   subpathway_presence_propagated_{50,75,100}.tsv     long tables for 49_
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(stringr) })

P <- function(var, default) Sys.getenv(var, unset = default)
CHECKPOINT_RDS <- P("CHECKPOINT_RDS", "/disk252/qamar/metagenome_analysis/annotation/results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds")
EGGNOG_F       <- P("EGGNOG_F", "/disk252/qamar/metagenome_analysis/annotation/eggnog_output/nr_catalog.emapper.annotations")
S36H_DIR       <- P("S36H_DIR", "/disk252/qamar/metagenome_analysis/phase_3_analysis/results_phase3/s36_methanogenesis")
BASE_36H_F     <- P("BASE_36H_F", file.path(S36H_DIR, "36h_genome_subpathway_completeness.tsv"))
TESTS_36H_F    <- P("TESTS_36H_F", file.path(S36H_DIR, "36h_subpathway_hosttier_tests.tsv"))
OUT_DIR        <- P("OUT_DIR", "/disk252/qamar/metagenome_analysis/annotation/results_phase3/s50_subpathway_propagated")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(CHECKPOINT_RDS, EGGNOG_F, BASE_36H_F, TESTS_36H_F)) if (!file.exists(f)) stop("Missing: ", f)
set.seed(1)   # the omnibus Fisher test uses simulated p-values

hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

# ---- Marker sets: copied verbatim from 36h_subpathway_hosttier_stats.R ------
subpathways <- list(
  propionate_succinate       = c("K01847", "K01848", "K01849", "K01960", "K01026"),
  butyrate_kinase_route      = c("K00248", "K00929", "K00634"),
  butyrate_transferase_route = c("K00248", "K01034", "K01035"),
  acetate_production         = c("K00925", "K00625"),
  cobalamin_biosynthesis     = c("K00589", "K02224", "K00768", "K02233"),
  thiamine_biosynthesis      = c("K03147", "K00788", "K00949"),
  urea_utilization           = c("K01428", "K01429", "K01430"),
  ammonia_assimilation       = c("K01915", "K00265", "K00266"),
  nitrate_reduction          = c("K00370", "K00362"),
  bcaa_biosynthesis          = c("K01652", "K00053"),
  lysine_biosynthesis        = c("K01714", "K00928"),
  tryptophan_biosynthesis_SINGLE_MARKER_ONLY = c("K01656"),
  methanogenesis             = c("K00399", "K00401", "K00402")
)
all_markers <- unique(unlist(subpathways))

# =============================================================================
# 1. Representative genes carrying a marker KO
# =============================================================================
hr("1  Marker KOs on NR representative genes (eggNOG)")
eg <- fread(cmd = sprintf("grep -v '^##' %s | cut -f1,12", shQuote(EGGNOG_F)),
            sep = "\t", header = TRUE, quote = "", colClasses = "character")
setnames(eg, c("rep_gene_id", "KEGG_ko"))
eg <- eg[KEGG_ko != "-" & !is.na(KEGG_ko) & grepl(paste(all_markers, collapse = "|"), KEGG_ko)]
rep_ko <- unique(eg[, .(ko = sub("^ko:", "", unlist(strsplit(KEGG_ko, ",", fixed = TRUE)))), by = rep_gene_id][ko %in% all_markers])
rm(eg); invisible(gc())
cat(sprintf("Representative genes with >=1 marker KO: %s (%s rep-KO pairs)\n",
            format(uniqueN(rep_ko$rep_gene_id), big.mark = ","), format(nrow(rep_ko), big.mark = ",")))

# =============================================================================
# 2. Propagate to every member gene -> genome x KO presence
# =============================================================================
hr("2  Propagating through MMseqs2 cluster membership")
cl <- readRDS(CHECKPOINT_RDS); setDT(cl)
cl <- cl[rep_gene_id %chin% unique(rep_ko$rep_gene_id), .(rep_gene_id, genome)]
invisible(gc())
setkey(cl, rep_gene_id); setkey(rep_ko, rep_gene_id)
gk <- unique(merge(cl, rep_ko, by = "rep_gene_id", allow.cartesian = TRUE)[, .(genome_id = genome, ko)])
rm(cl); invisible(gc())
cat(sprintf("Genome x marker-KO presence rows: %s across %s genomes\n",
            format(nrow(gk), big.mark = ","), format(uniqueN(gk$genome_id), big.mark = ",")))
fwrite(gk, file.path(OUT_DIR, "50_marker_ko_presence_propagated.tsv"), sep = "\t")

# =============================================================================
# 3. Completeness per sub-pathway (same arithmetic as 36h)
# =============================================================================
hr("3  Sub-pathway completeness")
base <- fread(BASE_36H_F, sep = "\t", header = TRUE)
old  <- copy(base)
old_cols <- grep("_completeness$", names(base), value = TRUE)
base[, (old_cols) := NULL]

n_unmatched <- sum(!base$join_id %in% unique(gk$genome_id))
cat(sprintf("36h genomes: %d | with >=1 propagated marker KO: %d | with none: %d\n",
            nrow(base), nrow(base) - n_unmatched, n_unmatched))
if (!any(base$join_id %in% gk$genome_id)) stop("No genome IDs match between 36h table and checkpoint — check ID format.")

for (sp in names(subpathways)) {
  kos <- subpathways[[sp]]
  hits <- gk[ko %in% kos, .(n_found = uniqueN(ko)), by = genome_id]
  setnames(hits, "genome_id", "join_id")
  base <- merge(base, hits, by = "join_id", all.x = TRUE)
  base[is.na(n_found), n_found := 0L]
  base[[paste0(sp, "_completeness")]] <- base$n_found / length(kos)
  base[, n_found := NULL]
}
fwrite(base, file.path(OUT_DIR, "50_genome_subpathway_completeness_propagated.tsv"), sep = "\t")

# Sanity: propagation can only ADD genes, so completeness must never fall.
chk <- merge(old[, c("join_id", old_cols), with = FALSE], base[, c("join_id", old_cols), with = FALSE],
             by = "join_id", suffixes = c(".old", ".new"))
n_decrease <- sum(sapply(old_cols, function(c) sum(chk[[paste0(c, ".new")]] < chk[[paste0(c, ".old")]] - 1e-9)))
cat(sprintf("Genome x sub-pathway cells where completeness DECREASED after propagation: %d (must be 0)\n", n_decrease))
if (n_decrease > 0) warning("Completeness decreased for some genomes — genome IDs in the checkpoint and the 36e cache may differ. Investigate before using.")

# =============================================================================
# 4. Two-stage tests: copied verbatim from 36h (regenerates Table 5)
# =============================================================================
hr("4  Host-tier tests (36h framework, propagated data)")
TIER_HOST   <- "Host-specific (1 host species)"
TIER_SHARED <- "Shared (2-3 host species)"
TIER_CORE   <- "Core / broad (>3 host species)"
base[, sharing_category := factor(sharing_category, levels = c(TIER_HOST, TIER_SHARED, TIER_CORE))]
thresholds <- c(0.5, 0.75, 1.0)

run_tests_one_category <- function(d, col, cat_name, thresholds) {
  rbindlist(lapply(thresholds, function(thr) {
    d[, carries := get(col) >= thr]
    tab <- table(d$sharing_category, d$carries)
    if (ncol(tab) < 2 || nrow(tab) < 2 || nrow(d) == 0) {
      return(data.table(category = cat_name, threshold = thr, n_genomes_tested = nrow(d),
                        omnibus_test = NA_character_, omnibus_p = NA_real_, min_expected = NA_real_,
                        core_pct = NA_real_, hostspec_pct = NA_real_,
                        OR_core_vs_hostspec = NA_real_, OR_ci_low = NA_real_, OR_ci_high = NA_real_,
                        directional_p = NA_real_, note = "degenerate table"))
    }
    csq <- suppressWarnings(chisq.test(tab)); min_exp <- min(csq$expected)
    if (min_exp < 5) {
      ft <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
      omni_test <- "Fisher (simulated; expected<5)"; omni_p <- ft$p.value
    } else { omni_test <- "Chi-square"; omni_p <- csq$p.value }
    sub <- d[sharing_category %in% c(TIER_CORE, TIER_HOST)]
    sub[, tier2 := factor(as.character(sharing_category), levels = c(TIER_HOST, TIER_CORE))]
    tab2 <- table(sub$tier2, sub$carries)
    if (ncol(tab2) < 2 || nrow(tab2) < 2) { or <- ci_l <- ci_h <- dir_p <- NA_real_ } else {
      ft2 <- fisher.test(tab2); or <- unname(ft2$estimate); ci_l <- ft2$conf.int[1]; ci_h <- ft2$conf.int[2]; dir_p <- ft2$p.value
    }
    data.table(category = cat_name, threshold = thr, n_genomes_tested = nrow(d),
               omnibus_test = omni_test, omnibus_p = omni_p, min_expected = round(min_exp, 1),
               core_pct = round(100 * mean(d[sharing_category == TIER_CORE]$carries), 2),
               hostspec_pct = round(100 * mean(d[sharing_category == TIER_HOST]$carries), 2),
               OR_core_vs_hostspec = round(or, 3), OR_ci_low = round(ci_l, 3), OR_ci_high = round(ci_h, 3),
               directional_p = dir_p, note = "")
  }))
}

results <- rbindlist(lapply(names(subpathways), function(sp) {
  col <- paste0(sp, "_completeness")
  d <- if (sp == "methanogenesis") base[is_methanogen == TRUE] else copy(base)
  run_tests_one_category(d, col, sp, thresholds)
}))
results[!is.na(omnibus_p), omnibus_p_adj := p.adjust(omnibus_p, method = "BH")]
results[!is.na(directional_p), directional_p_adj := p.adjust(directional_p, method = "BH")]
results[, direction := fifelse(is.na(OR_ci_low) | is.na(OR_ci_high), "untested",
                        fifelse(OR_ci_low > 1, "core-enriched",
                        fifelse(OR_ci_high < 1, "host-specific-enriched", "no directional claim (CI spans 1)")))]
setcolorder(results, c("category", "threshold", "n_genomes_tested", "omnibus_test", "omnibus_p", "omnibus_p_adj",
                       "min_expected", "core_pct", "hostspec_pct", "OR_core_vs_hostspec", "OR_ci_low", "OR_ci_high",
                       "directional_p", "directional_p_adj", "direction", "note"))
fwrite(results, file.path(OUT_DIR, "50_subpathway_hosttier_tests_propagated.tsv"), sep = "\t")

# =============================================================================
# 5. Old vs new, side by side
# =============================================================================
hr("5  36h (unpropagated) vs propagated")
old_tests <- fread(TESTS_36H_F)
cmp <- merge(old_tests[, .(category, threshold, core_pct_old = core_pct, host_pct_old = hostspec_pct,
                           OR_old = OR_core_vs_hostspec, dir_old = direction)],
             results[, .(category, threshold, core_pct_new = core_pct, host_pct_new = hostspec_pct,
                         OR_new = OR_core_vs_hostspec, p_adj_new = directional_p_adj, dir_new = direction)],
             by = c("category", "threshold"))
cmp[, direction_changed := dir_old != dir_new]
fwrite(cmp, file.path(OUT_DIR, "50_old_vs_new_prevalence.tsv"), sep = "\t")
options(width = 220)
print(cmp[threshold == 0.5][order(category)], digits = 3)
cat(sprintf("\nDirection calls that changed (any threshold): %d of %d\n", sum(cmp$direction_changed, na.rm = TRUE), nrow(cmp)))
if (any(cmp$direction_changed, na.rm = TRUE)) print(cmp[direction_changed == TRUE])

# =============================================================================
# 6. Long presence tables for 49_ (methanogenesis and single-marker tryptophan excluded)
# =============================================================================
hr("6  Presence tables for 49_")
keep <- setdiff(names(subpathways), c("methanogenesis", "tryptophan_biosynthesis_SINGLE_MARKER_ONLY"))
long <- melt(base[, c("join_id", paste0(keep, "_completeness")), with = FALSE], id.vars = "join_id",
             variable.name = "subpathway", value.name = "compl")
long[, subpathway := sub("_completeness$", "", subpathway)]
for (t in c(0.5, 0.75, 1)) {
  f <- file.path(OUT_DIR, sprintf("subpathway_presence_propagated_%d.tsv", round(t * 100)))
  fwrite(long[compl >= t - 1e-9, .(genome_id = join_id, subpathway)], f, sep = "\t")
  cat("Written:", f, "\n")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_s50.txt"))
hr("DONE")
