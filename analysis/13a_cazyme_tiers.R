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
#  SCRIPT  : 13a_cazyme_tiers.R  -  CAZy family presence (propagated calls) and tier tests
#  AUTHOR  : QRQ
#  INPUT   : annotation/dbcan_output_restored/overview.tsv
#            annotation/results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s54_cazyme_propagated/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(data.table) })

P <- function(var, default) Sys.getenv(var, unset = default)
A <- paste0(URGG_PROJECT, "/annotation")
OVERVIEW_F     <- P("OVERVIEW_F",     file.path(A, "dbcan_output_restored/overview.tsv"))
CHECKPOINT_RDS <- P("CHECKPOINT_RDS", file.path(A, "results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds"))
SHARING_F      <- P("SHARING_F",      file.path(A, "results_phase3/s33_cluster_sharing/representative_host_range.csv"))
OLD_ENRICH_F   <- P("OLD_ENRICH_F",   file.path(A, "results_phase3/s34_cazyme_core_specific/table_cazyme_enrichment_stats_with_CI.csv"))
OUT_DIR        <- P("OUT_DIR",        file.path(A, "results_phase3/s54_cazyme_propagated"))
FAMILY_MODE    <- P("FAMILY_MODE", "all_domains")
MIN_TOOLS      <- as.integer(P("MIN_TOOLS", "1"))
MIN_TOTAL_FAMILY_PREVALENCE <- as.integer(P("MIN_TOTAL_FAMILY_PREVALENCE", "30"))
MIN_GENOMES_PER_TIER        <- 20
SIM_B <- 10000
set.seed(20260930)

stopifnot(FAMILY_MODE %in% c("all_domains", "first"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(OVERVIEW_F, CHECKPOINT_RDS, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
log_con <- file(file.path(OUT_DIR, "54_run_report.txt"), open = "wt")
sink(log_con, split = TRUE)
hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",")
cat(sprintf("FAMILY_MODE = %s | MIN_TOOLS = %d | MIN_TOTAL_FAMILY_PREVALENCE = %d\n",
            FAMILY_MODE, MIN_TOOLS, MIN_TOTAL_FAMILY_PREVALENCE))

TIER_HOST   <- "Host-specific (1 host species)"
TIER_SHARED <- "Shared (2-3 host species)"
TIER_CORE   <- "Core / broad (>3 host species)"
FAM_RE <- "^(GH|GT|PL|CE|AA|CBM)[0-9]+"

hr("1  Family calls on representative genes (restored IDs)")
ov <- fread(OVERVIEW_F, sep = "\t", quote = "", colClasses = "character")
need <- c("Gene ID", "dbCAN_hmm", "dbCAN_sub", "DIAMOND", "#ofTools")
if (length(setdiff(need, names(ov)))) stop("overview.tsv lacks: ", paste(setdiff(need, names(ov)), collapse = ", "),
                                             "\nColumns: ", paste(names(ov), collapse = ", "))
setnames(ov, "Gene ID", "rep_gene_id")
cat("Columns:", paste(names(ov), collapse = " | "), "\n")
cat(sprintf("Overview rows: %s | unique gene IDs: %s | IDs with '~': %.4f%%\n", fmt(nrow(ov)),
            fmt(uniqueN(ov$rep_gene_id)), 100 * mean(grepl("~", ov$rep_gene_id, fixed = TRUE))))
if (uniqueN(ov$rep_gene_id) != nrow(ov)) stop("Duplicate gene IDs in the overview table; check the deduplication step.")
ov[, n_tools := suppressWarnings(as.integer(`#ofTools`))]
cat("#ofTools distribution:\n"); print(ov[, .N, by = n_tools][order(n_tools)])
ov <- ov[!is.na(n_tools) & n_tools >= MIN_TOOLS]
ov[, idx := .I]

parse_tool <- function(x, tool, rank) {
  x[is.na(x) | x %in% c("-", "", "N")] <- NA_character_
  s <- strsplit(x, "+", fixed = TRUE)
  s[is.na(x)] <- list(character(0))
  d <- data.table(idx = rep(seq_along(s), lengths(s)), pos = sequence(lengths(s)), tok = trimws(unlist(s)))
  if (FAMILY_MODE == "first") d <- d[pos == 1L]
  d[, tok := sub("\\(.*$", "", tok)]
  d[, fam := fifelse(grepl(FAM_RE, tok), sub(paste0("(", FAM_RE, ").*$"), "\\1", tok), NA_character_)]
  bad <- d[is.na(fam) & tok != "", .N, by = tok][order(-N)]
  list(calls = unique(d[!is.na(fam), .(idx, cazy_family = fam, tool = tool, rank = rank)]), bad = bad)
}
p1 <- parse_tool(ov$dbCAN_hmm, "dbCAN_hmm", 1L)
p2 <- parse_tool(ov$dbCAN_sub, "dbCAN_sub", 2L)
p3 <- parse_tool(ov$DIAMOND,   "DIAMOND",   3L)
bad <- rbindlist(list(p1$bad[, tool := "dbCAN_hmm"], p2$bad[, tool := "dbCAN_sub"], p3$bad[, tool := "DIAMOND"]))
cat("\nTokens rejected as non-CAZy (top 15; accessions/SLH/cohesin etc. are expected here):\n")
print(head(bad[order(-N)], 15))

calls <- rbindlist(list(p1$calls, p2$calls, p3$calls))
calls[, best := min(rank), by = idx]
calls <- unique(calls[rank == best, .(idx, cazy_family, tool)])
calls[, rep_gene_id := ov$rep_gene_id[idx]]
rep_fam <- calls[, .(rep_gene_id, cazy_family, tool)]
cat(sprintf("\nRepresentatives with >=1 family: %s of %s overview rows | rep x family pairs: %s | families: %d\n",
            fmt(uniqueN(rep_fam$rep_gene_id)), fmt(nrow(ov)), fmt(nrow(rep_fam)), uniqueN(rep_fam$cazy_family)))
cat("Calls by source tool:\n"); print(rep_fam[, .(reps = uniqueN(rep_gene_id)), by = tool])
fwrite(rep_fam, file.path(OUT_DIR, "54_rep_gene_cazy_families.tsv"), sep = "\t")

rep_sub <- if ("Substrate" %in% names(ov)) {
  s <- ov[!is.na(Substrate) & !Substrate %in% c("-", ""), .(rep_gene_id, Substrate)]
  s[, .(substrate = trimws(unlist(strsplit(Substrate, "[,;]")))), by = rep_gene_id][substrate != ""]
} else NULL
rm(ov, p1, p2, p3, calls); invisible(gc())

hr("2  Propagation through MMseqs2 clusters")
cl <- readRDS(CHECKPOINT_RDS); setDT(cl)
if (!all(c("rep_gene_id", "genome") %in% names(cl))) stop("Checkpoint lacks rep_gene_id/genome columns.")
cl <- cl[, .(rep_gene_id, genome)]
reps_all <- unique(cl$rep_gene_id)
hit <- unique(rep_fam$rep_gene_id) %chin% reps_all
cat(sprintf("Representatives with a family call found in the cluster table: %.4f%%\n", 100 * mean(hit)))
if (mean(hit) < 0.999) stop("Restored overview IDs do not match the cluster table - do not proceed.")
rm(reps_all); invisible(gc())

cl <- cl[rep_gene_id %chin% unique(rep_fam$rep_gene_id)]
setkey(cl, rep_gene_id); setkey(rep_fam, rep_gene_id)
gene_fam <- merge(cl, rep_fam[, .(rep_gene_id, cazy_family)], by = "rep_gene_id", allow.cartesian = TRUE)
prop <- gene_fam[, .(n_genes = .N), by = .(genome, cazy_family)]
cat(sprintf("Member genes with a family: %s | genome x family pairs (propagated): %s\n",
            fmt(uniqueN(gene_fam[, .(rep_gene_id, genome)])), fmt(nrow(prop))))
rm(gene_fam); invisible(gc())

rep_only <- unique(rep_fam[, .(genome = sub("~.*$", "", rep_gene_id), cazy_family)])
cat(sprintf("Genome x family pairs, representative-only (restored IDs): %s\n", fmt(nrow(rep_only))))

sub_prop <- if (!is.null(rep_sub)) {
  x <- merge(cl, rep_sub, by = "rep_gene_id", allow.cartesian = TRUE)
  unique(x[, .(genome, substrate)])
} else NULL
rm(cl); invisible(gc())

hr("3  Host-sharing tiers")
sh <- fread(SHARING_F, select = c("representative", "sharing_category"))
setnames(sh, "representative", "genome")
miss_lab <- setdiff(c(TIER_HOST, TIER_SHARED, TIER_CORE), unique(sh$sharing_category))
if (length(miss_lab)) stop("Tier labels not found: ", paste(miss_lab, collapse = "; "))
tier_tot <- sh[, .(n_genomes_total = .N), by = sharing_category]
print(tier_tot)
n_unknown <- prop[!genome %chin% sh$genome, uniqueN(genome)]
cat(sprintf("Genomes in propagated table that are not representatives: %d (must be 0)\n", n_unknown))
if (n_unknown > 0) warning("Non-representative genome IDs in the cluster table - check before using.")
prop <- prop[genome %chin% sh$genome]; rep_only <- rep_only[genome %chin% sh$genome]
cat(sprintf("Genomes with >=1 CAZy family: rep-only %s | propagated %s | of %s representatives\n",
            fmt(uniqueN(rep_only$genome)), fmt(uniqueN(prop$genome)), fmt(nrow(sh))))

tier_tests <- function(pres, label) {
  x <- merge(unique(pres[, .(genome, cazy_family)]), sh, by = "genome")
  byt <- x[, .(n_genomes_with_family = .N), by = .(cazy_family, sharing_category)]
  grid <- CJ(cazy_family = unique(byt$cazy_family), sharing_category = tier_tot$sharing_category)
  byt <- merge(grid, byt, by = c("cazy_family", "sharing_category"), all.x = TRUE)
  byt[is.na(n_genomes_with_family), n_genomes_with_family := 0L]
  byt <- merge(byt, tier_tot, by = "sharing_category")
  byt[, pct_genomes := round(100 * n_genomes_with_family / n_genomes_total, 2)]
  one <- function(f) {
    t <- byt[cazy_family == f]
    n_tot <- sum(t$n_genomes_with_family)
    res <- data.table(cazy_family = f, n_genomes = n_tot, test = NA_character_, statistic = NA_real_, p_value = NA_real_,
                      OR_core_vs_hostspecific = NA_real_, OR_ci_low = NA_real_, OR_ci_high = NA_real_, p_fisher = NA_real_)
    if (any(t$n_genomes_total < MIN_GENOMES_PER_TIER)) { res$test <- "skipped_low_n"; return(res) }
    if (n_tot < MIN_TOTAL_FAMILY_PREVALENCE)          { res$test <- "skipped_low_prevalence"; return(res) }
    m <- cbind(t$n_genomes_with_family, t$n_genomes_total - t$n_genomes_with_family)
    e <- suppressWarnings(chisq.test(m)$expected)
    ct <- if (any(e < 5)) chisq.test(m, simulate.p.value = TRUE, B = SIM_B) else chisq.test(m)
    res$test <- if (any(e < 5)) "chisq_simulated" else "chisq"
    res$statistic <- unname(ct$statistic); res$p_value <- ct$p.value
    a <- t[sharing_category == TIER_CORE]; b <- t[sharing_category == TIER_HOST]
    ft <- fisher.test(matrix(c(a$n_genomes_with_family, a$n_genomes_total - a$n_genomes_with_family,
                               b$n_genomes_with_family, b$n_genomes_total - b$n_genomes_with_family), nrow = 2, byrow = TRUE))
    res$OR_core_vs_hostspecific <- unname(ft$estimate); res$OR_ci_low <- ft$conf.int[1]
    res$OR_ci_high <- ft$conf.int[2]; res$p_fisher <- ft$p.value
    res
  }
  st <- rbindlist(lapply(unique(byt$cazy_family), one))
  tested <- st$test %in% c("chisq", "chisq_simulated")
  st[, p_adj_BH := NA_real_]; st[tested, p_adj_BH := p.adjust(p_value, "BH")]
  st[, direction := fcase(tested & p_adj_BH < 0.05 & OR_ci_low > 1, "core_broad_enriched",
                          tested & p_adj_BH < 0.05 & OR_ci_high < 1, "host_specific_enriched",
                          tested & p_adj_BH < 0.05, "omnibus_only (CI spans 1)",
                          tested, "ns", default = "not_tested")]
  wide <- dcast(byt, cazy_family ~ sharing_category, value.var = "pct_genomes")
  st <- merge(st, wide, by = "cazy_family")
  cat(sprintf("\n[%s] families detected: %d | tested: %d | BH p < 0.05: %d | core/broad-enriched: %d | host-specific-enriched: %d | omnibus only: %d\n",
              label, nrow(st), sum(tested), sum(st$p_adj_BH < 0.05, na.rm = TRUE),
              sum(st$direction == "core_broad_enriched"), sum(st$direction == "host_specific_enriched"),
              sum(st$direction == "omnibus_only (CI spans 1)")))
  list(stats = st[order(p_adj_BH)], by_tier = byt, wide = wide)
}

hr("4  Descriptive tier tests (unadjusted; defines the 'raw' tier-associated set)")
res_rep  <- tier_tests(rep_only, "rep-only, restored IDs")
res_prop <- tier_tests(prop,     "PROPAGATED (primary)")
st <- res_prop$stats

fwrite(res_prop$by_tier, file.path(OUT_DIR, "table_s54_cazyme_family_by_tier_long.csv"))
fwrite(res_prop$wide,    file.path(OUT_DIR, "table_s54_cazyme_family_by_tier_wide.csv"))
fwrite(st,               file.path(OUT_DIR, "table_s54_cazyme_enrichment_stats_with_CI.csv"))
fwrite(prop,             file.path(OUT_DIR, "54_cazy_family_presence_propagated.tsv"), sep = "\t")
fwrite(rep_only,         file.path(OUT_DIR, "54_cazy_family_presence_repOnly.tsv"), sep = "\t")
tested_fams <- st[test %in% c("chisq", "chisq_simulated"), cazy_family]
glmm_in <- prop[cazy_family %chin% tested_fams, .(genome, cazy_family)]
fwrite(glmm_in, file.path(OUT_DIR, "54_cazy_presence_for_glmm.tsv"), sep = "\t")

cat("\nOrientation / key families (propagated):\n")
key <- c("GH190", "CBM98", "GT117", "GH189", "GT121", "GT26")
print(st[cazy_family %in% key, .(cazy_family, n_genomes, `Core / broad (>3 host species)`, `Host-specific (1 host species)`,
                                 OR = round(OR_core_vs_hostspecific, 2), lo = round(OR_ci_low, 2), hi = round(OR_ci_high, 2),
                                 p_adj_BH = signif(p_adj_BH, 3), direction)])
rk <- st[direction %in% c("core_broad_enriched", "host_specific_enriched")]
rk[, sort_key := fifelse(direction == "core_broad_enriched", -OR_core_vs_hostspecific, OR_core_vs_hostspecific)]
rk[, rank_in_direction := frank(sort_key, ties.method = "min"), by = direction][, sort_key := NULL]
fwrite(rk[order(direction, rank_in_direction)], file.path(OUT_DIR, "table_s54_directional_families_ranked.csv"))
for (dd in c("core_broad_enriched", "host_specific_enriched")) {
  cat(sprintf("\nTop 15 %s families by OR (all carrier counts; see n_genomes):\n", dd))
  print(head(rk[direction == dd][order(rank_in_direction)], 15)[,
        .(rank = rank_in_direction, cazy_family, n_genomes, core = `Core / broad (>3 host species)`,
          shared = `Shared (2-3 host species)`, host = `Host-specific (1 host species)`,
          OR = round(OR_core_vs_hostspecific, 2), lo = round(OR_ci_low, 2), hi = round(OR_ci_high, 2))])
}
cat("\nRank of key families within their direction:\n")
print(rk[cazy_family %in% key, .(cazy_family, direction, rank_in_direction, n_genomes,
                                 OR = round(OR_core_vs_hostspecific, 2))][order(direction, rank_in_direction)])
cat(sprintf("(core/broad-enriched: %d families; host-specific-enriched: %d)\n",
            rk[direction == "core_broad_enriched", .N], rk[direction == "host_specific_enriched", .N]))

cat("\nTop 20 families by adjusted p (propagated):\n")
print(head(st[, .(cazy_family, n_genomes, core = `Core / broad (>3 host species)`, host = `Host-specific (1 host species)`,
                  OR = round(OR_core_vs_hostspecific, 2), lo = round(OR_ci_low, 2), hi = round(OR_ci_high, 2),
                  p_adj_BH = signif(p_adj_BH, 3), direction)], 20))

hr("5  Comparison with the s34 tier table")
if (file.exists(OLD_ENRICH_F)) {
  old <- fread(OLD_ENRICH_F)
  old <- old[, .(cazy_family, old_test = test, old_n_genomes = n_genomes, old_p_adj_BH = p_adj_BH,
                 old_OR = OR_core_vs_hostspecific, old_lo = OR_ci_low, old_hi = OR_ci_high)]
  old[, old_direction := fcase(old_p_adj_BH < 0.05 & old_OR > 1, "core_broad_enriched",
                               old_p_adj_BH < 0.05 & old_OR < 1, "host_specific_enriched",
                               old_test %in% c("chisq", "chisq_simulated"), "ns", default = "not_tested")]
  cmp <- merge(old, st[, .(cazy_family, n_genomes, test, p_adj_BH, OR = OR_core_vs_hostspecific, OR_ci_low, OR_ci_high, direction)],
               by = "cazy_family", all = TRUE)
  cmp[is.na(old_direction), old_direction := "absent"]; cmp[is.na(direction), direction := "absent"]
  fwrite(cmp, file.path(OUT_DIR, "table_s54_old_vs_new_families.csv"))
  o_valid <- cmp[grepl(FAM_RE, cazy_family)]
  cat(sprintf("OLD: detected %d (valid names %d) | tested %d | significant %d (core %d, host %d)\n",
              nrow(old), sum(grepl(FAM_RE, old$cazy_family)), sum(old$old_test %in% c("chisq", "chisq_simulated")),
              sum(old$old_p_adj_BH < 0.05, na.rm = TRUE), sum(old$old_direction == "core_broad_enriched"),
              sum(old$old_direction == "host_specific_enriched")))
  cat(sprintf("NEW: detected %d | tested %d | significant %d (core %d, host %d, omnibus only %d)\n",
              nrow(st), length(tested_fams), sum(st$p_adj_BH < 0.05, na.rm = TRUE), sum(st$direction == "core_broad_enriched"),
              sum(st$direction == "host_specific_enriched"), sum(st$direction == "omnibus_only (CI spans 1)")))
  cat("\nOld direction (rows) x new direction (columns):\n")
  print(o_valid[, table(old_direction, direction)])
  cat("\nFamilies that REVERSED direction:\n")
  print(cmp[(old_direction == "core_broad_enriched" & direction == "host_specific_enriched") |
            (old_direction == "host_specific_enriched" & direction == "core_broad_enriched"),
            .(cazy_family, old_n_genomes, n_genomes, old_OR = round(old_OR, 2), OR = round(OR, 2))])
  cat("\nLargest change in carrier genomes (new - old), top 10:\n")
  print(head(cmp[!is.na(old_n_genomes) & !is.na(n_genomes)][, d := n_genomes - old_n_genomes][order(-abs(d)),
             .(cazy_family, old_n_genomes, n_genomes, d, old_OR = round(old_OR, 2), OR = round(OR, 2))], 10))
} else cat("Old table not found at ", OLD_ENRICH_F, " - comparison skipped.\n", sep = "")

if (!is.null(sub_prop)) {
  hr("6  Substrate carriage by tier (propagated)")
  sub_prop <- sub_prop[genome %chin% sh$genome & !substrate %in% c("-", "", "N/A", "NA")]
  fwrite(sub_prop, file.path(OUT_DIR, "54_cazy_substrate_presence_propagated.tsv"), sep = "\t")
  s <- merge(sub_prop, sh, by = "genome")
  sb <- s[, .(n_genomes = .N), by = .(substrate, sharing_category)]
  sb <- merge(sb, tier_tot, by = "sharing_category")
  sb[, pct_genomes := round(100 * n_genomes / n_genomes_total, 2)]
  fwrite(sb[order(sharing_category, -pct_genomes)], file.path(OUT_DIR, "table_s54_substrate_by_tier.csv"))
  sres <- tier_tests(sub_prop[, .(genome, cazy_family = substrate)], "substrates")$stats
  setnames(sres, "cazy_family", "substrate")
  fwrite(sres, file.path(OUT_DIR, "table_s54_substrate_tier_tests.csv"))
  options(width = 200)
  print(sres[order(-`Core / broad (>3 host species)`),
             .(substrate, n_genomes, core = `Core / broad (>3 host species)`, shared = `Shared (2-3 host species)`,
               host = `Host-specific (1 host species)`, OR = round(OR_core_vs_hostspecific, 2), lo = round(OR_ci_low, 2),
               hi = round(OR_ci_high, 2), p_adj_BH = signif(p_adj_BH, 3), direction)])
}

hr("DONE")
cat("Outputs in: ", OUT_DIR, "\n", sep = "")
cat("\nTaxonomy-adjusted models for the tested families (15a):\n")
lab49 <- if (MIN_TOOLS == 1 && FAMILY_MODE == "all_domains") "cazy_prop" else sprintf("cazy_prop_min%dtools_%s", MIN_TOOLS, FAMILY_MODE)
cat(sprintf("  GENOME_COL=genome TRAIT_COL=cazy_family N_CORES=24 ADJ_KEEP=unadjusted,full \\\n    Rscript 15a_trait_presence_glmm.R %s %s 2>&1 | tee 15a_%s.log\n",
            file.path(OUT_DIR, "54_cazy_presence_for_glmm.tsv"), lab49, lab49))
sink(); close(log_con)
