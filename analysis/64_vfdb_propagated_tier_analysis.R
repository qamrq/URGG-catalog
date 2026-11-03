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
#  SCRIPT  : 64_vfdb_propagated_tier_analysis.R  -  Virulence factors per genome (propagated) and tier tests
#  AUTHOR  : QRQ
#  INPUT   : annotation/diamond_vfdb/diamond_vfdb_filtered.tsv
#            annotation/results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s64_vfdb_propagated/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
A  <- paste0(URGG_PROJECT, "/annotation")
R3 <- file.path(A, "results_phase3")
VF_F           <- P("VF_F", file.path(A, "diamond_vfdb/diamond_vfdb_filtered.tsv"))
CHECKPOINT_RDS <- P("CHECKPOINT_RDS", file.path(R3, "s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds"))
SHARING_F      <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
OUT_DIR        <- P("OUT_DIR", file.path(R3, "s64_vfdb_propagated"))
MIN_ID         <- as.numeric(P("MIN_ID", "80"))
MAX_E          <- as.numeric(P("MAX_E", "1e-5"))
MIN_N          <- 30
SIM_B          <- 10000
set.seed(20261002)
options(width = 200)
for (f in c(VF_F, CHECKPOINT_RDS, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "64_run_report.txt"), split = TRUE)
hr  <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",")
TIER_HOST <- "Host-specific (1 host species)"; TIER_SHARED <- "Shared (2-3 host species)"; TIER_CORE <- "Core / broad (>3 host species)"

sh <- fread(SHARING_F, select = c("representative", "sharing_category"))
setnames(sh, "representative", "genome")
stopifnot(all(c(TIER_HOST, TIER_SHARED, TIER_CORE) %in% sh$sharing_category))

hr("1  VFDB hits on cluster representatives")
vf <- fread(VF_F, sep = "\t", header = FALSE, quote = "",
            col.names = c("rep_gene_id", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"))
cat(sprintf("Rows: %s | unique representative genes: %s\n", fmt(nrow(vf)), fmt(uniqueN(vf$rep_gene_id))))
vf <- vf[pident >= MIN_ID & evalue <= MAX_E]
vf <- vf[order(rep_gene_id, -bitscore)][, .SD[1], by = rep_gene_id]
vf[, vfg := sub("^(VFG[0-9]+).*$", "\\1", sseqid)]
m <- regmatches(vf$stitle, regexec("\\[([^][]+) \\((VF[0-9]+)\\) - ([^][]+) \\((VFC[0-9]+)\\)\\]", vf$stitle))
vf[, vf_name     := vapply(m, function(x) if (length(x) == 5) paste0(x[2], " (", x[3], ")") else NA_character_, "")]
vf[, vf_category := vapply(m, function(x) if (length(x) == 5) paste0(x[4], " (", x[5], ")") else NA_character_, "")]
cat(sprintf("Parsed VF name for %.1f%% and VF category for %.1f%% of hits\n",
            100 * mean(!is.na(vf$vf_name)), 100 * mean(!is.na(vf$vf_category))))
if (mean(is.na(vf$vf_category)) > 0.05) {
  cat("Unparsed stitle examples (check the VFDB header format):\n"); print(head(vf[is.na(vf_category), stitle], 5))
}
fwrite(vf, file.path(OUT_DIR, "64_vf_rep_hits_parsed.tsv"), sep = "\t")

cl <- readRDS(CHECKPOINT_RDS); setDT(cl)
cl <- cl[, .(rep_gene_id, genome)]
reps_in_map <- vf$rep_gene_id %chin% cl$rep_gene_id
cat(sprintf("Hit IDs found as rep_gene_id in the cluster map: %s of %s (%.2f%%) — must be ~100%%\n",
            fmt(sum(reps_in_map)), fmt(nrow(vf)), 100 * mean(reps_in_map)))
if (mean(reps_in_map) < 0.999) {
  cat("Unmatched examples (truncated or reformatted IDs?):\n"); print(head(vf[!reps_in_map, rep_gene_id], 10))
  cat("Lengths of unmatched IDs:\n"); print(summary(nchar(vf[!reps_in_map, rep_gene_id])))
}

hr("2  Representative-only vs propagated genome-level calls")
cl <- cl[rep_gene_id %chin% vf$rep_gene_id]
invisible(gc())
prop <- merge(cl, vf[, .(rep_gene_id, vfg, vf_name, vf_category)], by = "rep_gene_id", allow.cartesian = TRUE)
old  <- vf[, .(rep_gene_id, genome = sub("~.*$", "", rep_gene_id), vfg, vf_name, vf_category)]
n_bad <- prop[!genome %chin% sh$genome, uniqueN(genome)]
cat(sprintf("Non-representative genome IDs in the propagated table: %d (must be 0)\n", n_bad))
if (old[!genome %chin% sh$genome, .N] > 0)
  warning("Rep-gene ID prefixes (text before '~') do not all match representative genome IDs; the rep-only comparison may be off.")
prop <- prop[genome %chin% sh$genome]; old <- old[genome %chin% sh$genome]
cat(sprintf("Genome x VF-gene rows: rep-only %s -> propagated %s\n",
            fmt(uniqueN(old[, .(genome, vfg)])), fmt(uniqueN(prop[, .(genome, vfg)]))))
long_prop <- unique(rbind(prop[, .(genome, level = "vf_gene", trait = vfg)],
                          prop[!is.na(vf_name), .(genome, level = "vf", trait = vf_name)],
                          prop[!is.na(vf_category), .(genome, level = "vf_category", trait = vf_category)],
                          prop[, .(genome, level = "any", trait = "any_VF")]))
long_old  <- unique(rbind(old[!is.na(vf_name), .(genome, level = "vf", trait = vf_name)],
                          old[!is.na(vf_category), .(genome, level = "vf_category", trait = vf_category)],
                          old[, .(genome, level = "any", trait = "any_VF")]))
fwrite(long_prop, file.path(OUT_DIR, "64_vf_presence_propagated.tsv"), sep = "\t")

tot <- sh[, .(n_total = .N), by = sharing_category]
car <- merge(tot, merge(long_old[trait == "any_VF"][sh, on = "genome", nomatch = 0][, .(rep_only = .N), by = sharing_category],
                        long_prop[trait == "any_VF"][sh, on = "genome", nomatch = 0][, .(propagated = .N), by = sharing_category],
                        by = "sharing_category", all = TRUE), by = "sharing_category")
car[, `:=`(pct_rep_only = round(100 * rep_only / n_total, 1), pct_propagated = round(100 * propagated / n_total, 1))]
car <- rbind(car, data.table(sharing_category = "ALL", n_total = sum(car$n_total), rep_only = sum(car$rep_only),
                             propagated = sum(car$propagated)), fill = TRUE)
car[sharing_category == "ALL", `:=`(pct_rep_only = round(100 * rep_only / n_total, 1), pct_propagated = round(100 * propagated / n_total, 1))]
cat("\nGenomes with >=1 VF hit, by tier (the draft's ~26.5% is the rep-only figure):\n"); print(car)
fwrite(car, file.path(OUT_DIR, "64_vf_carriage_old_vs_new_by_tier.csv"))

tier_tests <- function(pres, label) {
  x <- merge(unique(pres[, .(genome, trait)]), sh, by = "genome")
  byt <- x[, .(n = .N), by = .(trait, sharing_category)]
  grid <- CJ(trait = unique(byt$trait), sharing_category = tot$sharing_category)
  byt <- merge(grid, byt, by = c("trait", "sharing_category"), all.x = TRUE)[is.na(n), n := 0L]
  byt <- merge(byt, tot, by = "sharing_category")
  one <- function(tr) {
    t <- byt[trait == tr]; setkey(t, sharing_category)
    r <- data.table(trait = tr, n_genomes = sum(t$n),
                    pct_core = round(100 * t[TIER_CORE]$n / t[TIER_CORE]$n_total, 2),
                    pct_shared = round(100 * t[TIER_SHARED]$n / t[TIER_SHARED]$n_total, 2),
                    pct_host = round(100 * t[TIER_HOST]$n / t[TIER_HOST]$n_total, 2),
                    p_chisq = NA_real_, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_)
    if (r$n_genomes < MIN_N) return(r)
    mm <- cbind(t$n, t$n_total - t$n)
    e <- suppressWarnings(chisq.test(mm)$expected)
    r$p_chisq <- (if (any(e < 5)) chisq.test(mm, simulate.p.value = TRUE, B = SIM_B) else chisq.test(mm))$p.value
    ft <- fisher.test(matrix(c(t[TIER_CORE]$n, t[TIER_CORE]$n_total - t[TIER_CORE]$n,
                               t[TIER_HOST]$n, t[TIER_HOST]$n_total - t[TIER_HOST]$n), 2, byrow = TRUE))
    r$OR <- unname(ft$estimate); r$OR_lo <- ft$conf.int[1]; r$OR_hi <- ft$conf.int[2]
    r
  }
  st <- rbindlist(lapply(unique(byt$trait), one))
  st[, p_adj_BH := p.adjust(p_chisq, "BH")]
  st[, direction := fcase(is.na(p_chisq), "not tested",
                          p_adj_BH < 0.05 & OR_lo > 1, "core/broad", p_adj_BH < 0.05 & OR_hi < 1, "host-specific",
                          p_adj_BH < 0.05, "omnibus only", default = "ns")]
  cat(sprintf("[%s] traits: %d | tested: %d | core/broad: %d | host-specific: %d | omnibus only: %d\n", label, nrow(st),
              sum(!is.na(st$p_chisq)), sum(st$direction == "core/broad"), sum(st$direction == "host-specific"),
              sum(st$direction == "omnibus only")))
  st[order(p_adj_BH)]
}
hr("3  VF categories and any-VF by tier")
cat_new <- tier_tests(long_prop[level %in% c("vf_category", "any")], "propagated")
cat_old <- tier_tests(long_old[level %in% c("vf_category", "any")], "rep-only")
cats <- merge(cat_new, cat_old[, .(trait, pct_core_old = pct_core, pct_host_old = pct_host, OR_old = OR, direction_old = direction)],
              by = "trait", all.x = TRUE)
print(cats[order(p_adj_BH), .(trait, n_genomes, pct_core, pct_shared, pct_host, OR = signif(OR, 3), OR_lo = signif(OR_lo, 3),
                             OR_hi = signif(OR_hi, 3), p_adj_BH = signif(p_adj_BH, 3), direction, OR_old = signif(OR_old, 3), direction_old)])
fwrite(cats, file.path(OUT_DIR, "table_s64_vf_category_tier_tests.csv"))

hr("4  Individual VFs by tier (>= 30 carriers)")
vf_new <- tier_tests(long_prop[level == "vf"], "propagated")
vf_old <- tier_tests(long_old[level == "vf"], "rep-only")
vfs <- merge(vf_new, vf_old[, .(trait, OR_old = OR, direction_old = direction)], by = "trait", all.x = TRUE)
cat("\nDirection, rep-only (rows) vs propagated (columns):\n"); print(table(old = vfs$direction_old, new = vfs$direction, useNA = "ifany"))
cat("\nTop 15 by adjusted p (propagated):\n")
print(vfs[order(p_adj_BH)][1:min(15, .N), .(trait, n_genomes, pct_core, pct_host, OR = signif(OR, 3), p_adj_BH = signif(p_adj_BH, 3), direction, direction_old)])
fwrite(vfs, file.path(OUT_DIR, "table_s64_vf_tier_tests.csv"))

hr("5  Long tables for the taxonomy-adjusted models (49_)")
fwrite(long_prop[level %in% c("vf_category", "any"), .(genome, trait)], file.path(OUT_DIR, "64_vf_category_presence_for_glmm.tsv"), sep = "\t")
fwrite(long_prop[level == "vf" & trait %in% vfs[!is.na(p_chisq), trait], .(genome, trait)],
       file.path(OUT_DIR, "64_vf_presence_for_glmm.tsv"), sep = "\t")
cat("Next:\n")
cat(sprintf("  GENOME_COL=genome TRAIT_COL=trait N_CORES=24 ADJ_KEEP=unadjusted,full Rscript 49_trait_presence_taxonomy_glmm.R %s vf_category\n",
            file.path(OUT_DIR, "64_vf_category_presence_for_glmm.tsv")))
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
