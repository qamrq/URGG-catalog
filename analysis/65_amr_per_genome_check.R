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
#  SCRIPT  : 65_amr_per_genome_check.R  -  AMR gene presence per genome and tier tests
#  AUTHOR  : QRQ
#  INPUT   : annotation/amr_results/, amrfinderplus_merged.tsv
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s65_amr_check/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
A  <- paste0(URGG_PROJECT, "/annotation")
R3 <- file.path(A, "results_phase3")
AMR_DIR   <- P("AMR_DIR", file.path(A, "amr_results"))
MERGED_F  <- P("MERGED_F", file.path(A, "amrfinderplus_merged.tsv"))
SHARING_F <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
OUT_DIR   <- P("OUT_DIR", file.path(R3, "s65_amr_check"))
MIN_N <- 30; SIM_B <- 10000; set.seed(20261002); options(width = 200)
for (f in c(MERGED_F, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "65_run_report.txt"), split = TRUE)
hr  <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",")
TIER_HOST <- "Host-specific (1 host species)"; TIER_SHARED <- "Shared (2-3 host species)"; TIER_CORE <- "Core / broad (>3 host species)"
sh <- fread(SHARING_F, select = c("representative", "sharing_category")); setnames(sh, "representative", "genome")

hr("1  Run completeness")
if (dir.exists(AMR_DIR)) {
  fi <- data.table(path = list.files(AMR_DIR, pattern = "_amr\\.tsv$", full.names = TRUE))
  fi[, `:=`(genome = sub("_amr\\.tsv$", "", basename(path)), size = file.size(path))]
  st <- merge(sh, fi[, .(genome, size)], by = "genome", all.x = TRUE)
  st[, status := fcase(is.na(size), "no result file (not run or failed)", size == 0, "empty file (failed)", default = "ok")]
  cat(sprintf("Result files in %s: %s | representatives: %s\n", AMR_DIR, fmt(nrow(fi)), fmt(nrow(sh))))
  cat(sprintf("Result files that are NOT representatives: %s (GENOME_DIR held other genomes too)\n",
              fmt(sum(!fi$genome %chin% sh$genome))))
  print(st[, .N, by = .(status, sharing_category)][order(status, sharing_category)])
  fwrite(st, file.path(OUT_DIR, "65_amr_run_completeness.csv"))
  run_ok <- st[status == "ok", genome]
} else {
  warning("AMR_DIR not found; completeness judged from the merged table only (genomes with zero hits cannot be distinguished from genomes not run).")
  run_ok <- NULL
}

hr("2  AMR rows: per-genome files vs merged table")
mg <- fread(MERGED_F, sep = "\t", quote = "", colClasses = "character")
cat(sprintf("Merged table: %s rows | %s genomes | modified %s\n", fmt(nrow(mg)), fmt(uniqueN(mg$genome)),
            format(file.mtime(MERGED_F), "%Y-%m-%d %H:%M")))
if (exists("fi") && nrow(fi)) {
  cat(sprintf("Result files modified between %s and %s\n", format(min(file.mtime(fi$path)), "%Y-%m-%d"),
              format(max(file.mtime(fi$path)), "%Y-%m-%d")))
  cat(sprintf("Reading %s result files (header-only files = successful runs with no hits) ...\n", fmt(sum(fi$size > 0))))
  am <- rbindlist(lapply(fi[size > 0, path], function(p) {
    d <- fread(p, sep = "\t", quote = "", colClasses = "character")
    if (nrow(d)) d[, genome := sub("_amr\\.tsv$", "", basename(p))]
    d
  }), fill = TRUE)
  if (!"genome" %in% names(am)) am[, genome := character()]
  am <- am[!is.na(genome)]
  hit_files <- unique(am$genome)
  no_file <- st[status != "ok", genome]
  mg_g <- unique(mg$genome)
  mg_nofile <- mg_g[mg_g %chin% no_file]
  cat(sprintf("Genomes with >=1 row in their own result file: %s\n", fmt(length(hit_files))))
  cat(sprintf("Merged-table genomes that HAVE a result file: %s | that have NO result file: %s\n",
              fmt(sum(mg_g %chin% fi$genome)), fmt(length(mg_nofile))))
  cat(sprintf("Genomes with hits in files but absent from the merged table: %s\n", fmt(sum(!hit_files %chin% mg$genome))))
  if (length(mg_nofile)) {
    cat(sprintf(paste0("NOTE: %s merged-table genomes have no result file in amr_results/. Their rows in the merged table are\n",
                       "      evidence of a completed run, so they are counted as run, with hits from the merged table.\n"), fmt(length(mg_nofile))))
    if (length(hit_files) == 0 && length(mg_nofile) == length(no_file))
      cat("      Every result file still in amr_results/ is header-only and the missing files are exactly the merged-table\n",
          "     genomes: the files WITH hits were moved or deleted after the merge (19 Jul). Treating header-only files as\n",
          "     zero-hit runs; confirm with the checks in the reply before publishing.\n")
    am <- rbind(am, mg[genome %chin% mg_nofile], fill = TRUE)
    run_ok <- union(run_ok, mg_nofile)
    st[genome %chin% mg_nofile, status := "ok (merged table only)"]
    fwrite(st, file.path(OUT_DIR, "65_amr_run_completeness.csv"))
  }
  cat(sprintf("Genomes counted as run: %s of %s representatives\n", fmt(length(run_ok)), fmt(nrow(sh))))
} else {
  am <- mg
}
nm <- names(am)
col <- function(...) { x <- intersect(c(...), nm); if (!length(x)) stop("None of these columns found: ", paste(c(...), collapse = ", ")); x[1] }
C_TYPE <- col("Element type", "Type"); C_SYM <- col("Gene symbol", "Element symbol")
C_CLASS <- col("Class"); C_SUB <- col("Subclass"); C_METHOD <- col("Method")
cat(sprintf("Rows: %s | genomes: %s\n", fmt(nrow(am)), fmt(uniqueN(am$genome))))
print(am[, .N, by = C_TYPE])
amr <- am[get(C_TYPE) == "AMR" & genome %chin% sh$genome,
          .(genome, gene = get(C_SYM), class = get(C_CLASS), subclass = get(C_SUB), method = get(C_METHOD))]
cat("\nAMR rows by detection method (PARTIAL / PARTIAL_CONTIG_END are truncated hits, common in MAGs):\n")
print(amr[, .N, by = method][order(-N)])

hr("3  AMR carriage and tier tests")
den <- if (is.null(run_ok)) sh else sh[genome %chin% run_ok]
cat(sprintf("Denominator: %s representatives with a successful run\n", fmt(nrow(den))))
long <- unique(rbind(amr[, .(genome, level = "class", trait = class)],
                     amr[, .(genome, level = "gene", trait = gene)],
                     amr[, .(genome, level = "any", trait = "any_AMR")]))
long <- long[genome %chin% den$genome]
fwrite(long, file.path(OUT_DIR, "65_amr_presence.tsv"), sep = "\t")
tot <- den[, .(n_total = .N), by = sharing_category]
car <- merge(tot, long[trait == "any_AMR"][den, on = "genome", nomatch = 0][, .(n_any_AMR = .N), by = sharing_category],
             by = "sharing_category", all.x = TRUE)[is.na(n_any_AMR), n_any_AMR := 0L]
car[, pct := round(100 * n_any_AMR / n_total, 2)]
print(car); fwrite(car, file.path(OUT_DIR, "65_amr_carriage_by_tier.csv"))

tier_tests <- function(pres) {
  if (!nrow(pres)) { cat("No rows to test.\n"); return(data.table(trait = character(), p_chisq = numeric())) }
  x <- merge(unique(pres[, .(genome, trait)]), den, by = "genome")
  byt <- x[, .(n = .N), by = .(trait, sharing_category)]
  byt <- merge(CJ(trait = unique(byt$trait), sharing_category = tot$sharing_category), byt, all.x = TRUE)[is.na(n), n := 0L]
  byt <- merge(byt, tot, by = "sharing_category")
  st <- rbindlist(lapply(unique(byt$trait), function(tr) {
    t <- byt[trait == tr]; setkey(t, sharing_category)
    r <- data.table(trait = tr, n_genomes = sum(t$n), pct_core = round(100 * t[TIER_CORE]$n / t[TIER_CORE]$n_total, 2),
                    pct_shared = round(100 * t[TIER_SHARED]$n / t[TIER_SHARED]$n_total, 2),
                    pct_host = round(100 * t[TIER_HOST]$n / t[TIER_HOST]$n_total, 2),
                    p_chisq = NA_real_, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_)
    if (r$n_genomes < MIN_N) return(r)
    mm <- cbind(t$n, t$n_total - t$n); e <- suppressWarnings(chisq.test(mm)$expected)
    r$p_chisq <- (if (any(e < 5)) chisq.test(mm, simulate.p.value = TRUE, B = SIM_B) else chisq.test(mm))$p.value
    ft <- fisher.test(matrix(c(t[TIER_CORE]$n, t[TIER_CORE]$n_total - t[TIER_CORE]$n, t[TIER_HOST]$n, t[TIER_HOST]$n_total - t[TIER_HOST]$n), 2, byrow = TRUE))
    r$OR <- unname(ft$estimate); r$OR_lo <- ft$conf.int[1]; r$OR_hi <- ft$conf.int[2]; r }))
  st[, p_adj_BH := p.adjust(p_chisq, "BH")]
  st[, direction := fcase(is.na(p_chisq), "not tested", p_adj_BH < 0.05 & OR_lo > 1, "core/broad",
                          p_adj_BH < 0.05 & OR_hi < 1, "host-specific", p_adj_BH < 0.05, "omnibus only", default = "ns")]
  st[order(p_adj_BH)]
}
cls <- tier_tests(long[level %in% c("class", "any")])
print(cls[, .(trait, n_genomes, pct_core, pct_shared, pct_host, OR = signif(OR, 3), OR_lo = signif(OR_lo, 3),
              OR_hi = signif(OR_hi, 3), p_adj_BH = signif(p_adj_BH, 3), direction)])
fwrite(cls, file.path(OUT_DIR, "table_s65_amr_class_tier_tests.csv"))
genes <- tier_tests(long[level == "gene"])
cat("\nAMR genes carried by >= 30 genomes:\n")
if (nrow(genes)) print(genes[!is.na(p_chisq), .(trait, n_genomes, pct_core, pct_host, OR = signif(OR, 3), p_adj_BH = signif(p_adj_BH, 3), direction)])
fwrite(genes, file.path(OUT_DIR, "table_s65_amr_gene_tier_tests.csv"))

fwrite(long[level %in% c("class", "any"), .(genome, trait)], file.path(OUT_DIR, "65_amr_class_presence_for_glmm.tsv"), sep = "\t")
cat("\nNext (taxonomy-adjusted):\n")
cat(sprintf("  GENOME_COL=genome TRAIT_COL=trait N_CORES=24 ADJ_KEEP=unadjusted,full Rscript 49_trait_presence_taxonomy_glmm.R %s amr_class\n",
            file.path(OUT_DIR, "65_amr_class_presence_for_glmm.tsv")))
if (!is.null(run_ok) && length(run_ok) < nrow(sh))
  cat("NOTE: 49_ treats genomes absent from the table as AMR-negative. Representatives without a successful run are\n",
      "      counted as negative there; rerun them (05_ is resumable) before the taxonomy-adjusted models.\n")
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
