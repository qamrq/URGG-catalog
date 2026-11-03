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
#  SCRIPT  : 25a_virulence_amr_plasmids.R  -  Virulence factors, AMR genes and plasmid-borne AMR by tier
#  AUTHOR  : QRQ
#  INPUT   : annotation/diamond_vfdb/diamond_vfdb_filtered.tsv, amr_results/, amrfinderplus_merged.tsv, platon_merged.tsv
#            annotation/results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s64_vfdb_propagated/, s65_amr_check/, s78_para126/
#  USAGE   : Rscript 25a_virulence_amr_plasmids.R [vfdb amr plasmid_amr]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("vfdb")) local({
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
cat(sprintf("  GENOME_COL=genome TRAIT_COL=trait N_CORES=24 ADJ_KEEP=unadjusted,full Rscript 15a_trait_presence_glmm.R %s vf_category\n",
            file.path(OUT_DIR, "64_vf_category_presence_for_glmm.tsv")))
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
})

if (run_part("amr")) local({
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
cat(sprintf("  GENOME_COL=genome TRAIT_COL=trait N_CORES=24 ADJ_KEEP=unadjusted,full Rscript 15a_trait_presence_glmm.R %s amr_class\n",
            file.path(OUT_DIR, "65_amr_class_presence_for_glmm.tsv")))
if (!is.null(run_ok) && length(run_ok) < nrow(sh))
  cat("NOTE: 49_ treats genomes absent from the table as AMR-negative. Representatives without a successful run are\n",
      "      counted as negative there; rerun them (05_ is resumable) before the taxonomy-adjusted models.\n")
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
})

if (run_part("plasmid_amr")) local({
suppressPackageStartupMessages(library(data.table))
P  <- function(var, default) Sys.getenv(var, unset = default)
A  <- paste0(URGG_PROJECT, "/annotation")
R3 <- file.path(A, "results_phase3")
SHARING_F <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
AMR_F     <- P("AMR_F", file.path(R3, "s65_amr_check/65_amr_presence.tsv"))
MERGED_F  <- P("MERGED_F", file.path(A, "amrfinderplus_merged.tsv"))
PLATON_F  <- P("PLATON_F", file.path(A, "platon_merged.tsv"))
OUT_DIR   <- P("OUT_DIR", file.path(R3, "s78_para126"))
for (f in c(SHARING_F, AMR_F, MERGED_F, PLATON_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "78_run_report.txt"), split = TRUE); options(width = 200)
TIER_HOST <- "Host-specific (1 host species)"; TIER_CORE <- "Core / broad (>3 host species)"

sh <- fread(SHARING_F, select = c("representative", "sharing_category")); setnames(sh, "representative", "genome")
amr_g <- unique(fread(AMR_F)[trait == "any_AMR", genome])
pl <- fread(PLATON_F, sep = "\t", quote = "")
amrcol <- grep("^# ?AMRs$", names(pl), value = TRUE)[1]
mg <- fread(MERGED_F, sep = "\t", quote = "")
tcol <- intersect(c("Type", "Element type"), names(mg))[1]; ccol <- intersect(c("Contig id", "Contig ID", "contig"), names(mg))[1]
am <- unique(mg[get(tcol) == "AMR", .(genome, contig = get(ccol))])
plc <- unique(pl[, .(genome, contig = ID)])
hit <- merge(am, plc, by = c("genome", "contig"))
cat(sprintf("AMR-type gene contigs: %d | on Platon plasmid contigs: %d (in %d genomes)\n", nrow(am), nrow(hit), uniqueN(hit$genome)))
if (!nrow(hit)) cat("WARNING: no AMR contig matched a Platon contig; check contig naming:\n  AMR e.g.", am$contig[1], "\n  Platon e.g.", plc$contig[1], "\n")

sh[, `:=`(plasmid = genome %chin% pl$genome, any_AMR = genome %chin% amr_g,
          amr_on_plasmid = genome %chin% hit$genome,
          platon_amr = if (!is.na(amrcol)) genome %chin% pl[get(amrcol) > 0, genome] else NA)]
tier_row <- function(v, label) {
  a <- sh[sharing_category == TIER_CORE, get(v)]; b <- sh[sharing_category == TIER_HOST, get(v)]
  ft <- fisher.test(matrix(c(sum(a), sum(!a), sum(b), sum(!b)), 2, byrow = TRUE))
  data.table(outcome = label, n_all = sum(sh[[v]]), pct_all = round(100 * mean(sh[[v]]), 2),
             pct_core = round(100 * mean(a), 2), pct_host = round(100 * mean(b), 2),
             OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = signif(ft$p.value, 2))
}
cat("\n(a, b) Core/broad vs host-specific\n")
tt <- rbind(tier_row("plasmid", "Plasmid contig (Platon)"),
            tier_row("amr_on_plasmid", "AMR-type gene (AMRFinderPlus) on a Platon plasmid contig"),
            if (!is.na(amrcol)) tier_row("platon_amr", "Platon '# AMRs' > 0 (draft definition, for comparison)"))
print(tt)
cat("Draft: plasmid 47.6 vs 42.0%, OR 1.25 (1.03-1.52), padj 0.058; AMR on plasmid 2.42 vs 1.14%, OR 2.14 (1.05-3.95), padj 0.058\n")

cat("\n(c) AMR-type gene anywhere: plasmid-positive vs plasmid-negative genomes\n")
ft <- fisher.test(table(factor(sh$plasmid, c(TRUE, FALSE)), factor(sh$any_AMR, c(TRUE, FALSE))))
pc <- data.table(pct_amr_plasmid_pos = round(100 * mean(sh[plasmid == TRUE, any_AMR]), 2), pct_amr_plasmid_neg = round(100 * mean(sh[plasmid == FALSE, any_AMR]), 2),
                 n_plasmid_pos = sum(sh$plasmid), OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = signif(ft$p.value, 2))
print(pc); cat("Draft: 2.0 vs 0.6%, OR 3.38 (2.61-4.41), p = 2.9e-23\n")
fwrite(tt, file.path(OUT_DIR, "table_s78_plasmid_amr_by_tier.csv")); fwrite(pc, file.path(OUT_DIR, "table_s78_plasmid_vs_amr.csv"))
sink()
})
