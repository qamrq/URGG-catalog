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
#  SCRIPT  : 24a_thiamine.R  -  Thiamine guilds per genome, by tier and by host
#  AUTHOR  : QRQ
#  INPUT   : annotation/eggnog_output/nr_catalog.emapper.annotations
#            annotation/results_phase3/s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s63_thiamine_guild_propagated/, s66_thiamine_by_host/
#  USAGE   : Rscript 24a_thiamine.R [guilds hosts]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("guilds")) local({
suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
A  <- paste0(URGG_PROJECT, "/annotation")
R3 <- file.path(A, "results_phase3")
EGGNOG_F       <- P("EGGNOG_F", file.path(A, "eggnog_output/nr_catalog.emapper.annotations"))
CHECKPOINT_RDS <- P("CHECKPOINT_RDS", file.path(R3, "s36_cooperative_ko_core_specific/gene_level_cooperative_join_CHECKPOINT.rds"))
SHARING_F      <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
COMPLETENESS_F <- P("COMPLETENESS_F", paste0(URGG_PROJECT, "/phase_3_analysis/results_phase3/s36_methanogenesis/36h_genome_subpathway_completeness.tsv"))
COMPL_ID_COL   <- P("COMPL_ID_COL", "join_id")
COMPL_COL      <- P("COMPL_COL", "completeness")
GTDB_BAC_F     <- P("GTDB_BAC_F", file.path(A, "gtdbtk_output_v232/gtdbtk.bac120.summary.tsv"))
GTDB_AR_F      <- P("GTDB_AR_F", file.path(A, "gtdbtk_output_v232/gtdbtk.ar53.summary.tsv"))
LINEAGE_F      <- P("LINEAGE_F", file.path(R3, "s32_taxonomy_novelty/lineage_full_classification_table.csv"))
OUT_DIR        <- P("OUT_DIR", file.path(R3, "s63_thiamine_guild_propagated"))
THRESH         <- as.numeric(P("THRESH", "50"))
COMPL_MIN      <- as.numeric(P("COMPL_MIN", "90"))
MIN_N          <- 30
SIM_B          <- 10000
set.seed(20261002)
options(width = 200)

for (f in c(EGGNOG_F, CHECKPOINT_RDS, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "63_run_report.txt"), split = TRUE)
hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",")

TIER_HOST <- "Host-specific (1 host species)"; TIER_SHARED <- "Shared (2-3 host species)"; TIER_CORE <- "Core / broad (>3 host species)"

KO <- c(thiF = "K03148", thiS = "K03154", thiI = "K03151", thiH = "K03150", thiO = "K03153", thiG = "K03149",
        tenI = "K10810", thi4 = "K22699", thiC = "K03147", thiD = "K00941", thiE = "K00788", thiDE = "K14153",
        thiDN = "K21219", thiN_coupling = "K21220", thiL = "K00946", thiM = "K00878",
        thiN_TPK = "K00949", thiK = "K07251", thiB = "K02064", thiP = "K02063", thiQ = "K02062", thiT = "K16789")
MODULES <- list(
  M00127 = "K03148+K03154 K03151 K03150 K03149 K03147 ((K00941 K00788),K14153,K21219) K00946",
  M00895 = "K03148+K03154 K03151 K03153 K03149 K10810 K03147 K00941 K00788 K00946",
  M00896 = "K22699 K03147 ((K00941 (K00788,K21220)),K21219) K00946",
  M00899 = "K00878 ((K00941 K00788),K14153,K21219)")
MODULE_NAMES <- c(M00127 = "prokaryotic de novo, ThiH (tyrosine) thiazole",
                  M00895 = "prokaryotic de novo, ThiO (glycine) thiazole",
                  M00896 = "archaeal de novo, Thi4 thiazole",
                  M00899 = "HMP/HET salvage => TMP")

strip_outer <- function(s) {
  s <- trimws(s)
  while (startsWith(s, "(") && endsWith(s, ")")) {
    d <- 0; whole <- TRUE; ch <- strsplit(s, "")[[1]]
    for (i in seq_along(ch)) { d <- d + (ch[i] == "(") - (ch[i] == ")"); if (d == 0 && i < length(ch)) { whole <- FALSE; break } }
    if (!whole) break
    s <- trimws(substr(s, 2, nchar(s) - 1))
  }
  s
}
split_top <- function(s, sep) {
  ch <- strsplit(s, "")[[1]]; d <- 0; out <- character(); cur <- ""
  for (c in ch) {
    if (c == "(") d <- d + 1 else if (c == ")") d <- d - 1
    if (c == sep && d == 0) { out <- c(out, cur); cur <- "" } else cur <- paste0(cur, c)
  }
  trimws(c(out, cur))[nzchar(trimws(c(out, cur)))]
}
parse_def <- function(s) {
  s <- strip_outer(s)
  p <- split_top(s, ","); if (length(p) > 1) return(list(type = "or",  kids = lapply(p, parse_def)))
  p <- split_top(s, " "); if (length(p) > 1) return(list(type = "and", kids = lapply(p, parse_def)))
  s2 <- gsub("-", "+-", s, fixed = TRUE)
  p <- split_top(s2, "+")
  if (length(p) > 1) { p <- p[!startsWith(p, "-")]; return(list(type = "and", kids = lapply(p, parse_def))) }
  if (!grepl("^K[0-9]{5}$", s)) stop("Cannot parse KEGG definition token: ", s)
  list(type = "ko", ko = s)
}
eval_node <- function(n, M) switch(n$type,
  ko  = if (n$ko %in% colnames(M)) M[, n$ko] else rep(FALSE, nrow(M)),
  and = Reduce(`&`, lapply(n$kids, eval_node, M = M)),
  or  = Reduce(`|`, lapply(n$kids, eval_node, M = M)))
module_completeness <- function(def, M) {
  steps <- split_top(strip_outer(def), " ")
  ok <- sapply(steps, function(st) eval_node(parse_def(st), M))
  if (is.null(dim(ok))) ok <- matrix(ok, nrow = nrow(M))
  round(100 * rowMeans(ok), 1)
}

hr("1  Genomes and tiers")
sh <- fread(SHARING_F, select = c("representative", "sharing_category"))
setnames(sh, "representative", "genome")
stopifnot(all(c(TIER_HOST, TIER_SHARED, TIER_CORE) %in% sh$sharing_category))
print(sh[, .N, by = sharing_category])

hr("2  Thiamine KO calls (eggNOG), representative-only vs propagated")
eg <- fread(cmd = sprintf("grep -v '^##' %s | cut -f1,12", shQuote(EGGNOG_F)),
            sep = "\t", header = TRUE, quote = "", colClasses = "character")
setnames(eg, c("rep_gene_id", "KEGG_ko"))
eg <- eg[KEGG_ko != "-" & grepl(paste(KO, collapse = "|"), KEGG_ko)]
rep_ko <- unique(eg[, .(ko = sub("^ko:", "", unlist(strsplit(KEGG_ko, ",", fixed = TRUE)))), by = rep_gene_id][ko %in% KO])
rm(eg); invisible(gc())
cat(sprintf("Representative genes with >=1 thiamine KO: %s\n", fmt(uniqueN(rep_ko$rep_gene_id))))

old_gk <- unique(rep_ko[, .(genome = sub("~.*$", "", rep_gene_id), ko)])
cl <- readRDS(CHECKPOINT_RDS); setDT(cl)
cl <- cl[rep_gene_id %chin% rep_ko$rep_gene_id, .(rep_gene_id, genome)]
invisible(gc())
gk <- unique(merge(cl, rep_ko, by = "rep_gene_id", allow.cartesian = TRUE)[, .(genome, ko)])
rm(cl); invisible(gc())
n_bad <- gk[!genome %chin% sh$genome, uniqueN(genome)]
cat(sprintf("Genome x KO rows: rep-only %s -> propagated %s | non-representative genome IDs: %d (must be 0)\n",
            fmt(nrow(old_gk)), fmt(nrow(gk)), n_bad))
gk <- gk[genome %chin% sh$genome]
gk[, gene := names(KO)[match(ko, KO)]]
fwrite(gk, file.path(OUT_DIR, "63_thiamine_ko_presence_propagated.tsv"), sep = "\t")
cat("\n% of genomes carrying each KO (propagated; rep-only in brackets):\n")
pct <- data.table(gene = names(KO), ko = KO)
pct[, pct_prop := round(100 * sapply(ko, function(k) gk[ko == k, uniqueN(genome)]) / nrow(sh), 1)]
pct[, pct_rep  := round(100 * sapply(ko, function(k) old_gk[ko == k & genome %chin% sh$genome, uniqueN(genome)]) / nrow(sh), 1)]
print(pct)

M <- matrix(FALSE, nrow(sh), length(KO), dimnames = list(sh$genome, unname(KO)))
M[cbind(match(gk$genome, sh$genome), match(gk$ko, KO))] <- TRUE
H <- function(g) M[, KO[[g]]]

hr("3  KEGG module completeness (step-based)")
res <- data.table(genome = sh$genome, sharing_category = sh$sharing_category)
for (m in names(MODULES)) res[, (paste0(m, "_pct")) := module_completeness(MODULES[[m]], M)]
ar_ids <- if (file.exists(GTDB_AR_F)) fread(GTDB_AR_F, sep = "\t", quote = "", select = "user_genome")$user_genome else character()
res[, is_archaea := genome %chin% ar_ids]
res[, denovo_module_pct := pmax(M00127_pct, M00895_pct, fifelse(is_archaea, M00896_pct, 0))]
res[, denovo_call_threshold := fifelse(denovo_module_pct >= THRESH, "de novo module >= threshold", "none/incomplete")]
res[, thiazole_route := fcase(
  M[, KO[["thiH"]]] & M[, KO[["thiO"]]], "ThiH+ThiO",
  M[, KO[["thiH"]]], "ThiH (M00127, anaerobic-type)",
  M[, KO[["thiO"]]], "ThiO (M00895, aerobic-type)",
  M[, KO[["thi4"]]], "Thi4 (M00896)",
  default = "none")]
for (m in names(MODULES)) {
  v <- res[[paste0(m, "_pct")]]
  cat(sprintf("%s (%s): median %.1f%% | >=%g%%: %s genomes | 100%%: %s\n", m, MODULE_NAMES[[m]], median(v),
              THRESH, fmt(sum(v >= THRESH)), fmt(sum(v == 100))))
}
cat(sprintf("\nDe novo module (best of M00127/M00895; M00896 for archaea) >= %g%%: %s genomes\n", THRESH, fmt(sum(res$denovo_module_pct >= THRESH))))
cat("Thiazole route (ThiH-type expected for strict anaerobes):\n")
print(res[, .N, by = thiazole_route][order(-N)])

hr("4  Branch flags and guilds")
res[, `:=`(
  hmp_denovo        = H("thiC"),
  het_denovo        = (H("thiG") & (H("thiH") | H("thiO"))) | H("thi4"),
  coupling          = H("thiE") | H("thiDE") | H("thiDN") | H("thiN_coupling"),
  tmp_kinase_thiL   = H("thiL"),
  hmp_kinase        = H("thiD") | H("thiDE") | H("thiDN"),
  het_kinase_thiM   = H("thiM"),
  tpk_thiN          = H("thiN_TPK"),
  thiK_thiL         = H("thiK") & H("thiL"),
  transport         = H("thiT") | ((H("thiB") + H("thiP") + H("thiQ")) >= 2))]
res[, free_thiamine_activation := tpk_thiN | thiK_thiL]
res[, guild := fcase(
  hmp_denovo & het_denovo & coupling & tmp_kinase_thiL,         "1_de_novo",
  hmp_denovo,                                                   "2_partial_de_novo",
  coupling & (hmp_kinase | het_kinase_thiM),                    "3_precursor_salvage",
  free_thiamine_activation,                                     "4_thiamine_salvage",
  default =                                                     "5_none_detected")]
res[, old_thiC := H("thiC")][, old_K00949_labelled_thiL := H("thiN_TPK")][, old_thiE := H("thiE")]
res[, old_salvage_signature := old_K00949_labelled_thiL & !old_thiC]

if (file.exists(COMPLETENESS_F)) {
  cm <- unique(fread(COMPLETENESS_F, select = c(COMPL_ID_COL, COMPL_COL)), by = COMPL_ID_COL)
  setnames(cm, c("genome", "mag_completeness"))
  res <- merge(res, cm, by = "genome", all.x = TRUE)
  cat(sprintf("MAG completeness joined for %s of %s genomes\n", fmt(sum(!is.na(res$mag_completeness))), fmt(nrow(res))))
  cat("\nMedian MAG completeness by guild (partial de novo / none should be lower if incompleteness drives them):\n")
  print(res[, .(n = .N, median_completeness = median(mag_completeness, na.rm = TRUE)), by = guild][order(guild)])
} else { res[, mag_completeness := NA_real_]; warning("No completeness file: completeness sensitivity skipped.") }

tax <- rbindlist(lapply(c(GTDB_BAC_F, GTDB_AR_F), function(f)
  if (file.exists(f)) fread(f, sep = "\t", quote = "", select = c("user_genome", "classification")) else NULL))
if (nrow(tax)) {
  setnames(tax, c("genome", "classification"))
  tax[, phylum := sub("^p__", "", regmatches(classification, regexpr("p__[^;]*", classification)))]
  res <- merge(res, tax[, .(genome, phylum)], by = "genome", all.x = TRUE)
}
fwrite(res, file.path(OUT_DIR, "63_genome_thiamine_guilds.tsv"), sep = "\t")

cat("\nGuilds, whole catalogue:\n")
print(res[, .(n = .N, pct = round(100 * .N / nrow(res), 1)), by = guild][order(guild)])

hr("5  Guild and branch carriage by host-sharing tier")
tot <- res[, .(n_total = .N), by = sharing_category]
byt <- dcast(res[, .N, by = .(guild, sharing_category)], guild ~ sharing_category, value.var = "N", fill = 0)
byt_pct <- copy(byt)
for (t in tot$sharing_category) byt_pct[, (t) := round(100 * get(t) / tot[sharing_category == t, n_total], 1)]
cat("Guild % of genomes per tier:\n"); print(byt_pct[])
fwrite(merge(byt, byt_pct, by = "guild", suffixes = c("_n", "_pct")), file.path(OUT_DIR, "63_guild_by_tier.csv"))
cat("\nGuild composition differs by tier (chi-square, guild x tier):\n")
print(chisq.test(as.matrix(byt[, -1])))

traits_of <- function(d) {
  long <- rbind(d[, .(genome, sharing_category, trait = paste0("guild:", guild), has = TRUE)],
                melt(d[, .(genome, sharing_category, hmp_denovo, het_denovo, coupling, tmp_kinase_thiL, hmp_kinase,
                           het_kinase_thiM, tpk_thiN, free_thiamine_activation, transport,
                           denovo_module_ge_thresh = denovo_module_pct >= THRESH,
                           salvage_M00899_ge_thresh = M00899_pct >= THRESH)],
                     id.vars = c("genome", "sharing_category"), variable.name = "trait", value.name = "has")[has == TRUE])
  long[, trait := as.character(trait)]
  long
}
tier_tests <- function(d, label) {
  long <- traits_of(d); tt <- d[, .(n_total = .N), by = sharing_category]
  one <- function(tr) {
    cnt <- merge(tt, long[trait == tr, .(n = uniqueN(genome)), by = sharing_category], by = "sharing_category", all.x = TRUE)
    cnt[is.na(n), n := 0L]; setkey(cnt, sharing_category)
    r <- data.table(trait = tr, n_genomes = sum(cnt$n),
                    pct_core = round(100 * cnt[TIER_CORE]$n / cnt[TIER_CORE]$n_total, 1),
                    pct_shared = round(100 * cnt[TIER_SHARED]$n / cnt[TIER_SHARED]$n_total, 1),
                    pct_host = round(100 * cnt[TIER_HOST]$n / cnt[TIER_HOST]$n_total, 1),
                    p_chisq = NA_real_, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_)
    if (r$n_genomes < MIN_N || r$n_genomes > sum(cnt$n_total) - MIN_N) return(r)
    m <- cbind(cnt$n, cnt$n_total - cnt$n)
    e <- suppressWarnings(chisq.test(m)$expected)
    r$p_chisq <- (if (any(e < 5)) chisq.test(m, simulate.p.value = TRUE, B = SIM_B) else chisq.test(m))$p.value
    ft <- fisher.test(matrix(c(cnt[TIER_CORE]$n, cnt[TIER_CORE]$n_total - cnt[TIER_CORE]$n,
                               cnt[TIER_HOST]$n, cnt[TIER_HOST]$n_total - cnt[TIER_HOST]$n), 2, byrow = TRUE))
    r$OR <- unname(ft$estimate); r$OR_lo <- ft$conf.int[1]; r$OR_hi <- ft$conf.int[2]
    r
  }
  st <- rbindlist(lapply(sort(unique(long$trait)), one))
  st[, p_adj_BH := p.adjust(p_chisq, "BH")]
  st[, direction := fcase(is.na(p_chisq), "not tested",
                          p_adj_BH < 0.05 & OR_lo > 1, "core/broad", p_adj_BH < 0.05 & OR_hi < 1, "host-specific",
                          p_adj_BH < 0.05, "omnibus only", default = "ns")]
  cat(sprintf("\n[%s] %s genomes\n", label, fmt(nrow(d))))
  print(st[, .(trait, n_genomes, pct_core, pct_shared, pct_host, OR = signif(OR, 3), OR_lo = signif(OR_lo, 3),
               OR_hi = signif(OR_hi, 3), p_adj_BH = signif(p_adj_BH, 3), direction)])
  st
}
st_all <- tier_tests(res, "all genomes")
fwrite(st_all, file.path(OUT_DIR, "63_guild_tier_tests.csv"))
if (any(!is.na(res$mag_completeness))) {
  st_hq <- tier_tests(res[!is.na(mag_completeness) & mag_completeness >= COMPL_MIN], sprintf("MAG completeness >= %g%%", COMPL_MIN))
  fwrite(st_hq, file.path(OUT_DIR, "63_guild_tier_tests_complete90.csv"))
  cmp <- merge(st_all[, .(trait, dir_all = direction, OR_all = signif(OR, 3))],
               st_hq[, .(trait, dir_hq = direction, OR_hq = signif(OR, 3))], by = "trait")
  cat("\nDirection, all genomes vs high-completeness subset:\n"); print(cmp)
}

hr("6  Consistency: thiamine transporter carriage by guild")
tr <- res[, .(n = .N, pct_transport = round(100 * mean(transport), 1),
              pct_thiT = round(100 * mean(H("thiT")[match(genome, sh$genome)]), 1)), by = guild][order(guild)]
print(tr); fwrite(tr, file.path(OUT_DIR, "63_transport_by_guild.csv"))
cat("Expectation: salvage guilds (3, 4) carry transporters more often than de novo producers (1).\n")

hr("7  What the old calls were (thiC, thiE, K00949 mislabelled 'thiL')")
old <- res[, .N, by = .(old_salvage_signature, guild)][order(old_salvage_signature, guild)]
old[, pct_within := round(100 * N / sum(N), 1), by = old_salvage_signature]
print(old); fwrite(old, file.path(OUT_DIR, "63_old_signature_vs_guild.csv"))
cat(sprintf("Genomes with old 'thiL' (K00949 = thiN/TPK): %s | with true thiL (K00946): %s | both: %s\n",
            fmt(sum(res$old_K00949_labelled_thiL)), fmt(sum(res$tmp_kinase_thiL)),
            fmt(sum(res$old_K00949_labelled_thiL & res$tmp_kinase_thiL))))

if (file.exists(LINEAGE_F)) {
  hr("8  RUC lineage vs rest of catalogue")
  lin <- unique(fread(LINEAGE_F)$genome)
  res[, in_lineage := genome %chin% lin]
  cat(sprintf("Lineage genomes found: %d\n", sum(res$in_lineage)))
  long <- traits_of(res)
  rr <- rbindlist(lapply(sort(unique(long$trait)), function(t) {
    has <- res$genome %chin% long[trait == t, genome]
    a <- sum(has & res$in_lineage); b <- sum(!has & res$in_lineage); c <- sum(has & !res$in_lineage); d <- sum(!has & !res$in_lineage)
    ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
    data.table(trait = t, n_lineage = a, pct_lineage = round(100 * a / (a + b), 1), pct_rest = round(100 * c / (c + d), 1),
               OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = ft$p.value)
  }))
  rr[, p_adj_BH := signif(p.adjust(p, "BH"), 3)]
  print(rr[, -"p"]); fwrite(rr, file.path(OUT_DIR, "63_ruc_lineage_guilds.csv"))
}

hr("9  Long table for the taxonomy-adjusted models (15a)")
gl <- traits_of(res)[, .(genome, trait)]
fwrite(gl, file.path(OUT_DIR, "63_guild_presence_for_glmm.tsv"), sep = "\t")
cat(sprintf("Wrote %s rows (%d traits). Taxonomy-adjusted models:\n", fmt(nrow(gl)), uniqueN(gl$trait)))
cat(sprintf("  GENOME_COL=genome TRAIT_COL=trait N_CORES=24 ADJ_KEEP=unadjusted,compl,full \\\n  Rscript 15a_trait_presence_glmm.R %s thiamine_guild\n",
            file.path(OUT_DIR, "63_guild_presence_for_glmm.tsv")))
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
sink()
})

if (run_part("hosts")) local({
suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
A  <- paste0(URGG_PROJECT, "/annotation")
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

sh <- fread(SHARING_F)
setnames(sh, "representative", "genome")
if (!"host_species_list" %in% names(sh)) stop("representative_host_range.csv has no host_species_list column; columns are: ",
                                              paste(names(sh), collapse = ", "))
g <- fread(GUILD_F, select = c("genome", "guild", "hmp_denovo", "tpk_thiN", "transport", "tmp_kinase_thiL", "denovo_module_pct"))
g[, denovo_capable := hmp_denovo]
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
})
