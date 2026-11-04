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
#  SCRIPT  : 23a_pagel_correlated_evolution.R  -  Correlated evolution of host restriction and CAZy families (Pagel)
#  AUTHOR  : QRQ
#  INPUT   : subset tree from 19a (PAGEL_TREE_F)
#            annotation/results_phase3/s54_cazyme_propagated/
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s57_pagel_cazy/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(ape); library(phytools); library(data.table); library(parallel) })

P <- function(var, default) Sys.getenv(var, unset = default)
R3 <- paste0(URGG_PROJECT, "/annotation/results_phase3")
PAGEL_TREE_F <- P("PAGEL_TREE_F", "")
PRESENCE_F   <- P("PRESENCE_F", file.path(R3, "s54_cazyme_propagated/54_cazy_family_presence_propagated.tsv"))
SHARING_F    <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
OUT_DIR      <- P("OUT_DIR", file.path(R3, "s57_pagel_cazy"))
FAMILIES     <- trimws(strsplit(P("FAMILIES", "CBM98,CE2,GH123,GT30,GH38,GH1,GH190,GT117"), ",")[[1]])
NBOOT        <- as.integer(P("NBOOT", "1000"))
NSIM         <- as.integer(P("NSIM", "100"))
N_CORES      <- as.integer(P("N_CORES", as.character(min(24, max(1, detectCores() - 2)))))
MODE         <- P("MODE", "full")
MIN_TREE_CARRIERS <- as.integer(P("MIN_TREE_CARRIERS", "10"))
STATS_F      <- P("STATS_F", file.path(R3, "s54_cazyme_propagated/table_s54_cazyme_enrichment_stats_with_CI.csv"))
GLMM_F       <- P("GLMM_F", file.path(R3, "s49_trait_presence_taxonomy_glmm/table_s49_cazy_prop_presence_models.csv"))
DROP_TIPS    <- trimws(strsplit(P("DROP_TIPS", ""), ",")[[1]]); DROP_TIPS <- DROP_TIPS[nzchar(DROP_TIPS)]

if (!nzchar(PAGEL_TREE_F) || !file.exists(PAGEL_TREE_F)) stop("Set PAGEL_TREE_F to the bac120 IQ-TREE .treefile used for the Pagel analysis.")
for (f in c(PRESENCE_F, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG"); set.seed(20260930)
hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("phytools %s | families: %s | NBOOT = %d | NSIM = %d | cores = %d\n",
            packageVersion("phytools"), paste(FAMILIES, collapse = ", "), NBOOT, NSIM, N_CORES))

hr("1  Tree and data")
tr <- read.tree(PAGEL_TREE_F)
tr$node.label <- NULL
cat(sprintf("Tree tips: %d | rooted: %s | zero-length branches: %d\n", Ntip(tr), is.rooted(tr), sum(tr$edge.length == 0)))
sh <- fread(SHARING_F)
setnames(sh, "representative", "genome")
strip <- function(x) sub("\\.(fa|fna|fasta|faa)(\\.gz)?$", "", x)
if (mean(tr$tip.label %chin% sh$genome) < 0.9 && mean(strip(tr$tip.label) %chin% sh$genome) > 0.9) tr$tip.label <- strip(tr$tip.label)
miss <- setdiff(tr$tip.label, sh$genome)
if (length(miss)) { cat(sprintf("Dropping %d tips that are not catalogue representatives: %s\n", length(miss), paste(head(miss, 5), collapse = ", "))); tr <- drop.tip(tr, miss) }
if (length(DROP_TIPS)) { cat(sprintf("Dropping %d tips listed in DROP_TIPS\n", length(intersect(DROP_TIPS, tr$tip.label)))); tr <- drop.tip(tr, intersect(DROP_TIPS, tr$tip.label)) }
tips <- tr$tip.label
sh <- sh[genome %chin% tips]
cat(sprintf("Tips analysed: %d\n", length(tips)))
print(sh[, .N, by = sharing_category])
if (sum(tr$edge.length == 0)) tr$edge.length[tr$edge.length == 0] <- 1e-8

pres <- unique(fread(PRESENCE_F)[genome %chin% tips, .(genome, cazy_family)])
no_cazy <- setdiff(tips, pres$genome)
cat(sprintf("Tips with no CAZy family at all: %d %s\n", length(no_cazy), if (length(no_cazy)) paste0("(", paste(head(no_cazy, 5), collapse = ", "), ")") else ""))

x <- setNames(ifelse(sh$sharing_category[match(tips, sh$genome)] == "Host-specific (1 host species)", "1", "0"), tips)
cat(sprintf("Host restriction: %d host-specific (1), %d shared or core/broad (0)\n", sum(x == "1"), sum(x == "0")))

getv <- function(o, nm) { for (n in nm) if (!is.null(o[[n]])) return(o[[n]]); NA }
fit_once <- function(tree, x, y, pi_arg) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(if (is.null(pi_arg)) fitPagel(tree, x, y, model = "ARD", dep.var = "xy")
             else fitPagel(tree, x, y, model = "ARD", dep.var = "xy", pi = pi_arg), error = function(e) e),
    warning = function(v) { w <<- c(w, conditionMessage(v)); invokeRestart("muffleWarning") })
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f)))
  li <- getv(f, c("independent.logL", "logL.indep")); ld <- getv(f, c("dependent.logL", "logL.dep"))
  ok <- is.finite(li) && is.finite(ld) && (ld - li) > -1e-6 && !any(grepl("converg", w, ignore.case = TRUE))
  list(ok = ok, fit = f, li = li, ld = ld, LR = 2 * (ld - li), msg = paste(unique(w), collapse = " | "))
}
fit_pagel <- function(tree, x, y) {
  a <- fit_once(tree, x, y, NULL)
  if (a$ok) return(c(a, pi = "equal"))
  b <- fit_once(tree, x, y, "fitzjohn")
  if (b$ok) return(c(b, pi = "fitzjohn"))
  c(a, pi = "failed")
}
rates <- function(Q) {
  s <- rownames(Q)
  q <- function(from, to) if (all(c(from, to) %in% s)) Q[from, to] else NA_real_
  data.table(gain_broader = q("0|0", "0|1"), loss_broader = q("0|1", "0|0"),
             gain_hostrestricted = q("1|0", "1|1"), loss_hostrestricted = q("1|1", "1|0"),
             to_hostrestricted_absent = q("0|0", "1|0"), to_hostrestricted_present = q("0|1", "1|1"),
             to_broader_absent = q("1|0", "0|0"), to_broader_present = q("1|1", "0|1"))
}

one_family <- function(fam) {
  t0 <- Sys.time()
  y <- setNames(ifelse(tips %in% pres[cazy_family == fam, genome], "1", "0"), tips)
  row <- data.table(family = fam, n_tips = length(tips), n_present = sum(y == "1"),
                    pct_present_host = round(100 * mean(y[x == "1"] == "1"), 1),
                    pct_present_broader = round(100 * mean(y[x == "0"] == "1"), 1))
  if (min(table(factor(y, levels = c("0", "1")))) < 2) return(list(row = row[, status := "no_variation"]))
  obs <- fit_pagel(tr, x, y)
  if (obs$pi == "failed") return(list(row = row[, `:=`(status = "fit_failed", note = obs$msg)]))
  Qi <- getv(obs$fit, c("independent.Q")); Qd <- getv(obs$fit, c("dependent.Q"))
  # ---- parametric bootstrap under the independent model ----
  boot <- unlist(mclapply(seq_len(NBOOT), function(b) {
    sim <- tryCatch(as.character(sim.Mk(tr, Qi)), error = function(e) NULL)
    if (is.null(sim)) return(NA_real_)
    xs <- setNames(sub("\\|.*$", "", sim), tips); ys <- setNames(sub("^.*\\|", "", sim), tips)
    if (length(unique(xs)) < 2 || length(unique(ys)) < 2) return(NA_real_)
    r <- fit_pagel(tr, xs, ys)
    if (r$pi == "failed") NA_real_ else r$LR
  }, mc.cores = N_CORES, mc.set.seed = TRUE))
  n_valid <- sum(is.finite(boot))
  p_boot <- (1 + sum(boot[is.finite(boot)] >= obs$LR)) / (n_valid + 1)
  # ---- stochastic character mapping of the family ----
  sm <- tryCatch(make.simmap(tr, factor(y, levels = c("0", "1")), model = "ARD", nsim = NSIM, message = FALSE), error = function(e) e)
  if (!inherits(sm, "error")) {
    cnt <- describe.simmap(sm, plot = FALSE)$count
    gcol <- grep("^0,1$", colnames(cnt)); lcol <- grep("^1,0$", colnames(cnt))
    gains <- cnt[, gcol]; losses <- cnt[, lcol]
    sm_row <- data.table(simmap_mean_gains = mean(gains), simmap_gain_range = paste(range(gains), collapse = "-"),
                         simmap_mean_losses = mean(losses), simmap_loss_range = paste(range(losses), collapse = "-"),
                         simmap_pct_no_gain = 100 * mean(gains == 0))
  } else sm_row <- data.table(simmap_mean_gains = NA_real_, simmap_note = conditionMessage(sm))
  row <- cbind(row, data.table(root_prior = obs$pi, logL_indep = obs$li, logL_dep = obs$ld, LR = obs$LR,
                               p_chisq4_reference = pchisq(obs$LR, df = 4, lower.tail = FALSE),
                               n_valid = n_valid, p_boot = p_boot), rates(Qd), sm_row,
               data.table(minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), status = "ok",
                          fit_warnings = obs$msg))
  saveRDS(list(fit = obs$fit, boot_LR = boot, simmap_counts = if (!inherits(sm, "error")) cnt else NULL),
          file.path(OUT_DIR, sprintf("57_pagel_%s.rds", fam)))
  cat(sprintf("  %-6s LR = %6.2f | p_boot = %.4f (n_valid %d) | prior %s | %.1f min\n",
              fam, obs$LR, p_boot, n_valid, obs$pi, row$minutes))
  list(row = row)
}

if (MODE == "screen") {
  hr("2  SCREEN: one fitPagel per family (no bootstrap) - exploratory")
  cnt <- pres[, .(n_present = .N), by = cazy_family][n_present >= MIN_TREE_CARRIERS & (length(tips) - n_present) >= MIN_TREE_CARRIERS]
  cat(sprintf("Families with >= %d carriers and non-carriers on the tree: %d (of %d present)\n",
              MIN_TREE_CARRIERS, nrow(cnt), uniqueN(pres$cazy_family)))
  t0 <- Sys.time()
  scr <- rbindlist(mclapply(cnt$cazy_family, function(fam) {
    y <- setNames(ifelse(tips %in% pres[cazy_family == fam, genome], "1", "0"), tips)
    o <- fit_pagel(tr, x, y)
    row <- data.table(family = fam, n_present = sum(y == "1"),
                      pct_present_broader = round(100 * mean(y[x == "0"] == "1"), 1),
                      pct_present_host = round(100 * mean(y[x == "1"] == "1"), 1), root_prior = o$pi)
    if (o$pi == "failed") return(row)
    cbind(row, data.table(LR = o$LR, p_chisq4 = pchisq(o$LR, df = 4, lower.tail = FALSE)), rates(getv(o$fit, "dependent.Q")))
  }, mc.cores = N_CORES), fill = TRUE)
  scr[, p_chisq4_BH := p.adjust(p_chisq4, "BH")]
  if (file.exists(STATS_F)) scr <- merge(scr, fread(STATS_F)[, .(family = cazy_family, OR_unadjusted = OR_core_vs_hostspecific,
                                                                  direction_unadjusted = direction)], by = "family", all.x = TRUE)
  if (file.exists(GLMM_F)) scr <- merge(scr, fread(GLMM_F)[adjustment == "adjusted for log genome size + MAG completeness",
                                                            .(family = trait, OR_adjusted = OR_glmm, direction_adjusted = direction_glmm)],
                                        by = "family", all.x = TRUE)
  setorder(scr, p_chisq4)
  fwrite(scr, file.path(OUT_DIR, "table_s57_pagel_screen.csv"))
  cat(sprintf("Done in %.1f min. Fits failed: %d\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), scr[root_prior == "failed", .N]))
  cat(sprintf("BH p < 0.05 (asymptotic, anti-conservative): %d of %d\n", scr[p_chisq4_BH < 0.05, .N], nrow(scr)))
  options(width = 220)
  if ("direction_unadjusted" %in% names(scr)) {
    cat("\nScreen-positive (BH < 0.05) by tier-association status:\n")
    print(scr[, .(n = .N, screen_positive = sum(p_chisq4_BH < 0.05, na.rm = TRUE)),
              by = .(direction_unadjusted, direction_adjusted = if ("direction_adjusted" %in% names(scr)) direction_adjusted else NA)])
  }
  cat("\nTop 30 by LR:\n")
  print(head(scr[, intersect(c("family", "n_present", "pct_present_broader", "pct_present_host", "root_prior", "LR", "p_chisq4",
                               "p_chisq4_BH", "OR_unadjusted", "direction_unadjusted", "OR_adjusted", "direction_adjusted"), names(scr)),
                 with = FALSE][order(-LR)], 30))
  cat("\nCandidates for the full bootstrap test = screen-positive families not already in FAMILIES.\n")
  quit(save = "no")
}

hr("2  Pagel tests")
res <- rbindlist(lapply(FAMILIES, function(f) one_family(f)$row), fill = TRUE)
res[, p_boot_BH := p.adjust(p_boot, "BH")]
fwrite(res, file.path(OUT_DIR, "table_s57_pagel_summary.csv"))
options(width = 220)
hr("3  Summary")
print(res[, .(family, n_present, pct_present_broader, pct_present_host, root_prior, LR = round(LR, 2), p_boot = signif(p_boot, 3),
              p_boot_BH = signif(p_boot_BH, 3), n_valid, gain_broader = signif(gain_broader, 3), loss_broader = signif(loss_broader, 3),
              gain_hostrestricted = signif(gain_hostrestricted, 3), loss_hostrestricted = signif(loss_hostrestricted, 3))])
cat("\nStochastic mapping (family gains/losses across the tree):\n")
print(res[, intersect(c("family", "simmap_mean_gains", "simmap_gain_range", "simmap_mean_losses", "simmap_loss_range", "simmap_pct_no_gain"), names(res)), with = FALSE])
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_s57.txt"))
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
