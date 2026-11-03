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
#  SCRIPT  : 79_castor_coop_categories.R  -  Gain/loss rates of cooperative-trait categories (castor)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s41_castor_gain_loss/catalogue_16382_midpointrooted.tree
#            annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_category_presence_per_genome.csv
#  OUTPUT  : annotation/results_phase3/s79_castor_coop_categories/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(ape); library(castor); library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
R3 <- paste0(URGG_PROJECT, "/annotation/results_phase3")
ROOTED_TREE_F <- P("ROOTED_TREE_F", file.path(R3, "s41_castor_gain_loss/catalogue_16382_midpointrooted.tree"))
PRES_F  <- P("PRES_F", file.path(R3, "s36_cooperative_ko_core_specific/table_cooperative_category_presence_per_genome.csv"))
GENOME_COL <- P("GENOME_COL", "genome"); CAT_COL <- P("CAT_COL", "cooperative_categories")
OLD_F   <- P("OLD_F", file.path(R3, "s42_castor_gain_loss_scaled/table_castor_gain_loss_scaled.csv"))
OUT_DIR <- P("OUT_DIR", file.path(R3, "s79_castor_coop_categories"))
RATE_CAP <- 1000
for (f in c(ROOTED_TREE_F, PRES_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "79_run_report.txt"), split = TRUE); options(width = 200)
cat("castor", as.character(packageVersion("castor")), "| ape", as.character(packageVersion("ape")), "|", R.version.string, "\n")

tr <- read.tree(ROOTED_TREE_F); tips <- tr$tip.label
cat(sprintf("Tree: %d tips, rooted = %s\n", Ntip(tr), is.rooted(tr)))

pr <- fread(PRES_F)
if (!all(c(GENOME_COL, CAT_COL) %in% names(pr))) stop("PRES_F lacks ", GENOME_COL, "/", CAT_COL, "; columns: ", paste(names(pr), collapse = ", "))
flag <- setdiff(names(pr)[vapply(pr, function(v) is.logical(v) || all(toupper(as.character(na.omit(v))) %in% c("TRUE", "FALSE")), logical(1))], c(GENOME_COL, CAT_COL))
if (length(flag)) { cat("Presence flag column:", flag[1], "\n"); pr <- pr[toupper(as.character(get(flag[1]))) == "TRUE"] } else
  cat("No TRUE/FALSE column: every row counted as presence\n")
pr <- unique(pr[, .(genome = get(GENOME_COL), category = get(CAT_COL))])
cat(sprintf("Genomes with >= 1 category on the tree: %d of %d tips\n", pr[genome %chin% tips, uniqueN(genome)], length(tips)))

run_asr <- function(cat_name) {
  present <- pr[category == cat_name, genome]
  tip_states <- ifelse(tips %in% present, 2L, 1L)
  n1 <- sum(tip_states == 2L); n0 <- sum(tip_states == 1L)
  base <- data.table(category = cat_name, n_present = n1, n_absent = n0, prevalence = n1 / length(tips))
  if (n1 < 5 || n0 < 5) return(base[, status := "skipped_low_n"])
  set.seed(42)
  fit <- tryCatch(asr_mk_model(tr, tip_states, Nstates = 2, rate_model = "ARD", include_ancestral_likelihoods = FALSE, Ntrials = 10),
                  error = function(e) e)
  escalated <- FALSE
  if (!inherits(fit, "error") && max(fit$transition_matrix[1, 2], fit$transition_matrix[2, 1]) > RATE_CAP) {
    escalated <- TRUE
    fit2 <- tryCatch(asr_mk_model(tr, tip_states, Nstates = 2, rate_model = "ARD", include_ancestral_likelihoods = FALSE, Ntrials = 50),
                     error = function(e) e)
    if (!inherits(fit2, "error") && fit2$loglikelihood > fit$loglikelihood) fit <- fit2
  }
  if (inherits(fit, "error")) return(base[, status := "fit_failed"])
  g <- fit$transition_matrix[1, 2]; l <- fit$transition_matrix[2, 1]
  base[, `:=`(gain_rate = g, loss_rate = l, loglik = fit$loglikelihood, escalated = escalated,
              status = if (max(g, l) > RATE_CAP) "high_rate_check" else "ok")]
}
res <- rbindlist(lapply(sort(unique(pr$category)), run_asr), fill = TRUE)
res[, `:=`(loss_to_gain = loss_rate / gain_rate, expected_from_prevalence = (1 - prevalence) / prevalence)]
res[, fold_departure := loss_to_gain / expected_from_prevalence]
if (file.exists(OLD_F)) {
  old <- fread(OLD_F)
  tcol <- intersect(c("trait", "category_name"), names(old))[1]; rcol <- intersect(c("loss_to_gain_ratio", "loss_to_gain"), names(old))[1]
  if (!is.na(tcol) && !is.na(rcol)) {
    o <- old[tolower(get(tcol)) %chin% tolower(res$category), .(category = get(tcol), old_loss_to_gain = get(rcol))]
    res[, key := tolower(category)]
    res <- merge(res, unique(o[, .(key = tolower(category), old_loss_to_gain)], by = "key"), by = "key", all.x = TRUE)[, key := NULL]
  }
}
res[, prevalence := round(100 * prevalence, 2)]
fwrite(res, file.path(OUT_DIR, "table_s79_castor_coop_categories.csv"))
cat("\nCooperative categories (prevalence in % of tree tips):\n")
print(res[, intersect(c("category", "n_present", "prevalence", "gain_rate", "loss_rate", "loss_to_gain", "expected_from_prevalence",
                        "fold_departure", "old_loss_to_gain", "status"), names(res)), with = FALSE], digits = 3)
cat("\nv14 Table 17 (castor 1.7.10): antibiotic degradation / siderophore / biofilm 0.18-0.38; quorum sensing 22.5; extracellular 0.02\n")
cat("Read: fold_departure near 1 = loss:gain explained by prevalence alone.\n")
sink()
