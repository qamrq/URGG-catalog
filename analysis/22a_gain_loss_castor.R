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
#  SCRIPT  : 22a_gain_loss_castor.R  -  Gain/loss rates of CAZy families and cooperative categories (castor)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s41_castor_gain_loss/catalogue_16382_midpointrooted.tree
#            annotation/results_phase3/s54_cazyme_propagated/
#            annotation/results_phase3/s36_cooperative_ko_core_specific/
#  OUTPUT  : annotation/results_phase3/s56_castor_cazy_propagated/, s79_castor_coop_categories/
#  USAGE   : Rscript 22a_gain_loss_castor.R [cazy cooperative]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("cazy")) local({
suppressPackageStartupMessages({ library(ape); library(castor); library(data.table); library(parallel) })

P <- function(var, default) Sys.getenv(var, unset = default)
`%||%` <- function(a, b) if (is.null(a)) b else a
R3 <-paste0(URGG_PROJECT, "/annotation/results_phase3")
ROOTED_TREE_F <- P("ROOTED_TREE_F", file.path(R3, "s41_castor_gain_loss/catalogue_16382_midpointrooted.tree"))
S54_DIR       <- P("S54_DIR", file.path(R3, "s54_cazyme_propagated"))
PRESENCE_F    <- P("PRESENCE_F", file.path(S54_DIR, "54_cazy_family_presence_propagated.tsv"))
STATS_F       <- P("STATS_F", file.path(S54_DIR, "table_s54_cazyme_enrichment_stats_with_CI.csv"))
GLMM_F        <- P("GLMM_F", file.path(R3, "s49_trait_presence_taxonomy_glmm/table_s49_cazy_prop_presence_models.csv"))
SHARING_F     <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
OLD_CASTOR_F  <- P("OLD_CASTOR_F", file.path(R3, "s42_castor_gain_loss_scaled/table_castor_gain_loss_scaled.csv"))
OUT_DIR       <- P("OUT_DIR", file.path(R3, "s56_castor_cazy_propagated"))
N_CORES       <- as.integer(P("N_CORES", as.character(max(1, detectCores() - 2))))
RATE_CAP      <- 1000
KEY <- c("GH190", "CBM98", "GT117", "GH189", "CE2", "GH123", "GT30", "GH38", "GH1")

for (f in c(ROOTED_TREE_F, PRESENCE_F, STATS_F, SHARING_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "table_s56_castor_summary.txt"), split = TRUE)
hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

# ---------------- Tree + tier composition ----------------
hr("1  Tree")
tr <- read.tree(ROOTED_TREE_F)
tips <- tr$tip.label
cat(sprintf("Tips: %d | rooted: %s\n", Ntip(tr), is.rooted(tr)))
sh <- fread(SHARING_F, select = c("representative", "sharing_category"))
comp <- sh[representative %chin% tips, .N, by = sharing_category]
comp <- merge(comp, sh[, .(catalogue = .N), by = sharing_category], by = "sharing_category")
comp[, pct_of_tier_on_tree := round(100 * N / catalogue, 1)]
print(comp); fwrite(comp, file.path(OUT_DIR, "table_s56_tree_tier_composition.csv"))
if (length(setdiff(tips, sh$representative))) warning(length(setdiff(tips, sh$representative)), " tree tips are not catalogue representatives.")

hr("2  Families")
st <- fread(STATS_F)
fam <- st[direction %in% c("core_broad_enriched", "host_specific_enriched"),
          .(cazy_family, tier_direction = direction, OR_unadjusted = OR_core_vs_hostspecific)]
if (file.exists(GLMM_F)) {
  g <- fread(GLMM_F)[adjustment == "adjusted for log genome size + MAG completeness",
                     .(cazy_family = trait, OR_adjusted = OR_glmm, direction_adjusted = direction_glmm)]
  fam <- merge(fam, g, by = "cazy_family", all.x = TRUE)
}
cat(sprintf("Unadjusted tier-associated families: %d (core/broad %d, host-specific %d)\n", nrow(fam),
            fam[tier_direction == "core_broad_enriched", .N], fam[tier_direction == "host_specific_enriched", .N]))
extra <- setdiff(KEY, fam$cazy_family)
pres <- fread(PRESENCE_F)
pres <- unique(pres[genome %chin% tips, .(genome, cazy_family)])
plist <- split(pres$genome, pres$cazy_family)

# ---------------- castor fit (verbatim from 42_) ----------------
run_asr <- function(trait_name, present_genomes) {
  tip_states <- ifelse(tips %in% present_genomes, 2L, 1L)
  n_present <- sum(tip_states == 2L); n_absent <- sum(tip_states == 1L)
  base <- data.table(trait = trait_name, n_present = n_present, n_absent = n_absent)
  if (n_present < 5 || n_absent < 5) return(base[, `:=`(gain_rate = NA_real_, loss_rate = NA_real_, loglik = NA_real_, status = "skipped_low_n")])
  set.seed(42)
  fit <- tryCatch(asr_mk_model(tr, tip_states, Nstates = 2, rate_model = "ARD",
                               include_ancestral_likelihoods = FALSE, Ntrials = 10), error = function(e) e)
  escalated <- FALSE
  if (!inherits(fit, "error") && max(fit$transition_matrix[1, 2], fit$transition_matrix[2, 1]) > RATE_CAP) {
    escalated <- TRUE
    fit2 <- tryCatch(asr_mk_model(tr, tip_states, Nstates = 2, rate_model = "ARD",
                                  include_ancestral_likelihoods = FALSE, Ntrials = 50), error = function(e) e)
    if (!inherits(fit2, "error") && fit2$loglikelihood > fit$loglikelihood) fit <- fit2
  }
  if (inherits(fit, "error")) return(base[, `:=`(gain_rate = NA_real_, loss_rate = NA_real_, loglik = NA_real_, status = "fit_failed")])
  gain <- fit$transition_matrix[1, 2]; loss <- fit$transition_matrix[2, 1]
  base[, `:=`(gain_rate = gain, loss_rate = loss, loglik = fit$loglikelihood, escalated = escalated,
              status = if (max(gain, loss) > RATE_CAP) "high_rate_check" else "ok")]
}

hr("3  castor fits")
todo <- c(fam$cazy_family, extra)
cat(sprintf("Fitting %d families on %d cores ...\n", length(todo), N_CORES))
res <- rbindlist(mclapply(todo, function(f) run_asr(f, plist[[f]] %||% character(0)), mc.cores = N_CORES), fill = TRUE)
res[, loss_to_gain := fifelse(gain_rate > 0, loss_rate / gain_rate, NA_real_)]
res <- merge(res, fam, by.x = "trait", by.y = "cazy_family", all.x = TRUE)
res[is.na(tier_direction), tier_direction := "not_tier_associated"]
fwrite(res[order(tier_direction, -loss_to_gain)], file.path(OUT_DIR, "table_s56_castor_cazy_gain_loss.csv"))
print(res[, .N, by = .(tier_direction, status)])

# ---------------- Comparisons ----------------
hr("4  Loss:gain, core/broad- vs host-specific-enriched (unadjusted set)")
ok <- res[status == "ok" & is.finite(loss_to_gain) & loss_to_gain > 0 & tier_direction != "not_tier_associated"]
print(ok[, .(n = .N, median_loss_gain = median(loss_to_gain), median_gain = median(gain_rate),
             median_loss = median(loss_rate), n_gain_biased = sum(loss_to_gain < 1)), by = tier_direction])
cat(sprintf("Reliable estimates: %d of %d core/broad-enriched, %d of %d host-specific-enriched\n",
            ok[tier_direction == "core_broad_enriched", .N], fam[tier_direction == "core_broad_enriched", .N],
            ok[tier_direction == "host_specific_enriched", .N], fam[tier_direction == "host_specific_enriched", .N]))
if (ok[, uniqueN(tier_direction)] == 2) {
  wt <- wilcox.test(log(loss_to_gain) ~ tier_direction, data = ok)
  cat(sprintf("Wilcoxon rank-sum on log(loss:gain): W = %.1f, p = %.3g\n", wt$statistic, wt$p.value))
}
excl <- res[tier_direction != "not_tier_associated" & status != "ok", .(trait, tier_direction, status, gain_rate, loss_rate)]
cat(sprintf("\nExcluded (implausible after escalation / failed / low n): %d\n", nrow(excl))); print(excl)
cat("\nGain-biased tier-associated families (loss:gain < 1):\n")
print(ok[loss_to_gain < 1, .(trait, tier_direction, gain_rate, loss_rate, loss_to_gain = round(loss_to_gain, 2), n_present)])

if ("direction_adjusted" %in% names(res)) {
  hr("5  Families tier-associated WITHIN lineages (49 full model)")
  print(res[direction_adjusted %in% c("Core-enriched", "Host-specific-enriched"),
            .(trait, OR_unadjusted = round(OR_unadjusted, 2), OR_adjusted = round(OR_adjusted, 2), direction_adjusted,
              gain_rate, loss_rate, loss_to_gain = round(loss_to_gain, 2), n_present, status)])
}

hr("6  Key families (Table 4 CAZy rows)")
k <- res[trait %in% KEY, .(trait, tier_direction, gain_rate = signif(gain_rate, 4), loss_rate = signif(loss_rate, 4),
                           loss_to_gain = signif(loss_to_gain, 4), n_present, n_absent, status)]
if (file.exists(OLD_CASTOR_F)) {
  old <- fread(OLD_CASTOR_F)[category == "CAZy_family", .(trait, old_gain = gain_rate, old_loss = loss_rate,
                                                        old_loss_gain = loss_to_gain_ratio, old_n_present = n_present)]
  k <- merge(k, old, by = "trait", all.x = TRUE)
}
print(k)
sink()
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
})

if (run_part("cooperative")) local({
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
})
