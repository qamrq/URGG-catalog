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
#  SCRIPT  : 48_cooperative_mcmcglmm_phylo.R  -  Phylogenetic sensitivity of the cooperative result (MCMCglmm)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s41_castor_gain_loss/catalogue_16382_midpointrooted.tree
#            annotation/results_phase3/s47_cooperative_taxonomy_glmm/table_s47_per_genome_counts.csv
#  OUTPUT  : annotation/results_phase3/s48_cooperative_mcmcglmm/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({
  library(ape); library(MCMCglmm); library(coda); library(dplyr); library(readr); library(parallel)
})

P <- function(var, default) Sys.getenv(var, unset = default)
args   <- commandArgs(trailingOnly = TRUE)
TARGET <- if (length(args) >= 1) args[1] else "pooled"
MODE   <- if (length(args) >= 2) args[2] else "pilot"

TREE_F       <- P("TREE_F", paste0(URGG_PROJECT, "/annotation/results_phase3/s41_castor_gain_loss/catalogue_16382_midpointrooted.tree"))
S47_DIR      <- P("S47_DIR", paste0(URGG_PROJECT, "/annotation/results_phase3/s47_cooperative_taxonomy_glmm"))
PER_GENOME_F <- P("PER_GENOME_F", file.path(S47_DIR, "table_s47_per_genome_counts.csv"))
PRESENCE_WIDE_F <- P("PRESENCE_WIDE_F", paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific/table_cooperative_category_presence_per_genome_WIDE.csv"))
OUT_DIR      <- P("OUT_DIR", paste0(URGG_PROJECT, "/annotation/results_phase3/s48_cooperative_mcmcglmm"))
EXPECTED_TIPS <- as.integer(P("EXPECTED_TIPS", "16382"))

NITT   <- as.integer(P("NITT",   "103000"))
BURNIN <- as.integer(P("BURNIN", "3000"))
THIN   <- as.integer(P("THIN",   "100"))
N_CHAINS <- as.integer(P("N_CHAINS", "2"))

COOP_CATEGORIES <- c("siderophores", "antibiotic_degradation", "biofilm",
                     "quorum_sensing", "secretion_system", "Extracellular")
if (!TARGET %in% c("pooled", COOP_CATEGORIES)) stop("Unknown target: ", TARGET)
if (!MODE %in% c("pilot", "full")) stop("MODE must be 'pilot' or 'full'")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(TREE_F, PER_GENOME_F, if (TARGET != "pooled") PRESENCE_WIDE_F)) if (!file.exists(f)) stop("Missing: ", f)

hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

hr(sprintf("Target: %s | mode: %s", TARGET, MODE))
dat <- read_csv(PER_GENOME_F, show_col_types = FALSE) %>%
  mutate(tier = factor(tier, levels = c("HostSpecific", "Shared", "CoreBroad")))
if (TARGET != "pooled") {
  pres <- read_csv(PRESENCE_WIDE_F, show_col_types = FALSE) %>% select(genome, all_of(TARGET))
  dat <- dat %>% inner_join(pres, by = "genome") %>% mutate(y = as.integer(as.logical(.data[[TARGET]])))
}

tr <- read.tree(TREE_F)
cat(sprintf("Tree tips: %d (expected %d)\n", Ntip(tr), EXPECTED_TIPS))
if (!is.na(EXPECTED_TIPS) && Ntip(tr) != EXPECTED_TIPS) warning("Unexpected tip count.")
tr$node.label <- NULL
dat <- dat %>% filter(genome %in% tr$tip.label)
if (TARGET == "pooled") dat <- dat %>% filter(n_ko_genes > 0)
tr <- keep.tip(tr, dat$genome)
if (!is.binary(tr)) tr <- multi2di(tr, random = FALSE)
eps <- 1e-6 * mean(tr$edge.length[tr$edge.length > 0])
tr$edge.length[tr$edge.length < eps] <- eps
rtt <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
tr$edge.length <- tr$edge.length / mean(rtt)
cat(sprintf("Genomes analysed (on tree, with data): %d\n", nrow(dat))); print(table(dat$tier))

Ainv <- inverseA(tr, nodes = "ALL", scale = FALSE)$Ainv
keep_cols <- intersect(c("genome", "tier", "n_coop", "n_ko_genes", "y"), names(dat))
dat <- as.data.frame(dat[, keep_cols]); dat$animal <- factor(dat$genome)

if (TARGET == "pooled") {
  fixed  <- cbind(n_coop, n_ko_genes - n_coop) ~ tier
  family <- "multinomial2"
  prior  <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
} else {
  fixed  <- y ~ tier
  family <- "categorical"
  prior  <- list(B = list(mu = rep(0, 3), V = diag(3) * (1 + pi^2 / 3)),
                 R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))
}

run_chain <- function(seed, nitt, burnin, thin) {
  set.seed(seed)
  MCMCglmm(fixed, random = ~ animal, ginverse = list(animal = Ainv), family = family,
           data = dat, prior = prior, nitt = nitt, burnin = burnin, thin = thin,
           pr = FALSE, verbose = FALSE)
}

if (MODE == "pilot") {
  PN <- 1100
  t0 <- Sys.time(); m <- run_chain(1, PN, 100, 10); el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  per_iter <- el / PN
  cat(sprintf("\nPilot: %d iterations in %.1f s (%.3f s/iter)\n", PN, el, per_iter))
  cat(sprintf("Projected full run (nitt = %d, %d chains in parallel): ~%.1f hours\n",
              NITT, N_CHAINS, per_iter * NITT / 3600))
  cat("Pilot posterior means (NOT for reporting — chain far too short):\n"); print(colMeans(m$Sol))
  cat("\nIf the projection is acceptable, run:  Rscript 48_cooperative_mcmcglmm_phylo.R", TARGET, "full\n")
  quit(save = "no", status = 0)
}

hr(sprintf("Full run: %d chains x nitt %d (burnin %d, thin %d)", N_CHAINS, NITT, BURNIN, THIN))
t0 <- Sys.time()
chains <- mclapply(seq_len(N_CHAINS), function(s) run_chain(s, NITT, BURNIN, THIN), mc.cores = N_CHAINS)
bad <- sapply(chains, function(x) inherits(x, "try-error") || !inherits(x, "MCMCglmm"))
if (any(bad)) stop("Chain(s) failed: ", paste(sapply(chains[bad], as.character), collapse = "\n"))
cat(sprintf("Finished in %.2f h\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))
saveRDS(chains, file.path(OUT_DIR, sprintf("mcmc_s48_%s.rds", TARGET)))

sol <- mcmc.list(lapply(chains, function(x) x$Sol))
vcv <- mcmc.list(lapply(chains, function(x) x$VCV))
both <- lapply(chains, function(x) cbind(as.matrix(x$Sol), as.matrix(x$VCV)))
varying <- apply(do.call(rbind, both), 2, var) > 0
both <- mcmc.list(lapply(both, function(m) as.mcmc(m[, varying, drop = FALSE])))
psrf <- gelman.diag(both, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
ess  <- effectiveSize(both)

sol_all <- as.matrix(sol); vcv_all <- as.matrix(vcv)
pmcmc <- function(x) min(1, 2 * min(mean(x > 0), mean(x < 0)))
tier_terms <- grep("^tier", colnames(sol_all), value = TRUE)
hpd <- HPDinterval(as.mcmc(sol_all[, tier_terms, drop = FALSE]))

VG <- vcv_all[, "animal"]; VR <- vcv_all[, "units"]
phylo_signal <- VG / (VG + VR + pi^2 / 3)

summ <- tibble(
  target = TARGET, family = family,
  contrast = recode(tier_terms, tierShared = "Shared vs HostSpecific", tierCoreBroad = "CoreBroad vs HostSpecific"),
  beta_mean = colMeans(sol_all[, tier_terms, drop = FALSE]),
  OR = exp(beta_mean),
  OR_hpd_low = exp(hpd[, 1]), OR_hpd_high = exp(hpd[, 2]),
  pMCMC = apply(sol_all[, tier_terms, drop = FALSE], 2, pmcmc),
  phylo_signal_mean = mean(phylo_signal),
  phylo_signal_hpd_low = HPDinterval(as.mcmc(phylo_signal))[1],
  phylo_signal_hpd_high = HPDinterval(as.mcmc(phylo_signal))[2],
  phylo_share_excl_link = mean(VG / (VG + VR)),
  n_genomes = nrow(dat), n_chains = N_CHAINS, nitt = NITT, burnin = BURNIN, thin = THIN,
  max_psrf = max(psrf, na.rm = TRUE), min_ess = min(ess, na.rm = TRUE)
)
write_csv(summ, file.path(OUT_DIR, sprintf("table_s48_%s_summary.csv", TARGET)))
print(as.data.frame(summ), digits = 4)

if (max(psrf, na.rm = TRUE) > 1.1 || min(ess, na.rm = TRUE) < 400) {
  warning("Convergence criteria not met (PSRF > 1.1 or ESS < 400). Rerun with larger NITT/THIN, e.g. ",
          "NITT=303000 THIN=300 Rscript 48_cooperative_mcmcglmm_phylo.R ", TARGET, " full")
} else cat("\nConvergence OK (all PSRF <= 1.1, all ESS >= 400).\n")

if (TARGET == "pooled") {
  suppressPackageStartupMessages(library(glmmTMB))
  g <- glmmTMB(cbind(n_coop, n_ko_genes - n_coop) ~ tier, data = dat, family = betabinomial())
  cf <- summary(g)$coefficients$cond
} else {
  g <- glm(y ~ tier, data = dat, family = binomial()); cf <- coef(summary(g))
}
r <- "tierCoreBroad"
np <- tibble(target = TARGET, model = ifelse(TARGET == "pooled", "beta-binomial GLM, no phylogeny", "logistic GLM, no phylogeny"),
             n_genomes = nrow(dat), OR = exp(cf[r, 1]),
             OR_ci_low = exp(cf[r, 1] - 1.959964 * cf[r, 2]), OR_ci_high = exp(cf[r, 1] + 1.959964 * cf[r, 2]),
             p = cf[r, 4])
write_csv(np, file.path(OUT_DIR, sprintf("table_s48_%s_same_subset_nonphylo.csv", TARGET)))
cat("\nSame 16k subset without phylogeny (CoreBroad vs HostSpecific):\n"); print(as.data.frame(np), digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, sprintf("sessionInfo_s48_%s.txt", TARGET)))
