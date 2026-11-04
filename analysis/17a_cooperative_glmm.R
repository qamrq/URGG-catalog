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
#  SCRIPT  : 17a_cooperative_glmm.R  -  Cooperative genes by tier within lineages (GLMM)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s36_cooperative_ko_core_specific/
#            annotation/eggnog_output/nr_catalog.emapper.annotations
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s47_cooperative_taxonomy_glmm/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tidyr)
  library(stringr); library(glmmTMB); library(parallel)
})

P <- function(var, default) Sys.getenv(var, unset = default)
S36 <- paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific")

CHECKPOINT_RDS  <- P("CHECKPOINT_RDS",  file.path(S36, "gene_level_cooperative_join_CHECKPOINT.rds"))
PREV_PROP_F     <- P("PREV_PROP_F",     file.path(S36, "table_cooperative_proportion_per_genome.csv"))
PRESENCE_WIDE_F <- P("PRESENCE_WIDE_F", file.path(S36, "table_cooperative_category_presence_per_genome_WIDE.csv"))
EGGNOG_F        <- P("EGGNOG_F",  paste0(URGG_PROJECT, "/annotation/eggnog_output/nr_catalog.emapper.annotations"))
SHARING_F       <- P("SHARING_F", paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv"))
GTDB_BAC_F      <- P("GTDB_BAC_F", paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232/gtdbtk.bac120.summary.tsv"))
GTDB_AR_F       <- P("GTDB_AR_F",  paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232/gtdbtk.ar53.summary.tsv"))
OUT_DIR         <- P("OUT_DIR", paste0(URGG_PROJECT, "/annotation/results_phase3/s47_cooperative_taxonomy_glmm"))
N_CORES         <- as.integer(P("N_CORES", "6"))
EXPECTED_KO_REPS <- as.numeric(P("EXPECTED_KO_REPS", "17702338"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(CHECKPOINT_RDS, PREV_PROP_F, PRESENCE_WIDE_F, EGGNOG_F, SHARING_F, GTDB_BAC_F, GTDB_AR_F)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

COOP_CATEGORIES <- c("siderophores", "antibiotic_degradation", "biofilm",
                     "quorum_sensing", "secretion_system", "Extracellular")
TIER_MAP <- c("Host-specific (1 host species)" = "HostSpecific",
              "Shared (2-3 host species)"      = "Shared",
              "Core / broad (>3 host species)" = "CoreBroad")
MIN_PRESENT_FOR_MODEL <- 30

hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")

COUNTS_CACHE <- file.path(OUT_DIR, "cache_s47_per_genome_counts.rds")
if (file.exists(COUNTS_CACHE)) {
  hr("STEPS 1-2  Loading cached per-genome counts")
  per_genome <- readRDS(COUNTS_CACHE)
  cat(sprintf("Loaded %d genomes from %s\n", nrow(per_genome), COUNTS_CACHE))
} else {
hr("STEP 1  KO-annotated representative genes")
ko_tab <- fread(cmd = sprintf("grep -v '^##' %s | cut -f1,12", shQuote(EGGNOG_F)),
                sep = "\t", header = TRUE, quote = "", colClasses = "character")
setnames(ko_tab, c("rep_gene_id", "KEGG_ko"))
ko_reps <- unique(ko_tab[!is.na(KEGG_ko) & KEGG_ko != "-" & KEGG_ko != "", rep_gene_id])
rm(ko_tab); invisible(gc())
cat(sprintf("Representative genes with >=1 KO: %s (expected %s)\n",
            format(length(ko_reps), big.mark = ","), format(EXPECTED_KO_REPS, big.mark = ",")))
if (!is.na(EXPECTED_KO_REPS) && length(ko_reps) != EXPECTED_KO_REPS) {
  warning("KO-annotated representative count differs from the expected value (",
          EXPECTED_KO_REPS, "). Check that column 12 of this eggNOG file is KEGG_ko before continuing.")
}

hr("STEP 2  Per-genome gene counts (all genes / KO genes / cooperative genes)")
cl <- readRDS(CHECKPOINT_RDS); setDT(cl)
need <- c("rep_gene_id", "is_cooperative", "genome")
if (!all(need %in% names(cl))) stop("Checkpoint lacks columns: ", paste(setdiff(need, names(cl)), collapse = ", "))

cl[, has_ko := rep_gene_id %chin% ko_reps]
rm(ko_reps); invisible(gc())

n_coop_no_ko <- cl[is_cooperative == TRUE & has_ko == FALSE, .N]
if (n_coop_no_ko > 0) stop(n_coop_no_ko, " cooperative genes have no KO - the checkpoint and eggNOG file are out of sync.")
if ("member_gene_id" %in% names(cl) && anyDuplicated(cl$member_gene_id)) {
  stop("Duplicate member_gene_id rows in checkpoint - genes would be double-counted.")
}

per_genome <- cl[, .(n_genes_total = .N,
                     n_ko_genes    = sum(has_ko),
                     n_coop        = sum(is_cooperative)), by = genome]
rm(cl); invisible(gc())
saveRDS(per_genome, COUNTS_CACHE)
}

prev <- read_csv(PREV_PROP_F, show_col_types = FALSE) %>%
  select(genome, n_cooperative_genes, n_genes_prodigal, sharing_category_prev = sharing_category)
chk <- inner_join(as_tibble(per_genome), prev, by = "genome")
n_bad <- sum(chk$n_genes_total != chk$n_genes_prodigal | chk$n_coop != chk$n_cooperative_genes)
cat(sprintf("Cross-check vs s36 table: %d genomes compared, %d mismatches\n", nrow(chk), n_bad))
if (n_bad > 0) stop("Per-genome counts do not reproduce the s36 table.")

hr("STEP 3  GTDB taxonomy")
read_gtdb <- function(f) fread(f, sep = "\t", quote = "", select = c("user_genome", "classification"),
                               colClasses = "character")
tax <- rbind(read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F))
if (anyDuplicated(tax$user_genome)) stop("Genome present in both bac120 and ar53 summaries - resolve before modelling.")

get_rank <- function(x, prefix) {
  v <- str_match(x, paste0("(?:^|;)", prefix, "__([^;]*)"))[, 2]
  ifelse(is.na(v) | v == "", NA_character_, v)
}
RANKS <- c(domain = "d", phylum = "p", class = "c", order = "o", family = "f", genus = "g")
for (r in names(RANKS)) tax[, (r) := get_rank(classification, RANKS[[r]])]

make_ids <- function(dt) {
  parent <- ifelse(is.na(dt$domain), paste0("unassigned:", dt$user_genome), dt$domain)
  for (r in c("phylum", "class", "order", "family", "genus")) {
    val <- dt[[r]]
    id  <- ifelse(is.na(val), paste0(parent, "|", r, "_unassigned:", dt$user_genome), paste0(parent, "|", val))
    dt[, (paste0(r, "_id")) := id]
    parent <- id
  }
  dt
}
tax <- make_ids(tax)
unassigned <- sapply(c("phylum", "class", "order", "family", "genus"), function(r) sum(is.na(tax[[r]])))
cat("Genomes with an unassigned rank (treated as singleton levels):\n"); print(unassigned)

hr("STEP 4  Analysis table")
sharing <- read_csv(SHARING_F, show_col_types = FALSE) %>%
  select(genome = representative, sharing_category, n_host_species)
unmapped <- setdiff(unique(sharing$sharing_category), names(TIER_MAP))
if (length(unmapped)) stop("Unrecognised sharing_category label(s): ", paste(unmapped, collapse = " | "))

dat <- sharing %>%
  inner_join(as_tibble(per_genome), by = "genome") %>%
  left_join(as_tibble(tax) %>% select(genome = user_genome, domain, ends_with("_id")), by = "genome") %>%
  mutate(tier = factor(unname(TIER_MAP[sharing_category]), levels = c("HostSpecific", "Shared", "CoreBroad")))

n_no_tax <- sum(is.na(dat$phylum_id))
cat(sprintf("Genomes in analysis: %d of %d representatives | without a GTDB-Tk row: %d\n",
            nrow(dat), nrow(sharing), n_no_tax))
if (nrow(dat) < nrow(sharing)) warning(nrow(sharing) - nrow(dat), " representatives lack per-genome counts.")
if (n_no_tax > 0) {
  if (n_no_tax / nrow(dat) > 0.05) stop(">5% of genomes have no GTDB-Tk row - check genome ID format (e.g. file extensions).")
  warning(n_no_tax, " genomes have no GTDB-Tk row; each is given singleton taxonomy levels.")
  for (r in c("phylum", "class", "order", "family", "genus")) {
    idc <- paste0(r, "_id")
    dat[[idc]] <- ifelse(is.na(dat[[idc]]), paste0("noGTDB:", dat$genome, "|", r), dat[[idc]])
  }
}
if (any(dat$n_ko_genes == 0)) {
  warning(sum(dat$n_ko_genes == 0), " genomes have zero KO-annotated genes and are dropped from the KO-denominator models.")
}
for (r in c("phylum_id", "class_id", "order_id", "family_id", "genus_id")) dat[[r]] <- factor(dat[[r]])
dat <- dat %>% mutate(log_genes_c = as.numeric(scale(log(n_genes_total))))

write_csv(dat %>% select(genome, sharing_category, tier, n_host_species, n_genes_total, n_ko_genes, n_coop,
                         domain, ends_with("_id")),
          file.path(OUT_DIR, "table_s47_per_genome_counts.csv"))
cat("Genomes per tier:\n"); print(table(dat$tier))

hr("STEP 5  KO-annotation rate by tier (the denominator confound)")
ann <- dat %>% mutate(ko_rate = n_ko_genes / n_genes_total) %>%
  group_by(tier) %>%
  summarise(n_genomes = n(), median_ko_rate = median(ko_rate), iqr_low = quantile(ko_rate, .25),
            iqr_high = quantile(ko_rate, .75), median_prop_coop_allgenes = median(n_coop / n_genes_total),
            median_prop_coop_kogenes = median(n_coop / pmax(n_ko_genes, 1)), .groups = "drop")
kw <- kruskal.test(n_ko_genes / n_genes_total ~ tier, data = dat)
ann$kruskal_p_ko_rate <- kw$p.value
print(as.data.frame(ann))
write_csv(ann, file.path(OUT_DIR, "table_s47_annotation_rate_by_tier.csv"))

RE_FULL  <- "(1 | phylum_id) + (1 | class_id) + (1 | order_id) + (1 | family_id) + (1 | genus_id)"
RE_NOGEN <- "(1 | phylum_id) + (1 | class_id) + (1 | order_id) + (1 | family_id)"

fit_tmb <- function(formula, data, family, cores = N_CORES) {
  warns <- character(0)
  try_fit <- function(ctrl) {
    withCallingHandlers(
      tryCatch(glmmTMB(formula, data = data, family = family, control = ctrl), error = function(e) e),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  }
  m <- try_fit(glmmTMBControl(parallel = cores))
  ok <- function(m) {
    if (inherits(m, "error") || !isTRUE(m$sdr$pdHess) || m$fit$convergence != 0) return(FALSE)
    se <- suppressWarnings(sqrt(diag(vcov(m)$cond)))
    all(is.finite(se))
  }
  if (!ok(m)) {
    m2 <- try_fit(glmmTMBControl(parallel = cores, optimizer = optim, optArgs = list(method = "BFGS")))
    if (ok(m2)) m <- m2
  }
  list(model = if (inherits(m, "error")) NULL else m, converged = ok(m),
       msg = if (inherits(m, "error")) conditionMessage(m) else paste(unique(warns), collapse = " | "))
}

fit_pair <- function(lhs, fx, fx0, data, family, cores = N_CORES) {
  for (re in list(list(f = RE_FULL,  lab = "phylum/class/order/family/genus"),
                  list(f = RE_NOGEN, lab = "phylum/class/order/family (genus dropped)"))) {
    full <- fit_tmb(as.formula(paste(lhs, "~", fx,  "+", re$f)), data, family, cores)
    null <- fit_tmb(as.formula(paste(lhs, "~", fx0, "+", re$f)), data, family, cores)
    if (full$converged && null$converged) break
  }
  full$re_structure <- re$lab
  list(full = full, null = null)
}

tier_contrasts <- function(m) {
  cf <- summary(m)$coefficients$cond
  rows <- grep("^tier", rownames(cf))
  tibble(term = rownames(cf)[rows], beta = cf[rows, "Estimate"], se = cf[rows, "Std. Error"],
         p_wald = cf[rows, "Pr(>|z|)"]) %>%
    mutate(contrast = recode(term, tierShared = "Shared vs HostSpecific", tierCoreBroad = "CoreBroad vs HostSpecific"),
           OR = exp(beta), OR_ci_low = exp(beta - 1.959964 * se), OR_ci_high = exp(beta + 1.959964 * se))
}

lrt <- function(m_full, m_null) {
  a <- anova(m_null, m_full)
  tibble(LRT_chisq = a$Chisq[2], LRT_df = a$`Chi Df`[2], LRT_p = a$`Pr(>Chisq)`[2])
}

hr("6a  Quasi-binomial GLM (s36 specification)")
q_full <- glm(cbind(n_coop, n_genes_total - n_coop) ~ tier, data = dat, family = quasibinomial())
pair   <- droplevels(filter(dat, tier %in% c("HostSpecific", "CoreBroad")))
q_pair <- glm(cbind(n_coop, n_genes_total - n_coop) ~ tier, data = pair, family = quasibinomial())
q_om   <- anova(q_full, test = "F")
b      <- coef(summary(q_pair))["tierCoreBroad", ]
reproduced <- tibble(
  model = "M0 quasi-binomial GLM (as published, two-tier refit)", denominator = "all genes", taxonomy_RE = "no",
  re_structure = NA_character_, contrast = "CoreBroad vs HostSpecific",
  OR = exp(b[["Estimate"]]), OR_ci_low = exp(b[["Estimate"]] - qt(.975, q_pair$df.residual) * b[["Std. Error"]]),
  OR_ci_high = exp(b[["Estimate"]] + qt(.975, q_pair$df.residual) * b[["Std. Error"]]),
  p_wald = b[["Pr(>|t|)"]], omnibus_test = "F", omnibus_stat = q_om$F[2], omnibus_df = 2, omnibus_p = q_om$`Pr(>F)`[2],
  n_genomes = nrow(dat), converged = TRUE, note = "t-based CI")
print(as.data.frame(reproduced %>% select(model, OR, OR_ci_low, OR_ci_high, omnibus_p)))

hr("STEP 6b  Beta-binomial models: 2 x 2 decomposition (denominator x taxonomy)")
dat_ko <- filter(dat, n_ko_genes > 0)
specs <- list(
  list(id = "M1", label = "Beta-binomial GLM",  denom = "all genes", tax = FALSE, data = dat,
       lhs = "cbind(n_coop, n_genes_total - n_coop)"),
  list(id = "M2", label = "Beta-binomial GLM",  denom = "KO genes",  tax = FALSE, data = dat_ko,
       lhs = "cbind(n_coop, n_ko_genes - n_coop)"),
  list(id = "M3", label = "Beta-binomial GLMM", denom = "all genes", tax = TRUE,  data = dat,
       lhs = "cbind(n_coop, n_genes_total - n_coop)"),
  list(id = "M4", label = "Beta-binomial GLMM (PRIMARY)", denom = "KO genes", tax = TRUE, data = dat_ko,
       lhs = "cbind(n_coop, n_ko_genes - n_coop)")
)

models <- list(); pooled_rows <- list(); re_rows <- list()
for (s in specs) {
  message(sprintf("Fitting %s (%s, %s, taxonomy RE = %s) ...", s$id, s$label, s$denom, s$tax))
  t0 <- Sys.time()
  if (s$tax) {
    fp <- fit_pair(s$lhs, "tier", "1", s$data, betabinomial(link = "logit"))
    full <- fp$full; null <- fp$null
  } else {
    full <- fit_tmb(as.formula(paste(s$lhs, "~ tier")), s$data, betabinomial(link = "logit"))
    full$re_structure <- "none"
    null <- fit_tmb(as.formula(paste(s$lhs, "~ 1")), s$data, betabinomial(link = "logit"))
  }
  cat(sprintf("  %s done in %.1f min | converged: full=%s null=%s | RE: %s\n", s$id,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), full$converged, null$converged, full$re_structure))
  if (nzchar(full$msg)) cat("  messages:", full$msg, "\n")
  models[[s$id]] <- list(full = full$model, null = null$model, spec = s[c("id", "label", "denom", "tax")])
  if (is.null(full$model)) { warning(s$id, " failed: ", full$msg); next }

  om <- if (!is.null(null$model)) lrt(full$model, null$model) else tibble(LRT_chisq = NA, LRT_df = NA, LRT_p = NA)
  pooled_rows[[s$id]] <- tier_contrasts(full$model) %>%
    transmute(model = paste(s$id, s$label), denominator = s$denom, taxonomy_RE = ifelse(s$tax, "yes", "no"),
              re_structure = full$re_structure, contrast, OR, OR_ci_low, OR_ci_high, p_wald,
              omnibus_test = "LRT", omnibus_stat = om$LRT_chisq, omnibus_df = om$LRT_df, omnibus_p = om$LRT_p,
              n_genomes = nrow(s$data), converged = full$converged && null$converged,
              note = ifelse(full$converged, "Wald 95% CI", paste("NOT CONVERGED:", full$msg)))
  if (s$tax) {
    vc <- VarCorr(full$model)$cond
    re_rows[[s$id]] <- tibble(model = s$id, level = names(vc),
                              variance = sapply(vc, function(x) attr(x, "stddev")^2),
                              n_levels = sapply(names(vc), function(g) nlevels(s$data[[g]])))
  }
}

pooled <- bind_rows(reproduced, bind_rows(pooled_rows))
write_csv(pooled, file.path(OUT_DIR, "table_s47_pooled_models.csv"))
write_csv(bind_rows(re_rows), file.path(OUT_DIR, "table_s47_pooled_random_effects.csv"))

cat("\n--- Core/broad vs host-specific, pooled cooperative proportion ---\n")
print(as.data.frame(pooled %>% filter(contrast == "CoreBroad vs HostSpecific") %>%
                      select(model, denominator, taxonomy_RE, OR, OR_ci_low, OR_ci_high, p_wald, omnibus_p, converged)),
      digits = 4)
cat("\nTaxonomic random-effect variances (latent logit scale):\n"); print(as.data.frame(bind_rows(re_rows)), digits = 4)

hr("STEP 7  Six cooperative categories: genome-level presence (logistic GLMM)")
pres <- read_csv(PRESENCE_WIDE_F, show_col_types = FALSE)
miss_cat <- setdiff(COOP_CATEGORIES, names(pres))
if (length(miss_cat)) stop("Presence table lacks categories: ", paste(miss_cat, collapse = ", "))
dat_p <- dat %>% inner_join(pres %>% select(genome, all_of(COOP_CATEGORIES)), by = "genome")
if (nrow(dat_p) != nrow(dat)) warning(nrow(dat) - nrow(dat_p), " genomes missing from the presence table.")

cat_jobs <- expand_grid(category = COOP_CATEGORIES,
                        adjust = c("unadjusted", "adjusted for log genome size"))

run_cat <- function(i) {
  cc <- cat_jobs$category[i]; adj <- cat_jobs$adjust[i]
  d <- dat_p %>% mutate(y = as.integer(as.logical(.data[[cc]])))
  n_pres <- sum(d$y); n_abs <- sum(1 - d$y)
  base <- tibble(category = cc, adjustment = adj, n_present = n_pres, n_absent = n_abs,
                 pct_core = 100 * mean(d$y[d$tier == "CoreBroad"]),
                 pct_shared = 100 * mean(d$y[d$tier == "Shared"]),
                 pct_host = 100 * mean(d$y[d$tier == "HostSpecific"]))
  if (min(n_pres, n_abs) < MIN_PRESENT_FOR_MODEL) return(base %>% mutate(status = "skipped_low_n"))
  fx <- if (adj == "unadjusted") "tier" else "tier + log_genes_c"
  fx0 <- if (adj == "unadjusted") "1" else "log_genes_c"

  g1 <- glm(as.formula(paste("y ~", fx)), data = d, family = binomial())
  g0 <- glm(as.formula(paste("y ~", fx0)), data = d, family = binomial())
  gb <- coef(summary(g1))["tierCoreBroad", ]
  naive <- tibble(OR_naive = exp(gb[["Estimate"]]),
                  OR_naive_low = exp(gb[["Estimate"]] - 1.959964 * gb[["Std. Error"]]),
                  OR_naive_high = exp(gb[["Estimate"]] + 1.959964 * gb[["Std. Error"]]),
                  LRT_p_naive = anova(g0, g1, test = "LRT")$`Pr(>Chi)`[2])

  fp <- fit_pair("y", fx, fx0, d, binomial(), cores = 1)
  full <- fp$full; null <- fp$null
  if (is.null(full$model)) return(bind_cols(base, naive) %>% mutate(status = paste("fit_failed:", full$msg)))
  tc <- tier_contrasts(full$model) %>% filter(contrast == "CoreBroad vs HostSpecific")
  om <- if (!is.null(null$model)) lrt(full$model, null$model) else tibble(LRT_chisq = NA, LRT_df = NA, LRT_p = NA)
  bind_cols(base, naive,
            tibble(OR = tc$OR, OR_ci_low = tc$OR_ci_low, OR_ci_high = tc$OR_ci_high, p_wald = tc$p_wald),
            om,
            tibble(re_structure = full$re_structure,
                   status = ifelse(full$converged && null$converged, "ok",
                                   paste0("check_convergence: full=", full$converged, " null=", null$converged,
                                          " ", full$msg, " ", null$msg))))
}

cat_res <- bind_rows(mclapply(seq_len(nrow(cat_jobs)), run_cat, mc.cores = min(N_CORES, nrow(cat_jobs))))
cat_res <- cat_res %>%
  group_by(adjustment) %>%
  mutate(LRT_p_BH = p.adjust(LRT_p, method = "BH"),
         p_wald_BH = p.adjust(p_wald, method = "BH"),
         LRT_p_naive_BH = p.adjust(LRT_p_naive, method = "BH"),
         direction = case_when(
           is.na(OR) ~ NA_character_,
           LRT_p_BH < 0.05 & OR_ci_low > 1 ~ "Core-enriched",
           LRT_p_BH < 0.05 & OR_ci_high < 1 ~ "Host-specific-enriched",
           TRUE ~ "ns")) %>%
  ungroup() %>% arrange(adjustment, category)

write_csv(cat_res, file.path(OUT_DIR, "table_s47_category_presence_models.csv"))
print(as.data.frame(cat_res %>% select(adjustment, category, pct_core, pct_host, OR_naive, OR, OR_ci_low, OR_ci_high,
                                       LRT_p_BH, direction, status)), digits = 3)

saveRDS(list(pooled_models = models, data = dat), file.path(OUT_DIR, "models_s47.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_s47.txt"))
hr("DONE")
cat("Outputs in: ", OUT_DIR, "\n", sep = "")
cat("Primary result = row 'M4 Beta-binomial GLMM (PRIMARY)' in table_s47_pooled_models.csv\n")
