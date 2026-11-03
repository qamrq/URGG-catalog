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
#  SCRIPT  : 49_trait_presence_taxonomy_glmm.R  -  Trait presence by tier within lineages (logistic GLMM)
#  AUTHOR  : QRQ
#  INPUT   : trait presence table (argument)
#            annotation/results_phase3/s47_cooperative_taxonomy_glmm/table_s47_per_genome_counts.csv
#  OUTPUT  : annotation/results_phase3/s49_trait_presence_taxonomy_glmm/table_s49_<label>_presence_models.csv
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tidyr); library(glmmTMB); library(parallel)
})

P <- function(var, default) Sys.getenv(var, unset = default)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 49_trait_presence_taxonomy_glmm.R <presence_table> <run_label>")
PRESENCE_F <- args[1]; RUN_LABEL <- args[2]

PER_GENOME_F <- P("PER_GENOME_F", paste0(URGG_PROJECT, "/annotation/results_phase3/s47_cooperative_taxonomy_glmm/table_s47_per_genome_counts.csv"))
OUT_DIR      <- P("OUT_DIR", paste0(URGG_PROJECT, "/annotation/results_phase3/s49_trait_presence_taxonomy_glmm"))
GENOME_COL   <- P("GENOME_COL", "genome")
TRAIT_COL    <- P("TRAIT_COL", "KO")
N_CORES      <- as.integer(P("N_CORES", "6"))
COMPLETENESS_F <- P("COMPLETENESS_F", paste0(URGG_PROJECT, "/phase_3_analysis/results_phase3/s36_methanogenesis/36h_genome_subpathway_completeness.tsv"))
COMPL_ID_COL   <- P("COMPL_ID_COL", "join_id")
COMPL_COL      <- P("COMPL_COL", "completeness")
MIN_PRESENT_FOR_MODEL <- 30
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(PRESENCE_F, PER_GENOME_F)) if (!file.exists(f)) stop("Missing: ", f)

TRAIT_LABELS <- c(K03147 = "thiC", K00788 = "thiE", K00949 = "thiN_TPK", K00946 = "thiL")

DERIVED <- list(`thiN_TPK+/thiC-` = function(w) if (all(c("thiN_TPK", "thiC") %in% names(w))) w$thiN_TPK & !w$thiC else NULL)

hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
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

hr(sprintf("Run: %s | presence table: %s", RUN_LABEL, basename(PRESENCE_F)))
dat <- read_csv(PER_GENOME_F, show_col_types = FALSE) %>%
  mutate(tier = factor(tier, levels = c("HostSpecific", "Shared", "CoreBroad")),
         log_genes_c = as.numeric(scale(log(n_genes_total))))
for (r in c("phylum_id", "class_id", "order_id", "family_id", "genus_id")) dat[[r]] <- factor(dat[[r]])

HAVE_COMPL <- nzchar(COMPLETENESS_F) && file.exists(COMPLETENESS_F)
if (HAVE_COMPL) {
  cm <- fread(COMPLETENESS_F, select = c(COMPL_ID_COL, COMPL_COL))
  setnames(cm, c("genome", "mag_completeness"))
  cm <- unique(cm, by = "genome")
  dat <- dat %>% left_join(as_tibble(cm), by = "genome")
  n_miss <- sum(is.na(dat$mag_completeness))
  cat(sprintf("MAG completeness joined for %d of %d genomes\n", nrow(dat) - n_miss, nrow(dat)))
  if (n_miss / nrow(dat) > 0.01) {
    warning(n_miss, " genomes lack a completeness value — completeness-adjusted models are skipped.")
    HAVE_COMPL <- FALSE
  } else {
    dat <- dat %>% mutate(mag_completeness = ifelse(is.na(mag_completeness), median(mag_completeness, na.rm = TRUE), mag_completeness),
                          compl_c = as.numeric(scale(mag_completeness)))
    print(dat %>% group_by(tier) %>% summarise(median_completeness = median(mag_completeness), .groups = "drop"))
  }
} else cat("No completeness file — completeness-adjusted models skipped.\n")

pl <- fread(PRESENCE_F)
if (!all(c(GENOME_COL, TRAIT_COL) %in% names(pl))) {
  stop("Presence table needs columns '", GENOME_COL, "' and '", TRAIT_COL, "'. Found: ", paste(names(pl), collapse = ", "),
       "\nSet GENOME_COL / TRAIT_COL env vars if they are named differently.")
}
pl <- unique(pl[, .(genome = get(GENOME_COL), trait = as.character(get(TRAIT_COL)))])
pl[trait %in% names(TRAIT_LABELS), trait := TRAIT_LABELS[trait]]
n_unknown <- pl[!genome %in% dat$genome, uniqueN(genome)]
if (n_unknown > 0) warning(n_unknown, " genomes in the presence table are not representatives in 47_'s table (ignored). ",
                           "If this is large, check genome ID parsing.")
pl <- pl[genome %in% dat$genome]

wide <- dcast(pl[, v := TRUE], genome ~ trait, value.var = "v", fill = FALSE)
wide <- merge(data.table(genome = dat$genome), wide, by = "genome", all.x = TRUE)
traits <- setdiff(names(wide), "genome")
for (t in traits) set(wide, which(is.na(wide[[t]])), t, FALSE)
for (nm in names(DERIVED)) { v <- DERIVED[[nm]](wide); if (!is.null(v)) { wide[, (nm) := v]; traits <- c(traits, nm) } }
cat(sprintf("Traits: %s\nGenomes (denominator): %d\n", paste(traits, collapse = ", "), nrow(wide)))

d_all <- dat %>% inner_join(as_tibble(wide), by = "genome")

ADJ_LEVELS <- c("unadjusted", "adjusted for log genome size",
                if (HAVE_COMPL) "adjusted for log genome size + MAG completeness")
ADJ_KEEP <- P("ADJ_KEEP", "")
if (nzchar(ADJ_KEEP)) {
  key <- c(unadjusted = "unadjusted", size = "adjusted for log genome size",
           compl = "adjusted for MAG completeness",
           full = "adjusted for log genome size + MAG completeness")
  want <- trimws(strsplit(ADJ_KEEP, ",")[[1]])
  if (length(setdiff(want, names(key)))) stop("ADJ_KEEP accepts only: unadjusted, size, compl, full")
  if (!HAVE_COMPL && any(c("compl", "full") %in% want)) warning("No completeness data — 'compl'/'full' skipped.")
  ADJ_LEVELS <- unname(key[want][HAVE_COMPL | !want %in% c("compl", "full")])
}
cat("Adjustment levels:", paste(ADJ_LEVELS, collapse = " | "), "\n")
jobs <- expand_grid(trait = traits, adjust = ADJ_LEVELS)
FX  <- c("unadjusted" = "tier", "adjusted for log genome size" = "tier + log_genes_c",
         "adjusted for MAG completeness" = "tier + compl_c",
         "adjusted for log genome size + MAG completeness" = "tier + log_genes_c + compl_c")
FX0 <- c("unadjusted" = "1", "adjusted for log genome size" = "log_genes_c",
         "adjusted for MAG completeness" = "compl_c",
         "adjusted for log genome size + MAG completeness" = "log_genes_c + compl_c")

run_one <- function(i) {
  tt <- jobs$trait[i]; adj <- jobs$adjust[i]
  d <- d_all %>% mutate(y = as.integer(.data[[tt]]))
  prev <- tapply(d$y, d$tier, mean) * 100
  base <- tibble(trait = tt, adjustment = adj, n_present = sum(d$y),
                 pct_core = prev[["CoreBroad"]], pct_shared = prev[["Shared"]], pct_host = prev[["HostSpecific"]])
  # Fisher core vs host (as currently reported in the manuscript)
  ch <- d %>% filter(tier %in% c("CoreBroad", "HostSpecific"))
  ft <- fisher.test(table(factor(ch$tier, levels = c("CoreBroad", "HostSpecific")), factor(ch$y, levels = c(1, 0))))
  base <- base %>% mutate(OR_fisher = unname(ft$estimate), OR_fisher_low = ft$conf.int[1],
                          OR_fisher_high = ft$conf.int[2], p_fisher = ft$p.value)
  if (min(sum(d$y), sum(1 - d$y)) < MIN_PRESENT_FOR_MODEL) return(base %>% mutate(status = "skipped_low_n"))
  fx  <- FX[[adj]]; fx0 <- FX0[[adj]]
  fp <- fit_pair("y", fx, fx0, d, binomial(), cores = 1)
  full <- fp$full; null <- fp$null
  if (is.null(full$model)) return(base %>% mutate(status = paste("fit_failed:", full$msg)))
  tc <- tier_contrasts(full$model) %>% filter(contrast == "CoreBroad vs HostSpecific")
  om <- if (!is.null(null$model)) lrt(full$model, null$model) else tibble(LRT_chisq = NA, LRT_df = NA, LRT_p = NA)
  bind_cols(base, tibble(OR_glmm = tc$OR, OR_glmm_low = tc$OR_ci_low, OR_glmm_high = tc$OR_ci_high, p_wald = tc$p_wald),
            om, tibble(re_structure = full$re_structure,
                       status = ifelse(full$converged && null$converged, "ok",
                                       paste0("check_convergence: full=", full$converged, " null=", null$converged))))
}

res <- bind_rows(mclapply(seq_len(nrow(jobs)), run_one, mc.cores = min(N_CORES, nrow(jobs))))
res <- res %>% group_by(adjustment) %>%
  mutate(p_fisher_BH = p.adjust(p_fisher, "BH"), LRT_p_BH = p.adjust(LRT_p, "BH"),
         direction_fisher = case_when(p_fisher_BH < 0.05 & OR_fisher_low > 1 ~ "Core-enriched",
                                      p_fisher_BH < 0.05 & OR_fisher_high < 1 ~ "Host-specific-enriched", TRUE ~ "ns"),
         direction_glmm = case_when(is.na(OR_glmm) ~ NA_character_,
                                    LRT_p_BH < 0.05 & OR_glmm_low > 1 ~ "Core-enriched",
                                    LRT_p_BH < 0.05 & OR_glmm_high < 1 ~ "Host-specific-enriched", TRUE ~ "ns")) %>%
  ungroup() %>% arrange(adjustment, trait)

out_f <- file.path(OUT_DIR, sprintf("table_s49_%s_presence_models.csv", RUN_LABEL))
write_csv(res, out_f)
options(width = 200)
show <- res %>% select(adjustment, trait, pct_core, pct_shared, pct_host, OR_fisher, p_fisher_BH,
                       OR_glmm, OR_glmm_low, OR_glmm_high, LRT_p_BH, direction_fisher, direction_glmm, status)
if (length(traits) > 40) {
  cat("\nDirection counts (Fisher = unadjusted descriptive; GLMM = taxonomy-controlled):\n")
  print(as.data.frame(res %>% count(adjustment, direction_fisher, direction_glmm)))
  cat("\nFit status:\n"); print(as.data.frame(res %>% count(adjustment, ok = status == "ok")))
  show <- show %>% filter(direction_glmm %in% c("Core-enriched", "Host-specific-enriched")) %>% arrange(adjustment, LRT_p_BH)
  cat("\nTraits significant after taxonomy control:\n")
}
print(as.data.frame(show), digits = 3)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, sprintf("sessionInfo_s49_%s.txt", RUN_LABEL)))
cat("\nWritten: ", out_f, "\n", sep = "")
