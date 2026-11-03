#!/usr/bin/env Rscript
###############################################################################
# 79_readbased_host_range_novelty.R
#
# QUESTION
#   Is "novel species are more often single-host" (69.4% vs 46.3% of
#   non-singleton clusters; OR 0.38 for multi-host, Fisher) partly a detection
#   artefact? Assembly needs a species to be abundant in a sample. If novel
#   species are rarer, they could look host-specific only because they were
#   never assembled from their other hosts.
#
# APPROACH
#   Re-derive host range from READS instead of assembly, using the existing
#   CoverM genome x sample matrix (35b: 1,369 samples x 25,451 representatives).
#   Nothing upstream is rerun; the abundance profiling is used exactly as
#   reported.
#
# DETECTION RULE
#   Species s is detected in host h if its relative abundance is > 0 in at
#   least MIN_SAMPLES_DETECT samples of h (primary = 2).
#   35b ran CoverM with its default --min-covered-fraction (10%), so a non-zero
#   value already means >=10% of the genome is covered by reads aligned at
#   >=95% identity over >=75% of their length, with competitive mapping against
#   all 25,451 representatives. A 30% breadth rule would need a
#   covered_fraction rerun of CoverM, which is deliberately NOT done here.
#
# STEPS
#   0. Reproduce the reported assembly-based numbers (69.4 / 46.3 / OR 0.38).
#   1. Read-based host range per species; who is undetected by reads?
#   2. Novel vs described, read-based single-host (Fisher), same subset as the
#      assembly-based estimate, plus 1/2/3-sample sensitivity and the novelty
#      depth gradient.
#   3. Logistic: single-host ~ novel + log10 abundance (for the assembly-based
#      AND the read-based outcome).
#   4. Agreement between assembly- and read-based classification
#      (supplementary table + Cohen's kappa).
#   5. Equal number of samples drawn from each host (repeated subsampling).
#
# OR CONVENTION
#   OR_single = odds of being single-host, novel vs described
#   (> 1 = novel more host-specific). The manuscript reports the inverse
#   (multi-host OR 0.38 = single-host OR 2.63); both are written out.
###############################################################################

suppressPackageStartupMessages(library(data.table))

## ============================ CONFIG ========================================
BASE <- "/disk252/qamar/metagenome_analysis"

ABUND_F     <- file.path(BASE, "annotation/results_phase3/s35_diversity/abundance/genome_x_sample_relative_abundance.tsv")
MANIFEST_F  <- file.path(BASE, "phase_4_SV_analysis/working/samples_manifest.tsv")
HOSTRANGE_F <- file.path(BASE, "annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv")

# *** CONFIRM *** GTDB-Tk summary files used for the novelty calls.
# Locate them with:  find /disk252/qamar -name "gtdbtk.*.summary.tsv" 2>/dev/null
GTDBTK_FILES <- c(
  file.path(BASE, "CONFIRM_PATH/gtdbtk.bac120.summary.tsv"),
  file.path(BASE, "CONFIRM_PATH/gtdbtk.ar53.summary.tsv")
)

OUT_DIR <- file.path(BASE, "annotation/results_phase3/s79_readbased_host_range")

# Columns in representative_host_range.csv. NA = auto-detect from common names;
# the script stops and lists the real columns if it cannot find one.
HR_ID_COL      <- "representative"
HR_NHOST_COL   <- NA   # number of host species among the cluster's genomes
HR_NGENOME_COL <- NA   # number of genomes in the cluster (singleton = 1)
HR_HOSTS_COL   <- NA   # optional: list of host species (host-identity check)

MIN_SAMPLES_DETECT   <- 2           # primary detection rule
MIN_SAMPLES_SENS     <- c(1, 2, 3)  # sensitivity
N_PER_HOST_SUBSAMPLE <- 40          # equal sampling: hosts with >= this many samples
N_SUBSAMPLE_REPS     <- 200
SEED                 <- 20261003

REPORTED <- list(pct_novel = 69.4, pct_described = 46.3, or_multi = 0.38)
## ============================================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
report_f <- file.path(OUT_DIR, "79_report.txt")
sink(report_f, split = TRUE)
set.seed(SEED)

say <- function(...) cat(sprintf(...), "\n", sep = "")
hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

strip_ext <- function(x) sub("\\.(fa|fasta|fna)(\\.gz)?$", "", x)
norm_host <- function(x) tolower(gsub("[[:space:]_]+", "_", trimws(x)))
split_hosts <- function(x) lapply(strsplit(as.character(x), "[;,|]"), function(v) {
  v <- norm_host(v); unique(v[nzchar(v) & !is.na(v)])
})

pick_col <- function(dt, given, candidates, what, required = TRUE) {
  if (!is.na(given)) {
    if (!given %in% names(dt)) stop(sprintf("Column '%s' (%s) not in %s.\nColumns present: %s",
                                            given, what, HOSTRANGE_F, paste(names(dt), collapse = ", ")))
    return(given)
  }
  hit <- intersect(candidates, names(dt))
  if (length(hit)) return(hit[1])
  if (required) stop(sprintf("Could not find the %s column in %s.\nColumns present: %s\n-> Set it in the CONFIG block.",
                             what, HOSTRANGE_F, paste(names(dt), collapse = ", ")))
  NA_character_
}

# Fisher's exact test, novel x single-host. OR_single > 1 = novel more host-specific.
fisher_row <- function(novel, single, label) {
  ok <- !is.na(novel) & !is.na(single)
  novel <- novel[ok]; single <- single[ok]
  tab <- table(factor(novel, c(TRUE, FALSE)), factor(single, c(TRUE, FALSE)))
  ft  <- fisher.test(tab)
  data.table(comparison = label,
             n_novel = sum(novel), n_described = sum(!novel),
             pct_single_novel = round(100 * mean(single[novel]), 1),
             pct_single_described = round(100 * mean(single[!novel]), 1),
             OR_single = unname(ft$estimate),
             OR_single_low = ft$conf.int[1], OR_single_high = ft$conf.int[2],
             OR_multi = 1 / unname(ft$estimate),
             OR_multi_low = 1 / ft$conf.int[2], OR_multi_high = 1 / ft$conf.int[1],
             p = ft$p.value)
}
fmt_or <- function(r) sprintf("%.1f%% vs %.1f%% single-host; OR_single %.2f (%.2f-%.2f) [= OR_multi %.2f]; p = %s; n = %d novel / %d described",
                              r$pct_single_novel, r$pct_single_described, r$OR_single, r$OR_single_low,
                              r$OR_single_high, r$OR_multi, format.pval(r$p, digits = 2), r$n_novel, r$n_described)

for (f in c(ABUND_F, MANIFEST_F, HOSTRANGE_F)) if (!file.exists(f)) stop("Missing input: ", f)
if (!any(file.exists(GTDBTK_FILES))) stop(
  "No GTDB-Tk summary file found. Set GTDBTK_FILES in the CONFIG block.\n",
  "Locate with:  find /disk252/qamar -name 'gtdbtk.*.summary.tsv' 2>/dev/null")

## ---------------------------------------------------------------------------
## A. Assembly-based host range (representative_host_range.csv)
## ---------------------------------------------------------------------------
hdr("A. Inputs")
hr <- fread(HOSTRANGE_F)
id_col    <- pick_col(hr, HR_ID_COL, c("representative", "genome", "rep_genome"), "representative ID")
ngen_col  <- pick_col(hr, HR_NGENOME_COL, c("n_genomes", "cluster_size", "n_members", "n_member_genomes",
                                            "n_genomes_in_cluster", "size"), "cluster size")
nhost_col <- pick_col(hr, HR_NHOST_COL, c("n_host_species", "n_hosts", "host_range", "n_host",
                                          "num_host_species"), "host count", required = FALSE)
hosts_col <- pick_col(hr, HR_HOSTS_COL, c("host_species_list", "host_list", "hosts", "host_species"),
                      "host list", required = FALSE)
say("representative_host_range.csv: %d rows | id = %s | cluster size = %s | host count = %s | host list = %s",
    nrow(hr), id_col, ngen_col, nhost_col, hosts_col)

asm <- data.table(rep_id = strip_ext(hr[[id_col]]), n_genomes = as.integer(hr[[ngen_col]]))
if (!is.na(hosts_col)) asm[, asm_hosts := split_hosts(hr[[hosts_col]])]
if (!is.na(nhost_col)) {
  asm[, asm_n_hosts := as.integer(hr[[nhost_col]])]
} else if (!is.na(hosts_col)) {
  asm[, asm_n_hosts := lengths(asm_hosts)]
  if (max(asm$asm_n_hosts) < 2) stop("Host list column '", hosts_col, "' never has >1 host; it is probably a ",
                                     "single primary-host column. Set HR_NHOST_COL explicitly.")
} else stop("Neither a host-count nor a host-list column was found; set HR_NHOST_COL.")
if (anyDuplicated(asm$rep_id)) stop("Duplicate representative IDs in ", HOSTRANGE_F)
if (any(asm$n_genomes < 1, na.rm = TRUE) || anyNA(asm$n_genomes)) stop("Cluster size column has NA or <1 values.")
say("Assembly host counts: %s", paste(capture.output(print(table(pmin(asm$asm_n_hosts, 5)))), collapse = " "))
say("(5 = 5 or more host species)")

## ---------------------------------------------------------------------------
## B. Novelty from GTDB-Tk classification of the representative
## ---------------------------------------------------------------------------
gt <- rbindlist(lapply(GTDBTK_FILES[file.exists(GTDBTK_FILES)], fread,
                       select = c("user_genome", "classification")))
gt <- unique(gt[, .(rep_id = strip_ext(user_genome), classification)], by = "rep_id")
gt[, `:=`(s_ok = grepl("s__[^;]+", classification),
          g_ok = grepl("g__[^;]+", classification),
          f_ok = grepl("f__[^;]+", classification))]
gt[, novelty_level := fifelse(s_ok, "described_species",
                      fifelse(g_ok, "novel_species",
                      fifelse(f_ok, "novel_genus", "novel_family_or_deeper")))]
gt[, novel := !s_ok]   # novel-species candidate = species field unassigned
say("GTDB-Tk genomes: %d | novelty levels: %s", nrow(gt),
    paste(names(table(gt$novelty_level)), table(gt$novelty_level), sep = "=", collapse = ", "))

## ---------------------------------------------------------------------------
## 0. Reproduce the reported assembly-based result (before touching reads)
## ---------------------------------------------------------------------------
hdr("0. Reproduction check: assembly-based host range, non-singleton clusters")
base0 <- merge(asm, gt[, .(rep_id, novel, novelty_level)], by = "rep_id")
say("Representatives with both host range and GTDB-Tk call: %d of %d", nrow(base0), nrow(asm))
if (nrow(base0) < 0.95 * nrow(asm)) stop("Fewer than 95% of representatives matched a GTDB-Tk call; check ID formats:\n",
                                         "  host range: ", paste(head(asm$rep_id, 3), collapse = ", "),
                                         "\n  GTDB-Tk:    ", paste(head(gt$rep_id, 3), collapse = ", "))
ns0 <- base0[n_genomes > 1]
r_broad  <- fisher_row(ns0$novel, ns0$asm_n_hosts == 1,
                       "Assembly, non-singletons; novel = species unassigned (incl. novel genus/family)")
sub_s    <- ns0[novelty_level %in% c("novel_species", "described_species")]
r_strict <- fisher_row(sub_s$novelty_level == "novel_species", sub_s$asm_n_hosts == 1,
                       "Assembly, non-singletons; novel = species-level only (genus known)")
repro <- rbind(r_broad, r_strict)
say("Broad definition : %s", fmt_or(r_broad))
say("Strict definition: %s", fmt_or(r_strict))
say("Reported         : %.1f%% vs %.1f%% single-host; OR_multi %.2f",
    REPORTED$pct_novel, REPORTED$pct_described, REPORTED$or_multi)
d_broad  <- abs(r_broad$pct_single_novel  - REPORTED$pct_novel) + abs(r_broad$pct_single_described  - REPORTED$pct_described)
d_strict <- abs(r_strict$pct_single_novel - REPORTED$pct_novel) + abs(r_strict$pct_single_described - REPORTED$pct_described)
NOVEL_DEF <- if (d_strict < d_broad) "strict" else "broad"
say("-> Using the '%s' novelty definition (closer to the reported numbers; deviation %.1f points).",
    NOVEL_DEF, min(d_broad, d_strict))
if (min(d_broad, d_strict) > 0.3) warning("Neither definition reproduces 69.4/46.3 within 0.3 points. ",
  "Check that HOSTRANGE_F and GTDBTK_FILES are the files behind the reported result before using anything below.")
fwrite(repro, file.path(OUT_DIR, "79_reproduction_check.csv"))

## ---------------------------------------------------------------------------
## C. Reads: CoverM matrix + sample hosts
## ---------------------------------------------------------------------------
hdr("C. Read-based detection (existing CoverM matrix, no rerun)")
man <- unique(fread(MANIFEST_F)[, .(sample_id, host = species)])
if (anyDuplicated(man$sample_id)) stop("A sample_id maps to more than one host in ", MANIFEST_F)

ab <- fread(ABUND_F)
if (!"Genome" %in% names(ab)) stop("No 'Genome' column in ", ABUND_F)
ab <- ab[Genome != "unmapped"]
gen_ids <- strip_ext(ab$Genome)
samp <- setdiff(names(ab), "Genome")
mat <- as.matrix(ab[, ..samp]); storage.mode(mat) <- "double"
rm(ab); invisible(gc())

host_of <- man$host[match(samp, man$sample_id)]
if (anyNA(host_of)) {
  warning(sum(is.na(host_of)), " matrix samples have no host in the manifest and are dropped.")
  mat <- mat[, !is.na(host_of), drop = FALSE]; samp <- samp[!is.na(host_of)]; host_of <- host_of[!is.na(host_of)]
}
say("Matrix: %d genomes x %d samples, %d hosts", nrow(mat), ncol(mat), length(unique(host_of)))

det  <- mat > 0
prev <- rowSums(det)
mean_all     <- rowMeans(mat)
mean_present <- ifelse(prev > 0, rowSums(mat) / prev, NA_real_)
rm(mat); invisible(gc())

hosts <- sort(unique(host_of))
H <- vapply(hosts, function(h) as.numeric(host_of == h), numeric(length(host_of)))
cnt <- det %*% H                                   # genomes x hosts: samples with detection
host_power <- data.table(host = hosts, n_samples = colSums(H),
                         n_species_detected_min2 = colSums(cnt >= MIN_SAMPLES_DETECT))
print(host_power[order(-n_samples)])
say("Hosts with < %d samples can never count under the primary rule: %s", MIN_SAMPLES_DETECT,
    paste(host_power[n_samples < MIN_SAMPLES_DETECT, host], collapse = ", "))
fwrite(host_power, file.path(OUT_DIR, "79_hosts_samples_and_detection.csv"))

sp <- data.table(rep_id = gen_ids, row_idx = seq_along(gen_ids), prevalence = prev,
                 mean_abund_all = mean_all, mean_abund_present = mean_present)
for (m in MIN_SAMPLES_SENS) set(sp, j = paste0("read_n_hosts_min", m), value = rowSums(cnt >= m))
w <- which(cnt >= MIN_SAMPLES_DETECT, arr.ind = TRUE)
read_hosts <- unname(split(norm_host(hosts[w[, 2]]), factor(w[, 1], levels = seq_len(nrow(cnt)))))
sp[, read_hosts := read_hosts]

n_mat <- nrow(sp)
sp <- merge(sp, base0, by = "rep_id")
say("Species in matrix: %d | with host range + GTDB-Tk call: %d", n_mat, nrow(sp))
if (nrow(sp) < 0.95 * n_mat) stop("Fewer than 95% of matrix genomes matched; check ID formats:\n",
                                  "  matrix: ", paste(head(gen_ids, 3), collapse = ", "),
                                  "\n  host range: ", paste(head(asm$rep_id, 3), collapse = ", "))

RN <- paste0("read_n_hosts_min", MIN_SAMPLES_DETECT)
sp[, `:=`(asm_single  = asm_n_hosts == 1,
          read_detect = get(RN) > 0,
          read_single = get(RN) == 1,
          singleton   = n_genomes == 1)]
sp_all_levels <- copy(sp)   # keeps every novelty level, for the depth gradient

# Apply the novelty definition that reproduced the reported result
if (NOVEL_DEF == "strict") {
  sp <- sp[novelty_level %in% c("novel_species", "described_species")]
  sp[, novel := novelty_level == "novel_species"]
}

## ---------------------------------------------------------------------------
## 1. Who is detected by reads at all?
## ---------------------------------------------------------------------------
hdr(sprintf("1. Read detection by novelty (rule: >0 abundance in >=%d samples of a host)", MIN_SAMPLES_DETECT))
det_tab <- sp[, .(n = .N,
                  pct_undetected = round(100 * mean(!read_detect), 1),
                  pct_read_single_of_all = round(100 * mean(read_single), 1),
                  pct_read_multi_of_all  = round(100 * mean(get(RN) >= 2), 1),
                  median_prevalence = as.numeric(median(prevalence)),
                  median_abund_present = signif(median(mean_abund_present, na.rm = TRUE), 3)),
              by = .(set = fifelse(singleton, "singleton", "non-singleton"), novel)][order(set, novel)]
print(det_tab)
fwrite(det_tab, file.path(OUT_DIR, "79_detection_by_novelty.csv"))
if (sp[!(singleton), sum(!read_detect)] == 0) {
  say("Every non-singleton species is detected by reads in at least one host.")
} else {
  r_undet <- fisher_row(sp[!(singleton)]$novel, !sp[!(singleton)]$read_detect,
                        "Undetected by reads (non-singletons)")
  say("Novel vs described, non-singletons, odds of being UNDETECTED by reads: OR %.2f (%.2f-%.2f), p = %s",
      r_undet$OR_single, r_undet$OR_single_low, r_undet$OR_single_high, format.pval(r_undet$p, digits = 2))
  say("(A large OR here is itself the detection-bias signature: novel species are harder to see.)")
}

## ---------------------------------------------------------------------------
## 2. Key test: novel vs described, read-based host range
## ---------------------------------------------------------------------------
hdr("2. Novel vs described: single-host by assembly and by reads")
core <- sp[!(singleton) & read_detect]   # same subset for both classifications
tests <- rbind(
  fisher_row(sp[!(singleton)]$novel, sp[!(singleton)]$asm_single,
             "A. Assembly-based, all non-singletons (= reported analysis)"),
  fisher_row(core$novel, core$asm_single,
             sprintf("B. Assembly-based, non-singletons detected by reads (same set as C)")),
  fisher_row(core$novel, core$read_single,
             sprintf("C. READ-BASED, non-singletons, >=%d samples per host [PRIMARY]", MIN_SAMPLES_DETECT))
)
for (m in setdiff(MIN_SAMPLES_SENS, MIN_SAMPLES_DETECT)) {
  col <- paste0("read_n_hosts_min", m)
  sub <- sp[!(singleton) & get(col) > 0]
  tests <- rbind(tests, fisher_row(sub$novel, sub[[col]] == 1,
                 sprintf("D. Read-based, non-singletons, >=%d sample(s) per host", m)))
}
allsp <- sp[(read_detect)]
tests <- rbind(tests, fisher_row(allsp$novel, allsp$read_single,
               "E. Read-based, ALL species incl. singletons (reads give singletons a host range)"))
for (i in seq_len(nrow(tests))) say("%s\n   %s", tests$comparison[i], fmt_or(tests[i]))
fwrite(tests, file.path(OUT_DIR, "79_novelty_host_specificity_tests.csv"))

grad <- sp_all_levels[!(singleton) & read_detect, .(n = .N,
             pct_single_assembly = round(100 * mean(asm_single), 1),
             pct_single_reads    = round(100 * mean(read_single), 1)),
           by = novelty_level][order(match(novelty_level,
             c("described_species", "novel_species", "novel_genus", "novel_family_or_deeper")))]
say("\nNovelty-depth gradient (non-singletons detected by reads):")
print(grad)
fwrite(grad, file.path(OUT_DIR, "79_novelty_depth_gradient.csv"))

## ---------------------------------------------------------------------------
## 3. Logistic models: does abundance explain the novelty effect?
## ---------------------------------------------------------------------------
hdr("3. Logistic regression: single-host ~ novel + log10 abundance (non-singletons detected by reads)")
say("Two abundance covariates:")
say("  mean_abund_all     = mean relative abundance over all %d samples (as proposed). It is partly a", length(samp))
say("                       CONSEQUENCE of host range (a one-host species has zeros in every other host),")
say("                       so it over-adjusts and biases the novelty OR toward 1.")
say("  mean_abund_present = mean abundance in samples where detected (abundance when present).")
say("                       Not mechanically tied to host breadth -> the cleaner detection control.")
fit_one <- function(outcome, rhs, label) {
  f <- as.formula(paste(outcome, "~", rhs))
  g <- glm(f, data = core, family = binomial)
  ci <- suppressMessages(confint.default(g))
  cf <- summary(g)$coefficients
  data.table(outcome = outcome, model = label, n = nobs(g),
             OR_novel = exp(cf["novelTRUE", 1]),
             OR_low = exp(ci["novelTRUE", 1]), OR_high = exp(ci["novelTRUE", 2]),
             p_novel = cf["novelTRUE", 4],
             OR_per_log10_abund = if (nrow(cf) > 2) exp(cf[3, 1]) else NA_real_,
             p_abund = if (nrow(cf) > 2) cf[3, 4] else NA_real_,
             AIC = AIC(g))
}
logit <- rbindlist(lapply(c("asm_single", "read_single"), function(y) rbind(
  fit_one(y, "novel", "M0 novelty only"),
  fit_one(y, "novel + log10(mean_abund_all)", "M1 + log10 mean abundance, all samples"),
  fit_one(y, "novel + log10(mean_abund_present)", "M2 + log10 abundance when present [recommended]")
)))
print(logit[, .(outcome, model, n, OR_novel = round(OR_novel, 2), OR_low = round(OR_low, 2),
                OR_high = round(OR_high, 2), p_novel = signif(p_novel, 2),
                OR_per_log10_abund = round(OR_per_log10_abund, 2))])
fwrite(logit, file.path(OUT_DIR, "79_logistic_novelty_abundance.csv"))
abund_cmp <- core[, .(median_abund_present = signif(median(mean_abund_present), 3),
                      median_abund_all = signif(median(mean_abund_all), 3)), by = novel]
say("\nMedian abundance by novelty (same subset):"); print(abund_cmp)

## ---------------------------------------------------------------------------
## 4. Agreement between assembly- and read-based host specificity
## ---------------------------------------------------------------------------
hdr("4. Agreement: assembly-based vs read-based (non-singletons)")
nsx <- sp[!(singleton)]
nsx[, asm_class  := fifelse(asm_single, "single-host", "multi-host")]
nsx[, read_class := fifelse(!read_detect, "undetected", fifelse(read_single, "single-host", "multi-host"))]
agree <- nsx[, .N, by = .(novelty = fifelse(novel, "novel", "described"), asm_class, read_class)]
agree_all <- nsx[, .N, by = .(asm_class, read_class)][, novelty := "all"]
agree <- rbind(agree_all, agree)[order(novelty, asm_class, read_class)]
agree[, pct_of_asm_class := round(100 * N / sum(N), 1), by = .(novelty, asm_class)]
print(dcast(agree, novelty + asm_class ~ read_class, value.var = "N", fill = 0))
fwrite(agree, file.path(OUT_DIR, "79_agreement_assembly_vs_reads.csv"))

kappa2 <- function(a, b) {
  t <- table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))); n <- sum(t)
  po <- sum(diag(t)) / n; pe <- sum(rowSums(t) * colSums(t)) / n^2
  c(agreement = 100 * po, kappa = (po - pe) / (1 - pe))
}
k_all <- kappa2(core$asm_single, core$read_single)
k_nov <- kappa2(core[(novel)]$asm_single, core[(novel)]$read_single)
k_des <- kappa2(core[!(novel)]$asm_single, core[!(novel)]$read_single)
say("Detected non-singletons: agreement %.1f%%, Cohen's kappa %.2f (novel %.1f%% / %.2f; described %.1f%% / %.2f)",
    k_all[1], k_all[2], k_nov[1], k_nov[2], k_des[1], k_des[2])
say("Assembly single-host species that reads place in >=2 hosts: novel %.1f%%, described %.1f%%",
    100 * core[(novel) & asm_single, mean(!read_single)], 100 * core[!(novel) & asm_single, mean(!read_single)])

if ("asm_hosts" %in% names(sp)) {
  hid <- core[(asm_single) & lengths(asm_hosts) == 1]
  hid[, asm_host := vapply(asm_hosts, `[`, "", 1)]
  hid[, in_reads := mapply(function(h, r) h %in% r, asm_host, read_hosts)]
  hid[, reads_only_asm_host := mapply(function(h, r) length(r) == 1 && r == h, asm_host, read_hosts)]
  say("Host identity, assembly single-host species detected by reads (n = %d): assembly host among read hosts %.1f%%; reads place it ONLY in the assembly host %.1f%%",
      nrow(hid), 100 * mean(hid$in_reads), 100 * mean(hid$reads_only_asm_host))
  if (mean(hid$in_reads) < 0.5) warning("Assembly host is rarely among read hosts: host labels in the two files ",
                                        "may not be harmonised (e.g. 'Ovis aries' vs 'sheep'). Check before quoting.")
}

supp_cols <- c("rep_id", "novelty_level", "novel", "n_genomes", "asm_n_hosts",
               paste0("read_n_hosts_min", MIN_SAMPLES_SENS), "prevalence",
               "mean_abund_all", "mean_abund_present")
supp <- sp[, ..supp_cols]
supp[, read_hosts_primary_rule := vapply(sp$read_hosts, paste, "", collapse = ";")]
supp[, `:=`(mean_abund_all = signif(mean_abund_all, 4), mean_abund_present = signif(mean_abund_present, 4))]
setnames(supp, c("rep_id", "n_genomes", "asm_n_hosts", "prevalence", "mean_abund_all", "mean_abund_present"),
         c("representative", "n_genomes_in_cluster", "assembly_n_hosts", "n_samples_detected",
           "mean_rel_abund_all_samples", "mean_rel_abund_when_present"))
fwrite(supp, file.path(OUT_DIR, "79_TableS_species_host_range_assembly_vs_reads.csv"))

## ---------------------------------------------------------------------------
## 5. Equal sampling effort per host
## ---------------------------------------------------------------------------
hdr(sprintf("5. Equal sampling: %d samples drawn per host, hosts with >= %d samples, %d replicates",
            N_PER_HOST_SUBSAMPLE, N_PER_HOST_SUBSAMPLE, N_SUBSAMPLE_REPS))
elig <- host_power[n_samples >= N_PER_HOST_SUBSAMPLE, host]
say("Hosts in the panel (%d): %s", length(elig), paste(elig, collapse = ", "))
ns_rows <- sp[!(singleton), row_idx]
ns_novel <- sp[!(singleton), novel]

panel_test <- function(cols, label) {
  Hs <- vapply(elig, function(h) as.numeric(host_of[cols] == h), numeric(length(cols)))
  c_sub <- det[ns_rows, cols, drop = FALSE] %*% Hs
  nh <- rowSums(c_sub >= MIN_SAMPLES_DETECT)
  keep <- nh > 0
  r <- fisher_row(ns_novel[keep], nh[keep] == 1, label)
  r[, n_detected := sum(keep)]
  r
}
panel_cols <- which(host_of %in% elig)
full_panel <- panel_test(panel_cols, "Panel hosts, all their samples (no subsampling)")
say("%s\n   %s", full_panel$comparison, fmt_or(full_panel))

idx_by_host <- lapply(elig, function(h) which(host_of == h))
reps <- rbindlist(lapply(seq_len(N_SUBSAMPLE_REPS), function(b) {
  cols <- unlist(lapply(idx_by_host, function(ix) ix[sample.int(length(ix), N_PER_HOST_SUBSAMPLE)]))
  panel_test(cols, paste0("rep", b))
}))
fwrite(reps, file.path(OUT_DIR, "79_equal_sampling_replicates.csv"))
eq_sum <- data.table(
  hosts = length(elig), samples_per_host = N_PER_HOST_SUBSAMPLE, replicates = N_SUBSAMPLE_REPS,
  median_pct_single_novel = median(reps$pct_single_novel),
  median_pct_single_described = median(reps$pct_single_described),
  median_OR_single = median(reps$OR_single),
  OR_single_2.5 = quantile(reps$OR_single, 0.025), OR_single_97.5 = quantile(reps$OR_single, 0.975),
  pct_reps_CI_above_1 = 100 * mean(reps$OR_single_low > 1),
  pct_reps_CI_below_1 = 100 * mean(reps$OR_single_high < 1),
  median_n_detected = median(reps$n_detected))
print(eq_sum)
fwrite(full_panel, file.path(OUT_DIR, "79_equal_sampling_full_panel.csv"))
fwrite(eq_sum, file.path(OUT_DIR, "79_equal_sampling_summary.csv"))

## ---------------------------------------------------------------------------
## VERDICT (explicit rule, mapped to the three outcomes in the analysis plan)
## ---------------------------------------------------------------------------
hdr("VERDICT")
A <- tests[1]; B <- tests[2]; C <- tests[3]
m2_asm  <- logit[outcome == "asm_single"  & grepl("^M2", model)]
m2_read <- logit[outcome == "read_single" & grepl("^M2", model)]
read_sig <- C$OR_single_low > 1
adj_sig  <- m2_asm$OR_low > 1 && m2_read$OR_low > 1
eq_sig   <- eq_sum$pct_reps_CI_above_1 >= 95
shrink   <- log(C$OR_single) / log(B$OR_single)   # share of the effect retained on the same subset

say("Assembly, same subset (B): OR_single %.2f | Read-based (C): OR_single %.2f -> %.0f%% of the log-OR retained",
    B$OR_single, C$OR_single, 100 * shrink)
say("Read-based CI excludes 1 (novel more host-specific):       %s", read_sig)
say("Novelty OR survives abundance-when-present (both outcomes): %s  [assembly %.2f (%.2f-%.2f); reads %.2f (%.2f-%.2f)]",
    adj_sig, m2_asm$OR_novel, m2_asm$OR_low, m2_asm$OR_high, m2_read$OR_novel, m2_read$OR_low, m2_read$OR_high)
say("Equal sampling: >=95%% of replicates with CI above 1:      %s  (%.0f%%)", eq_sig, eq_sum$pct_reps_CI_above_1)

verdict <- if (!read_sig || C$OR_single <= 1) {
  "DISAPPEARS: with read-based host range novelty no longer predicts host specificity. Title and abstract need to change."
} else if (adj_sig && eq_sig && shrink >= 0.5) {
  "HOLDS: novelty predicts host specificity with read-based host range, after abundance control and with equal sampling. Keep the title; add the Methods sentence, two Results sentences and the supplementary table."
} else if (adj_sig && eq_sig) {
  "SHRINKS BUT HOLDS: the effect is smaller with reads (<50% of the log-OR retained) but significant in every check. Keep the claim and report both estimates."
} else {
  "MIXED: the read-based effect is significant but at least one check fails (see above). Report every estimate and soften the claim; do not call it confirmed."
}
say("\n%s", verdict)
say("\nCaveats to state: (i) at 95%% identity some reads from close relatives can still map; the 10%% breadth floor and competitive mapping limit this. (ii) About three-quarters of samples are faecal, so absence in faeces is not absence in the rumen. (iii) Hosts with few samples have little power (see 79_hosts_samples_and_detection.csv). (iv) Three catalogue hosts have no reads, so species assembled only from them cannot be confirmed by reads.")
say("\nAll outputs in: %s", OUT_DIR)
sink()
