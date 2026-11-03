#!/usr/bin/env Rscript
###############################################################################
# 58_ruc_lineage_propagated.R
#
# Regenerates every number in the two RUC-lineage Results paragraphs from the
# CORRECTED annotations:
#   (a) CAZyme paragraph (old 35_ + substrate percentages): restored dbCAN IDs,
#       propagated calls, all domains (54_ outputs);
#   (b) thiamine paragraph (old 36j): propagated KO calls (50_ output).
#
# Tests: Fisher's exact, lineage vs rest of catalogue (as in the draft), for
# every family carried by >= 30 genomes catalogue-wide, BH across families;
# substrates reported as percentages in the lineage vs each tier; thiamine
# genes vs rest of catalogue AND vs host-specific genomes outside the lineage.
#
# OUTPUTS (OUT_DIR)
#   table_s58_ruc_family_fisher.csv       all tested families (lineage vs rest)
#   table_s58_ruc_substrates.csv          substrate % in lineage / core / shared / host
#   table_s58_ruc_thiamine.csv            thiC, thiE, thiL, thiL+/thiC- tests
#   table_s58_ruc_thiamine_by_host.csv    thiL and thiC % by host species within the lineage
###############################################################################

suppressPackageStartupMessages({ library(data.table) })
P <- function(var, default) Sys.getenv(var, unset = default)
R3 <- "/disk252/qamar/metagenome_analysis/annotation/results_phase3"
LINEAGE_F   <- P("LINEAGE_F", file.path(R3, "s32_taxonomy_novelty/lineage_full_classification_table.csv"))
SHARING_F   <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
FAM_F       <- P("FAM_F", file.path(R3, "s54_cazyme_propagated/54_cazy_family_presence_propagated.tsv"))
SUB_F       <- P("SUB_F", file.path(R3, "s54_cazyme_propagated/54_cazy_substrate_presence_propagated.tsv"))
KO_F        <- P("KO_F", file.path(R3, "s50_subpathway_propagated/50_marker_ko_presence_propagated.tsv"))
OUT_DIR     <- P("OUT_DIR", file.path(R3, "s58_ruc_lineage_propagated"))
MIN_CARRIERS <- 30
KEY_FAMS <- c("GH190", "GT121", "GT26")
KEY_SUBS <- c("peptidoglycan", "host glycan", "chitin")
THI <- c(thiC = "K03147", thiE = "K00788", thiL = "K00949")
for (f in c(LINEAGE_F, SHARING_F, FAM_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
hr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
options(width = 200)

sh <- fread(SHARING_F)
setnames(sh, "representative", "genome")
lin <- unique(fread(LINEAGE_F)$genome)
cat(sprintf("Lineage genomes: %d (expected 76) | in catalogue: %d\n", length(lin), sum(lin %chin% sh$genome)))
if (!all(lin %chin% sh$genome)) warning("Some lineage genomes are not catalogue representatives: ", paste(head(setdiff(lin, sh$genome)), collapse = ", "))
sh[, in_lineage := genome %chin% lin]
cat("Lineage genomes by sharing tier:\n"); print(sh[in_lineage == TRUE, .N, by = sharing_category])

fisher_row <- function(has, label, ref = !sh$in_lineage) {
  a <- sum(has & sh$in_lineage); b <- sum(!has & sh$in_lineage); c <- sum(has & ref); d <- sum(!has & ref)
  ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
  data.table(trait = label, n_lineage = a, pct_lineage = round(100 * a / (a + b), 1), n_ref = c,
             pct_ref = round(100 * c / (c + d), 1), OR = unname(ft$estimate), OR_lo = ft$conf.int[1],
             OR_hi = ft$conf.int[2], p = ft$p.value)
}

# =============================================================================
# (a) CAZy families
# =============================================================================
hr("(a) CAZy families: lineage vs rest of catalogue")
fam <- unique(fread(FAM_F)[, .(genome, cazy_family)])
cat(sprintf("Lineage genomes with >= 1 CAZy family: %d of %d\n", fam[genome %chin% lin, uniqueN(genome)], length(lin)))
carriers <- fam[, .N, by = cazy_family][N >= MIN_CARRIERS, cazy_family]
fl <- split(fam$genome, fam$cazy_family)
ff <- rbindlist(lapply(carriers, function(f) fisher_row(sh$genome %chin% fl[[f]], f)))
ff[, p_adj_BH := p.adjust(p, "BH")]
setorder(ff, p_adj_BH)
fwrite(ff, file.path(OUT_DIR, "table_s58_ruc_family_fisher.csv"))
cat(sprintf("Families tested: %d | BH p < 0.05: %d (enriched %d, depleted %d)\n", nrow(ff), ff[p_adj_BH < 0.05, .N],
            ff[p_adj_BH < 0.05 & OR > 1, .N], ff[p_adj_BH < 0.05 & OR < 1, .N]))
fmt <- function(d) d[, .(trait, pct_lineage, pct_ref, OR = signif(OR, 3), OR_lo = signif(OR_lo, 3), OR_hi = signif(OR_hi, 3), p_adj_BH = signif(p_adj_BH, 3))]
cat("\nDraft's key families:\n"); print(fmt(ff[trait %in% KEY_FAMS]))
cat("\nTop 15 enriched in the lineage (by OR, BH p < 0.05):\n"); print(fmt(head(ff[p_adj_BH < 0.05 & OR > 1][order(-OR)], 15)))
cat("\nTop 15 by adjusted p:\n"); print(fmt(head(ff, 15)))

# =============================================================================
# Substrates: lineage vs tiers (the draft compares lineage % with core/broad %)
# =============================================================================
if (file.exists(SUB_F)) {
  hr("Substrates: % of genomes, lineage vs tiers")
  sb <- unique(fread(SUB_F)[, .(genome, substrate)])
  grp <- rbind(sh[in_lineage == TRUE, .(genome, group = "RUC lineage")], sh[, .(genome, group = sharing_category)])
  tot <- grp[, .(n_total = .N), by = group]
  x <- merge(sb, grp, by = "genome", allow.cartesian = TRUE)[, .(n = .N), by = .(substrate, group)]
  x <- merge(x, tot, by = "group"); x[, pct := round(100 * n / n_total, 1)]
  w <- dcast(x, substrate ~ group, value.var = "pct", fill = 0)
  fwrite(w, file.path(OUT_DIR, "table_s58_ruc_substrates.csv"))
  print(w[substrate %in% KEY_SUBS]); cat("(all substrates in table_s58_ruc_substrates.csv)\n")
}

# =============================================================================
# (b) Thiamine
# =============================================================================
if (file.exists(KO_F)) {
  hr("(b) Thiamine genes (propagated KO calls)")
  ko <- fread(KO_F)
  gcol <- intersect(c("genome_id", "genome"), names(ko))[1]
  ko <- unique(ko[get("ko") %chin% THI, .(genome = get(gcol), ko)])
  has <- lapply(THI, function(k) sh$genome %chin% ko[ko == k, genome])
  has$`thiL+/thiC-` <- has$thiL & !has$thiC
  host_ref <- !sh$in_lineage & sh$sharing_category == "Host-specific (1 host species)"
  th <- rbind(
    rbindlist(lapply(names(has), function(n) fisher_row(has[[n]], n)))[, comparison := "vs rest of catalogue"],
    rbindlist(lapply(names(has), function(n) fisher_row(has[[n]], n, ref = host_ref)))[, comparison := "vs host-specific genomes outside the lineage"])
  th[, p_adj_BH := p.adjust(p, "BH"), by = comparison]
  fwrite(th, file.path(OUT_DIR, "table_s58_ruc_thiamine.csv"))
  print(th[, .(comparison, trait, n_lineage, pct_lineage, pct_ref, OR = signif(OR, 3), OR_lo = signif(OR_lo, 3),
               OR_hi = signif(OR_hi, 3), p = signif(p, 3), p_adj_BH = signif(p_adj_BH, 3))])
  if ("host_species_list" %in% names(sh)) {
    hs <- sh[in_lineage == TRUE, .(host = trimws(unlist(strsplit(host_species_list, "[;,|]")))), by = genome][host != ""]
    hs[, `:=`(thiL = genome %chin% ko[ko == THI["thiL"], genome], thiC = genome %chin% ko[ko == THI["thiC"], genome])]
    byh <- hs[, .(n_genomes = .N, pct_thiL = round(100 * mean(thiL), 1), pct_thiC = round(100 * mean(thiC), 1)), by = host][order(-n_genomes)]
    fwrite(byh, file.path(OUT_DIR, "table_s58_ruc_thiamine_by_host.csv"))
    cat(sprintf("\nthiL within host species (lineage): %.1f-%.1f%% across %d host species\n", min(byh$pct_thiL), max(byh$pct_thiL), nrow(byh)))
    print(byh)
  } else cat("No host_species_list column — per-host thiL skipped.\n")
} else cat("KO table not found (", KO_F, ") — thiamine part skipped.\n", sep = "")
cat("\nOutputs in ", OUT_DIR, "\n", sep = "")
