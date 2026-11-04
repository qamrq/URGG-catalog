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
#  SCRIPT  : 31a_target_tree_itol.R  -  iTOL datasets for the target-genome tree (Fig. 5a)
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s32_taxonomy_novelty/target_ml_tree_extrefs/
#            annotation/results_phase3/s32_taxonomy_novelty/lineage_full_classification_table.csv
#  OUTPUT  : iTOL dataset files
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(ape); library(dplyr); library(readr); library(stringr) })

S <- paste0(URGG_PROJECT, "/annotation/results_phase3")
env <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
PRES     <- env("PRES",     file.path(S, "s54_cazyme_propagated/54_cazy_family_presence_propagated.tsv"))
LABEL_F  <- env("LABEL_CSV", paste0(URGG_PROJECT, "/phase_3_analysis/scripts/results/lineage_full_classification_table.csv"))
WORK     <- env("WORK",     file.path(S, "s32_taxonomy_novelty/target_ml_tree_extrefs"))
KEY      <- env("KEY",      "raw")
if (!KEY %in% c("raw", "public")) stop("KEY must be raw or public")
TREE_F   <- env("TREE",     file.path(WORK, ifelse(KEY == "raw", "ml_tree_rooted_rawIDs.nwk", "ml_tree_rooted_publicID.nwk")))
TARGET   <- env("TARGET",   "CNP0001390__Sheep__Feces__T1900044680_TS20190146_bin.8")
OLD_DIR  <- env("OLD_DIR",  paste0(URGG_PROJECT, "/phase_3_analysis/scripts/results"))
TAG      <- ifelse(KEY == "raw", "RAW_ID", "PUBLIC_ID")
ROLES_F  <- env("ROLES",    file.path(WORK, "tip_roles.tsv"))
GTDB_TAX <- env("GTDB_TAX", "")
OUT_DIR  <- env("OUT_DIR",  file.path(WORK, "itol"))
CHECK    <- env("CHECK",    "die")
PRUNE76  <- env("PRUNE76",  "1")
FULL     <- env("FULL",     "1")
SUP_MIN  <- as.numeric(env("SUP_MIN", "95"))

FAM_PURPLE <- c("GH190", "GT121", "GT26")
FAM_OCHRE  <- c("PL6", "PL15", "PL21", "PL29", "PL34")
FAM_TEAL   <- c("AA6")
FAMS <- c(FAM_PURPLE, FAM_OCHRE, FAM_TEAL)
COL  <- c(setNames(rep("#7B3294", length(FAM_PURPLE)), FAM_PURPLE),
          setNames(rep("#C58B00", length(FAM_OCHRE)),  FAM_OCHRE),
          setNames(rep("#1B9E77", length(FAM_TEAL)),   FAM_TEAL))
EXPECT <- c(GH190 = 64, GT121 = 75, GT26 = 67, PL6 = 14, PL15 = 21, PL21 = 30, PL29 = 15, PL34 = 13, AA6 = 34)

stopf <- function(...) stop("CHECKPOINT FAILED: ", sprintf(...), call. = FALSE)
ok    <- function(...) cat("  [ok] ", sprintf(...), "\n", sep = "")

cat("== C1 inputs ==\n")
for (f in c(PRES, LABEL_F, TREE_F, ROLES_F)) if (!file.exists(f)) stopf("missing file: %s", f)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
ok("all inputs found; output dir %s", OUT_DIR)

cat("== C2 lineage label table ==\n")
L <- read_csv(LABEL_F, show_col_types = FALSE)
if (!all(c("genome", "public_id") %in% names(L))) stopf("label table needs columns genome, public_id")
if (nrow(L) != 76 || anyDuplicated(L$genome) || anyDuplicated(L$public_id))
  stopf("expected 76 unique genomes/public_ids, got %d rows (dup genome %d, dup id %d)",
        nrow(L), anyDuplicated(L$genome) > 0, anyDuplicated(L$public_id) > 0)
ok("76 unique genomes and public_ids")

cat("== C3 propagated dbCAN presence ==\n")
P <- read_tsv(PRES, col_types = cols_only(genome = col_character(), cazy_family = col_character()),
              progress = FALSE)
P <- P %>% filter(genome %in% L$genome) %>% distinct(genome, cazy_family)
miss <- setdiff(L$genome, P$genome)
if (length(miss)) stopf("%d of the 76 genomes have no CAZy call at all (first: %s). Check the genome-ID form in %s",
                        length(miss), miss[1], basename(PRES))
ok("all 76 genomes present in the presence table (%d genome-family pairs)", nrow(P))

M <- sapply(FAMS, function(f) as.integer(L$genome %in% P$genome[P$cazy_family == f]))
rownames(M) <- if (KEY == "raw") L$genome else L$public_id
cnt <- colSums(M)
chk <- data.frame(family = FAMS, n_carriers = as.integer(cnt), pct = round(100 * cnt / 76, 1),
                  expected = as.integer(EXPECT[FAMS]), match = as.integer(cnt) == as.integer(EXPECT[FAMS]))
write.csv(chk, file.path(OUT_DIR, "ruc_cazyme_counts_check.csv"), row.names = FALSE)
cat("== C5 carrier counts among the 76 lineage genomes ==\n"); print(chk, row.names = FALSE)
if (!all(chk$match)) {
  msg <- paste0("carrier counts differ from the rerun values for: ", paste(chk$family[!chk$match], collapse = ", "),
                ". The most likely causes are a different presence table (repOnly / min2tools) or a changed genome set.")
  if (CHECK == "die") stopf("%s Re-run with CHECK=warn to write the files anyway.", msg) else warning(msg)
} else ok("all nine families match the rerun values")

f_bin <- file.path(OUT_DIR, paste0("03b_cazyme_binary_RERUN_", TAG, ".txt"))
con <- file(f_bin, "w")
writeLines(c("DATASET_BINARY", "SEPARATOR TAB", "DATASET_LABEL\tCAZyme families (rerun)", "COLOR\t#000000",
             paste0("FIELD_SHAPES\t", paste(rep(2, length(FAMS)), collapse = "\t")),
             paste0("FIELD_LABELS\t", paste(FAMS, collapse = "\t")),
             paste0("FIELD_COLORS\t", paste(COL[FAMS], collapse = "\t")),
             "LEGEND_TITLE\tCAZyme family",
             paste0("LEGEND_SHAPES\t", paste(rep(2, 3), collapse = "\t")),
             "LEGEND_COLORS\t#7B3294\t#C58B00\t#1B9E77",
             "LEGEND_LABELS\tGH190 / GT121 / GT26\tPL6 / PL15 / PL21 / PL29 / PL34\tAA6",
             "DATA"), con)
writeLines(paste(rownames(M), apply(M, 1, paste, collapse = "\t"), sep = "\t"), con)
close(con)
ok("wrote %s (%d rows)", basename(f_bin), nrow(M))

cat("== C4 tree tips and roles ==\n")
tree <- read.tree(TREE_F)
R <- read_tsv(ROLES_F, show_col_types = FALSE)
if (!all(c("id", "role") %in% names(R))) stopf("tip_roles.tsv needs columns id, role")
lin_ids <- if (KEY == "raw") L$genome else L$public_id
tgt_tip <- if (KEY == "raw") TARGET else "TARGET"
if (!all(lin_ids %in% tree$tip.label)) stopf("%d of the 76 lineage names (%s) are not tips of %s (first missing: %s). Is the tree the %s one? (set KEY or TREE)",
      sum(!lin_ids %in% tree$tip.label), ifelse(KEY == "raw", "raw MAG names", "public ids"), basename(TREE_F),
      lin_ids[!lin_ids %in% tree$tip.label][1], ifelse(KEY == "raw", "rawIDs", "publicID"))
if (!tgt_tip %in% tree$tip.label) stopf("target tip '%s' not found in the tree", tgt_tip)
ok("all 76 lineage names and the target are tips of the tree (%d tips in total)", length(tree$tip.label))

tips <- tree$tip.label
role <- rep(NA_character_, length(tips))
role[tips %in% lin_ids] <- "lineage"
role[tips == tgt_tip]  <- "target"
rest <- is.na(role)
role[rest] <- R$role[match(tips[rest], R$id)]
if (anyNA(role)) stopf("%d tips have no role (first: %s)", sum(is.na(role)), tips[is.na(role)][1])
role[role == "ref_neighbour"] <- "reference"
ok("roles: %s", paste(names(table(role)), table(role), collapse = "; "))

grp <- role
if (nzchar(GTDB_TAX) && file.exists(GTDB_TAX)) {
  tx <- read_tsv(GTDB_TAX, col_names = c("acc", "tax"), col_types = "cc", progress = FALSE)
  ord <- str_match(tx$tax[match(tips, tx$acc)], "o__([^;]*)")[, 2]
  ref <- role == "reference"
  grp[ref] <- ifelse(is.na(ord[ref]) | ord[ref] == "", "reference", paste0("GTDB ", ord[ref]))
  ok("reference orders: %s", paste(names(table(grp[ref])), table(grp[ref]), collapse = "; "))
} else cat("  (GTDB_TAX not set: references are shown as one 'reference' class)\n")
grp[grp == "lineage"] <- "RUC lineage (catalogue)"
grp[grp == "target"]  <- "Unclassified sheep genome"
grp[grp == "outgroup"] <- "Outgroup"

pal <- c("RUC lineage (catalogue)" = "#1F78B4", "Unclassified sheep genome" = "#E31A1C",
         "GTDB UMGS1840" = "#B39DDB", "GTDB SK-Y3" = "#8C510A", "GTDB RGIG8785" = "#80B1D3",
         "reference" = "#999999", "Outgroup" = "#000000")
lv <- names(pal)[names(pal) %in% grp]
extra <- setdiff(unique(grp), names(pal))
if (length(extra)) { pal <- c(pal, setNames(rep("#BDBDBD", length(extra)), extra)); lv <- c(lv, extra) }
f_role <- file.path(OUT_DIR, paste0("09_tip_role_colorstrip_", TAG, ".txt"))
con <- file(f_role, "w")
writeLines(c("DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tTip role", "COLOR\t#000000",
             "COLOR_BRANCHES\t0", "LEGEND_TITLE\tTip role",
             paste0("LEGEND_SHAPES\t", paste(rep(1, length(lv)), collapse = "\t")),
             paste0("LEGEND_COLORS\t", paste(pal[lv], collapse = "\t")),
             paste0("LEGEND_LABELS\t", paste(lv, collapse = "\t")), "DATA"), con)
writeLines(paste(tips, pal[grp], grp, sep = "\t"), con)
close(con)
ok("wrote %s (%d rows)", basename(f_role), length(tips))

if (FULL == "1") {
  cat("== C7 full-tree layout files ==\n")
  nt <- Ntip(tree)
  one_tip <- function(n) { while (n > nt) n <- tree$edge[tree$edge[, 1] == n, 2][1]; tips[n] }
  lca_id <- function(n) {
    kids <- tree$edge[tree$edge[, 1] == n, 2]
    paste(one_tip(kids[1]), one_tip(kids[length(kids)]), sep = "|")
  }
  supv <- suppressWarnings(as.numeric(tree$node.label))
  if (is.null(tree$node.label) || all(is.na(supv))) stopf("tree has no numeric support values on its nodes")

  coll <- c("GTDB SK-Y3", "GTDB RGIG8785", "Outgroup")
  lay <- data.frame(group = c("RUC + GTDB UMGS1840", coll), n_tips = NA, monophyletic = NA, clade_UFBoot = NA, node_id = NA)
  members <- list("RUC + GTDB UMGS1840" = tips[grp %in% c("RUC lineage (catalogue)", "GTDB UMGS1840")])
  for (k in coll) members[[k]] <- tips[grp == k]
  for (i in seq_len(nrow(lay))) {
    m <- members[[lay$group[i]]]; lay$n_tips[i] <- length(m)
    if (length(m) >= 2) {
      mono <- is.monophyletic(tree, m); n <- getMRCA(tree, m)
      lay$monophyletic[i] <- mono; lay$clade_UFBoot[i] <- supv[n - nt]
      if (mono) lay$node_id[i] <- lca_id(n)
    } else if (length(m) == 1) { lay$monophyletic[i] <- TRUE; lay$node_id[i] <- m }
  }
  write.csv(lay, file.path(OUT_DIR, "full_tree_layout_check.csv"), row.names = FALSE)
  print(lay[, 1:4], row.names = FALSE)
  if (!nzchar(GTDB_TAX)) cat("  (GTDB_TAX not set: SK-Y3/RGIG8785 cannot be identified, so nothing is collapsed)\n")

  cl <- lay[lay$group %in% coll & lay$n_tips >= 2 & !is.na(lay$node_id), ]
  bad <- lay$group[lay$group %in% coll & lay$n_tips >= 2 & lay$monophyletic %in% FALSE]
  if (length(bad)) cat("  WARN: not collapsed because not monophyletic: ", paste(bad, collapse = ", "), "\n", sep = "")
  f10 <- file.path(OUT_DIR, paste0("10_collapse_groups_", TAG, ".txt"))
  writeLines(c("COLLAPSE", "DATA", cl$node_id), f10)
  ok("wrote %s (%d clades collapsed: %s)", basename(f10), nrow(cl), paste(cl$group, collapse = ", "))
  if (nrow(cl)) {
    cat(paste(cl$node_id, pal[cl$group], cl$group, sep = "\t"), file = f_role, sep = "\n", append = TRUE)
    ok("added %d collapsed-clade rows to %s", nrow(cl), basename(f_role))
  }

  f11 <- file.path(OUT_DIR, paste0("11_target_star_", TAG, ".txt"))
  writeLines(c("DATASET_SYMBOL", "SEPARATOR TAB", "DATASET_LABEL\tUnclassified sheep genome", "COLOR\t#E31A1C",
               "MAXIMUM_SIZE\t14", "DATA", paste(tgt_tip, 3, 14, "#E31A1C", 1, 1, "", sep = "\t")), f11)
  ok("wrote %s", basename(f11))

  nodes <- (nt + 2):(nt + tree$Nnode)
  keep <- nodes[!is.na(supv[nodes - nt]) & supv[nodes - nt] >= SUP_MIN]
  f12 <- file.path(OUT_DIR, paste0("12_support95_dots_", TAG, ".txt"))
  writeLines(c("DATASET_SYMBOL", "SEPARATOR TAB", paste0("DATASET_LABEL\tUFBoot >= ", SUP_MIN), "COLOR\t#000000",
               "MAXIMUM_SIZE\t4", "LEGEND_TITLE\tNode support", "LEGEND_SHAPES\t2", "LEGEND_COLORS\t#000000",
               paste0("LEGEND_LABELS\tUFBoot >= ", SUP_MIN), "DATA",
               vapply(keep, function(n) paste(lca_id(n), 2, 4, "#000000", 1, 1, "", sep = "\t"), "")), f12)
  ok("wrote %s (%d of %d internal nodes have UFBoot >= %g)", basename(f12), length(keep), length(nodes), SUP_MIN)

  lab <- tips
  ii <- match(tips, L$genome); lab[!is.na(ii)] <- L$public_id[ii[!is.na(ii)]]
  lab[tips == tgt_tip] <- "Unclassified sheep genome"
  isref <- grepl("^GTDB ", grp); lab[isref] <- sub("^(RS_|GB_)", "", tips[isref])
  iso <- grp == "Outgroup"; host <- str_split_fixed(tips[iso], "__", 3)[, 2]
  lab[iso] <- paste0("Outgroup MAG", ifelse(nzchar(host), paste0(" (", host, ")"), ""))
  extra_lab <- character(0)
  for (i in seq_len(nrow(cl))) extra_lab <- c(extra_lab, paste(cl$node_id[i],
      sprintf("%s (n = %d)", ifelse(cl$group[i] == "Outgroup", "Outgroup", cl$group[i]), cl$n_tips[i]), sep = "\t"))
  if (env("LABEL_OTHERS", "0") != "1") {
    blank <- tips == tgt_tip | grp %in% c("Outgroup", "GTDB SK-Y3", "GTDB RGIG8785")
    lab[blank] <- " "
    if (length(extra_lab)) extra_lab <- paste(sub("\t.*$", "", extra_lab), " ", sep = "\t")
    ok("labels hidden for %d tips (sheep genome, outgroup, SK-Y3, RGIG8785) and %d collapsed triangles", sum(blank), length(extra_lab))
  }
  f13 <- file.path(OUT_DIR, paste0("13_tip_labels_", TAG, ".txt"))
  writeLines(c("LABELS", "SEPARATOR TAB", "DATA", paste(tips, lab, sep = "\t"), extra_lab), f13)
  ok("wrote %s (%d tip labels + %d collapsed-clade labels)", basename(f13), length(tips), length(extra_lab))

  fam <- rep(NA_character_, length(tips))
  if ("classification" %in% names(L)) {
    ii <- match(tips, if (KEY == "raw") L$genome else L$public_id)
    fam[!is.na(ii)] <- str_match(L$classification[ii[!is.na(ii)]], "f__([^;]*)")[, 2]
  } else cat("  WARN: label table has no 'classification' column; RUC families left blank\n")
  if (nzchar(GTDB_TAX) && exists("tx")) {
    ru <- grp == "GTDB UMGS1840"
    fam[ru] <- str_match(tx$tax[match(tips[ru], tx$acc)], "f__([^;]*)")[, 2]
  }
  fam[!is.na(fam) & fam == ""] <- NA
  src <- ifelse(grepl("^GTDB ", grp), "GTDB reference", "catalog")
  cat("  family x source:\n"); print(table(family = fam, source = src, useNA = "no"))
  ruc_tab <- table(fam[grp == "RUC lineage (catalogue)"])
  if (!identical(as.integer(ruc_tab[c("HGM11507", "UMGS1840")]), c(44L, 32L)))
    warning("RUC family counts differ from HGM11507 44 / UMGS1840 32")
  fams <- sort(unique(na.omit(fam)))
  fcol <- c(HGM11507 = "#A6CEE3", UMGS1840 = "#E31A1C")
  extra_c <- c("#FFD92F", "#E78AC3", "#999933", "#FB9A99", "#CAB2D6", "#BC80BD", "#FFFFB3", "#D9D9D9")
  oth <- setdiff(fams, names(fcol)); fcol <- c(fcol, setNames(extra_c[seq_along(oth)], oth)); fcol <- fcol[fams]
  f15 <- file.path(OUT_DIR, paste0("15_family_colorstrip_", TAG, ".txt"))
  hasf <- !is.na(fam)
  writeLines(c("DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tFamily", "COLOR\t#000000", "COLOR_BRANCHES\t0",
               "LEGEND_TITLE\tFamily", paste0("LEGEND_SHAPES\t", paste(rep(1, length(fams)), collapse = "\t")),
               paste0("LEGEND_COLORS\t", paste(fcol, collapse = "\t")), paste0("LEGEND_LABELS\t", paste(fams, collapse = "\t")),
               "DATA", paste(tips[hasf], fcol[fam[hasf]], fam[hasf], sep = "\t")), f15)
  ok("wrote %s (%d tips: %d catalog + %d GTDB references; families: %s)", basename(f15), sum(hasf),
     sum(hasf & src == "catalog"), sum(hasf & src == "GTDB reference"), paste(fams, collapse = ", "))

  lcol <- ifelse(tips == tgt_tip, "#E31A1C", ifelse(grepl("^GTDB ", grp), "#8C8C8C", ifelse(grp == "Outgroup", "#555555", "#000000")))
  lsty <- ifelse(tips == tgt_tip, "bold", "normal")
  rows14 <- paste(tips, "label", lcol, lsty, 1, sep = "\t")
  if (nrow(cl)) rows14 <- c(rows14, paste(cl$node_id, "label", ifelse(cl$group == "Outgroup", "#555555", "#8C8C8C"), "normal", 1, sep = "\t"))
  f14 <- file.path(OUT_DIR, paste0("14_label_colors_", TAG, ".txt"))
  writeLines(c("TREE_COLORS", "SEPARATOR TAB", "DATA", rows14), f14)
  ok("wrote %s (black %d, grey %d, red %d, dark grey %d)", basename(f14), sum(lcol == "#000000"),
     sum(lcol == "#8C8C8C"), sum(lcol == "#E31A1C"), sum(lcol == "#555555"))
}

if (PRUNE76 == "1") {
  sub <- keep.tip(tree, lin_ids)
  f_sub <- paste0("ml_tree_76tips_", ifelse(KEY == "raw", "rawIDs", "publicID"), ".nwk")
  write.tree(sub, file.path(OUT_DIR, f_sub))
  ok("wrote %s (%d tips; support kept as node labels)", f_sub, length(sub$tip.label))
}

if (KEY == "raw") {
  cat("== C6 convert the earlier URGG-keyed strips to raw MAG names ==\n")
  map <- setNames(L$genome, L$public_id)
  old_files <- c("01_host_species_branch_colors_PUBLIC_ID.txt", "02_genus_colorstrip_PUBLIC_ID.txt",
                 "02_genus_colorstrip_FAMILY_LEVEL_PUBLIC_ID.txt", "05_host_species_colorstrip_PUBLIC_ID.txt",
                 "06_novelty_tier_colorstrip_PUBLIC_ID.txt", "07_seqcode_eligible_colorstrip_PUBLIC_ID.txt",
                 "08_animal_colorstrip_PUBLIC_ID.txt")
  for (f in old_files) {
    src <- file.path(OLD_DIR, f)
    if (!file.exists(src)) { cat("  (not found, skipped: ", src, ")\n", sep = ""); next }
    ln <- readLines(src, warn = FALSE); d <- which(ln == "DATA")
    if (length(d) != 1) { cat("  WARN: no single DATA line in ", f, "; skipped\n", sep = ""); next }
    body <- ln[(d + 1):length(ln)]; body <- body[nzchar(body)]
    id <- sub("\t.*$", "", body); hit <- id %in% names(map)
    body[hit] <- paste0(map[id[hit]], substring(body[hit], nchar(id[hit]) + 1))
    if (sum(hit) != 76) warning(sprintf("%s: %d of %d data rows matched a lineage public_id (expected 76)", f, sum(hit), length(id)))
    body <- body[hit]
    out <- file.path(OUT_DIR, sub("PUBLIC_ID", "RAW_ID", f))
    writeLines(c(ln[1:d], body), out)
    bad <- setdiff(sub("\t.*$", "", body), tree$tip.label)
    ok("%s -> %s (%d rows; %d not tips of the tree)", f, basename(out), length(body), length(bad))
  }
}

cat("\nDONE (KEY=", KEY, "). Upload to iTOL from ", OUT_DIR, ": the tree (", basename(TREE_F), " or the 76-tip file), then\n",
    "  09_tip_role_colorstrip_", TAG, ".txt, 03b_cazyme_binary_RERUN_", TAG, ".txt and the converted\n",
    "  01/02/05/06/07/08 files (", ifelse(KEY == "raw", "..._RAW_ID.txt", "the existing ..._PUBLIC_ID.txt files"), ").\n",
    "  Full-tree layout (FULL=1): 10 (collapse), 11 (star), 12 (UFBoot dots), 13 (labels), 14 (label colours), 15 (family incl. references; use instead of 02 FAMILY_LEVEL).\n",
    sep = "")
