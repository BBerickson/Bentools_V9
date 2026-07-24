#!/usr/bin/env Rscript
# Data-layer regression harness for BenTools v9.
#
# Purpose: exercise the perf-critical, deterministic data-layer functions in
# R_scripts/functions.R against the fixtures in test_files/ and emit a stable
# "fingerprint" (row counts, distinct-gene counts, rounded score sums, cluster
# sizes, ...). Run it BEFORE a refactor step to capture a golden fingerprint,
# then AFTER to confirm the numbers are unchanged.
#
# Usage:
#   Rscript tests/regression_harness.R            # print fingerprint
#   Rscript tests/regression_harness.R --save     # also write golden RDS
#   Rscript tests/regression_harness.R --check     # compare against golden RDS
#
# This covers the data layer only. The Shiny reactivity / UI is verified via the
# manual smoke path documented in tests/SMOKE_TEST.md.

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(valr)
  library(RColorBrewer)
  library(colorspace)
  library(fastcluster)
}))

# --- run from repo root regardless of cwd ---
args <- commandArgs(trailingOnly = TRUE)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- if (length(this_file)) normalizePath(file.path(dirname(this_file), "..")) else getwd()
setwd(repo_root)

# Shiny UI stubs so functions that pop modals on edge cases don't error headless.
showModal <- function(...) invisible(NULL)
modalDialog <- function(...) invisible(NULL)
removeModal <- function(...) invisible(NULL)
withProgress <- function(expr, ...) force(expr)
setProgress <- function(...) invisible(NULL)
incProgress <- function(...) invisible(NULL)

kBrewerList <- c("Set1", "Paired", "Dark2", "Spectral")
set.seed(1)  # PrepMetaFile assigns random colors; pin them

source("R_scripts/functions.R", local = TRUE)

MATRIX_FILES <- c("A_test.matrix.gz", "B_test.matrix.gz",
                  "C_test.matrix.gz", "D_test.matrix.gz")
GENELIST_FILE <- "test_genelist.txt"

# ---- Build a stable list_data the way the loader does (bind once) ----
build_list_data <- function() {
  ld_list <- list()
  binning <- NULL
  for (f in MATRIX_FILES) {
    path <- file.path("test_files", f)
    meta <- PrepMetaFile(path, f)
    bc <- tableTestbin(meta[1, ])
    if (is.null(binning)) binning <- bc$binning
    ld_list[[f]] <- LoadTableFile(meta[1, ], bc)
  }
  table_file <- bind_rows(ld_list) %>% filter(!is.na(set))
  sets <- unique(table_file$set)

  complete_full <- distinct(table_file, gene, chrom, start, end, strand)
  gene_file <- list(Complete = list(full = complete_full,
                                     info = tibble(loaded_info = "all loaded genes")))
  # meta_data: one row per (Complete, sample); onoff = set => all active
  meta_data <- tibble(
    gene_list = "Complete",
    count = paste("n =", n_distinct(complete_full$gene)),
    set = sets,
    mycol = "#1F78B4",
    onoff = sets,
    sub = " ",
    plot_legend = " ",
    group = "self"
  )
  list(
    table_file = table_file,
    gene_file = gene_file,
    meta_data = meta_data,
    meta_data_plot = list(binning = binning, rnaseq = FALSE)
  )
}

# round helper for stable float fingerprints
rs <- function(x) round(sum(x, na.rm = TRUE), 3)

fingerprint <- list()
ld <- build_list_data()

fingerprint$load <- list(
  n_sets = n_distinct(ld$table_file$set),
  n_rows = nrow(ld$table_file),
  n_genes = n_distinct(ld$table_file$gene),
  n_bins = n_distinct(ld$table_file$bin),
  score_sum = rs(ld$table_file$score),
  complete_genes = n_distinct(ld$gene_file$Complete$full$gene)
)
# column types matter downstream (e.g. Active_list_data's grouped mutates); guard them
fingerprint$load_types <- as.list(vapply(
  ld$table_file[c("chrom","start","end","gene","value","strand","bin","score","set")],
  function(x) class(x)[1], character(1)))

# ---- Load a real gene list as a second gene_file entry ----
ld2 <- tryCatch(
  LoadGeneFile(file.path("test_files", GENELIST_FILE), GENELIST_FILE, ld),
  error = function(e) { message("LoadGeneFile: ", conditionMessage(e)); ld }
)
if (!is.null(ld2) && length(ld2$gene_file) > length(ld$gene_file)) ld <- ld2
gl_name <- setdiff(names(ld$gene_file), "Complete")
fingerprint$genelists <- names(ld$gene_file)

# ---- FilterSepSize (pure) ----
fss <- FilterSepSize(distinct(ld$table_file, gene, chrom, start, end, strand),
                     separation = 500, minsize = 1000, maxsize = 0, stranded = FALSE)
fingerprint$filter_sep_size <- list(n_genes = n_distinct(fss$gene))

# ---- Active_list_data + ApplyMath (per-refresh hot path) ----
active <- Active_list_data(ld, group = FALSE, fulljoin = FALSE)
fingerprint$active <- list(n_rows = nrow(active), n_genes = n_distinct(active$gene),
                           score_sum = rs(active$score))
am <- ApplyMath(active, use_math = "mean")   # output score column is `value`
fingerprint$apply_math_mean <- list(n_rows = nrow(am), value_sum = rs(am$value))

# ---- Clustering ----
sample1 <- unique(ld$table_file$set)[1]
bin_hi <- max(ld$table_file$bin)
fc <- tryCatch(FindClusters(ld, "Complete", sample1, c(1, bin_hi), "pattern"),
               error = function(e) { message("FindClusters: ", conditionMessage(e)); NULL })
if (!is.null(fc$clust)) {
  # $clust$cm is an hclust; $clust$full holds the gene coords used
  merges <- if (inherits(fc$clust$cm, "hclust")) nrow(fc$clust$cm$merge) else NA_integer_
  fingerprint$find_clusters <- list(n_clust_genes = n_distinct(fc$clust$full$gene),
                                     hclust_merges = merges)
}

# ---- Groups ----
fg <- tryCatch(FindGroups(ld, "Complete", sample1, c(1, bin_hi)),
               error = function(e) { message("FindGroups: ", conditionMessage(e)); NULL })
if (!is.null(fg$groupies)) {
  fingerprint$find_groups <- list(n_group_genes = n_distinct(fg$groupies$full$gene))
} else if (!is.null(fg)) {
  fingerprint$find_groups <- "ran"
}

# stable per-column numeric aggregate (sorted names -> order-independent)
nsum <- function(df) {
  num <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
  if (!ncol(num)) return(list())
  as.list(round(vapply(num[order(names(num))], function(x) sum(x, na.rm = TRUE),
                       numeric(1)), 3))
}

# ---- FilterTop (lapply/inner_join accumulation; Phase 6 target) ----
ft <- tryCatch(FilterTop(ld, "Complete", sample1, c(1, 40), "1:40", 50, "top%"),
               error = function(e) { message("FilterTop: ", conditionMessage(e)); NULL })
if (!is.null(ft)) {
  new_list <- setdiff(names(ft$gene_file), names(ld$gene_file))
  fingerprint$filter_top <- list(
    new_list = new_list,
    n_genes = if (length(new_list)) n_distinct(ft$gene_file[[new_list[1]]]$full$gene) else 0L)
}
# "Middle%" (filter-between) exercises the count() path that must resolve to
# dplyr::count, not matrixStats::count
ftmid <- tryCatch(FilterTop(ld, "Complete", sample1, c(1, 40), "1:40", 50, "Middle%"),
                  error = function(e) { message("FilterTop(middle): ", conditionMessage(e)); NULL })
if (!is.null(ftmid)) {
  nlm <- setdiff(names(ftmid$gene_file), names(ld$gene_file))
  fingerprint$filter_top_middle <- list(
    n_genes = if (length(nlm)) n_distinct(ftmid$gene_file[[nlm[1]]]$full$gene) else 0L)
}
# multi-file FilterTop exercises the reduce(inner_join) path across samples
two_samples <- unique(ld$table_file$set)[1:2]
ftm <- tryCatch(FilterTop(ld, "Complete", two_samples, c(1, 40), "1:40", 50, "Top%"),
                error = function(e) { message("FilterTop(multi): ", conditionMessage(e)); NULL })
if (!is.null(ftm)) {
  nl <- setdiff(names(ftm$gene_file), names(ld$gene_file))
  fingerprint$filter_top_multi <- list(
    n_cols = if (length(nl)) ncol(ftm$gene_file[[nl[1]]]$full) else 0L,
    n_genes = if (length(nl)) n_distinct(ftm$gene_file[[nl[1]]]$full$gene) else 0L)
}

# ---- FilterPer (ratcheting while-loop; Phase 6 target) ----
fper <- tryCatch(FilterPer(ld, "Complete", sample1, c(1, 40), c(10, 90), "per%", "1:40"),
                 error = function(e) { message("FilterPer: ", conditionMessage(e)); NULL })
if (!is.null(fper$sortplot)) {
  fingerprint$filter_per <- c(list(n_rows = nrow(fper$sortplot)), nsum(fper$sortplot))
}

# ---- MakeGroupFile (grow-in-loop bind_rows; Phase 6 target) ----
mgf <- tryCatch(MakeGroupFile(ld, "mean"),
                error = function(e) { message("MakeGroupFile: ", conditionMessage(e)); NULL })
if (!is.null(mgf$table_file)) {
  new_sets <- setdiff(unique(mgf$table_file$set), unique(ld$table_file$set))
  fingerprint$make_group_file <- list(
    new_sets = new_sets,
    new_rows = sum(mgf$table_file$set %in% new_sets),
    new_score_sum = rs(mgf$table_file$score[mgf$table_file$set %in% new_sets]))
}

# ---- ApplyTtest / try_t_test (add_row bin loop; Phase 6 target) ----
att <- tryCatch(
  ApplyTtest(active, "by files", "-log10", "wilcox.test", "fdr", "two.sided", "FALSE", "FALSE"),
  error = function(e) { message("ApplyTtest: ", conditionMessage(e)); NULL })
if (!is.null(att) && is.data.frame(att)) {
  fingerprint$apply_ttest <- c(list(n_rows = nrow(att)), nsum(att))
}

# ---- Phase-5 decoupled functions (read list_data, not the global) ----
# sizes of the new gene lists a tool appends, keyed by list name
list_sizes <- function(res, pattern) {
  nm <- grep(pattern, names(res$gene_file), value = TRUE)
  setNames(as.integer(vapply(nm, function(n) n_distinct(res$gene_file[[n]]$full$gene),
                             integer(1))), sub("\n.*$", "", nm))[order(sub("\n.*$", "", nm))]
}

fa <- tryCatch(FilterAverage(ld, "Complete", sample1, c(1, 40), "1:40", "mean"),
               error = function(e) { message("FilterAverage: ", conditionMessage(e)); NULL })
if (!is.null(fa)) fingerprint$filter_average <- as.list(list_sizes(fa, "^Filter_all_bins"))

# multi-file FilterAverage exercises the reduce(full_join + merge) path
fam <- tryCatch(FilterAverage(ld, "Complete", unique(ld$table_file$set)[1:2], c(1, 40), "1:40", "mean"),
                error = function(e) { message("FilterAverage(multi): ", conditionMessage(e)); NULL })
if (!is.null(fam)) fingerprint$filter_average_multi <- as.list(list_sizes(fam, "^Filter_all_bins"))

if (!is.null(fc$clust)) {
  cnl <- tryCatch(ClusterNumList(fc, "Complete", sample1, "1:80", 3),
                  error = function(e) { message("ClusterNumList: ", conditionMessage(e)); NULL })
  if (!is.null(cnl)) fingerprint$cluster_num_list <- as.list(list_sizes(cnl, "^Cluster_"))
}

if (!is.null(fg$groupies)) {
  gnl <- tryCatch(GroupsNumList(fg, "Complete", sample1, "1:80", 3),
                  error = function(e) { message("GroupsNumList: ", conditionMessage(e)); NULL })
  if (!is.null(gnl)) fingerprint$groups_num_list <- as.list(list_sizes(gnl, "^Groups_"))
}

# ---- report ----
cat("\n===== BenTools data-layer fingerprint =====\n")
str(fingerprint, max.level = 3, digits.d = 6)

golden_path <- "tests/golden_fingerprint.rds"
if ("--save" %in% args) {
  saveRDS(fingerprint, golden_path)
  cat("\nSaved golden fingerprint ->", golden_path, "\n")
} else if ("--check" %in% args) {
  if (!file.exists(golden_path)) stop("No golden fingerprint; run with --save first.")
  golden <- readRDS(golden_path)
  ok <- identical(golden, fingerprint)
  if (ok) {
    cat("\nPASS: fingerprint matches golden.\n")
  } else {
    cat("\nFAIL: fingerprint differs from golden.\n")
    cat("Golden:\n"); str(golden, max.level = 3)
    quit(status = 1)
  }
}
