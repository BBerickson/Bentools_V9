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
