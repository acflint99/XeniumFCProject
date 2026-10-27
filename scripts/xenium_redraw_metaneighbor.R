#!/usr/bin/env Rscript
# Redraw completed MetaNeighbor heatmaps from AUROC tables only; no expression objects are loaded.
# Expected runtime <1-2 minutes, depending on filesystem latency and PDF rendering.
args <- commandArgs(trailingOnly = TRUE)
legacy_flag <- "--legacy-query-labels" %in% args
allowed <- c("input", "out")
opt <- list()
for (a in setdiff(args, c("--dry-run", "--legacy-query-labels"))) {
  if (!grepl("^--[^=]+=.+$", a)) stop("Use --input=COMPLETED_REFERENCE_DIR --out=NEW_PLOT_DIR [--dry-run].")
  key <- sub("=.*$", "", substring(a, 3L))
  if (!key %in% allowed || !is.null(opt[[key]])) stop("Unknown/repeated option: ", a)
  opt[[key]] <- sub("^[^=]+=", "", a)
}
if (!all(allowed %in% names(opt))) stop("Supply --input=COMPLETED_REFERENCE_DIR --out=NEW_PLOT_DIR.")
if (!requireNamespace("here", quietly = TRUE)) stop("Missing here in project renv.")
here::i_am("scripts/xenium_redraw_metaneighbor.R")
source(here::here("scripts", "R", "reference_similarity_plots.R"))
source(here::here("scripts", "R", "metaneighbor_results.R"))
resolve <- function(p) if (grepl("^/", p)) p else here::here(p)
input <- resolve(opt$input)
out <- resolve(opt$out)
if (file.exists(out)) stop("Output exists; choose a new plot directory: ", out)
if (!file.exists(file.path(input, "COMPLETED.txt"))) stop("Input must be a completed single-reference MetaNeighbor result directory.")
mode_file <- file.path(input, "comparison_mode.csv")
legacy_query_labels <- FALSE
if (file.exists(mode_file)) {
  mode <- read.csv(mode_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(mode) != 1L || !"legacy_query_labels" %in% names(mode)) stop("Invalid comparison_mode.csv.")
  legacy_query_labels <- isTRUE(mode$legacy_query_labels[[1L]])
}
if (legacy_flag && !legacy_query_labels) legacy_query_labels <- TRUE
files <- list.files(input, pattern = "_(corrected_broad_label|comb_subcluster)_all_AUROCs\\.csv$", full.names = TRUE)
if (length(files) != 2L) stop("Input must contain exactly two current *_all_AUROCs.csv files.")
if ("--dry-run" %in% args) {
  cat("Redraw:", length(files), "AUROC CSV inputs; no query/reference objects loaded. Query-label mode:", if (legacy_query_labels) "legacy" else "current", "\n")
  quit(save = "no")
}

make_scores <- function(file) {
  d <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("query_label", "reference_label", "query_recovery_AUROC", "reference_recovery_AUROC", "mean_bidirectional_AUROC")
  if (!all(required %in% names(d)) || anyNA(d[, required]) || any(!is.finite(as.matrix(d[, required[3:5]])))) stop("Incomplete AUROC CSV: ", file)
  level <- if (grepl("_comb_subcluster_", basename(file), fixed = TRUE)) "comb_subcluster" else "corrected_broad_label"
  ql <- reference_similarity_order(unique(d$query_label), level, legacy_query_labels = legacy_query_labels)
  rl <- unique(d$reference_label)
  index <- cbind(match(d$query_label, ql), match(d$reference_label, rl))
  if (nrow(d) != length(ql) * length(rl) || anyNA(index) || anyDuplicated(data.frame(index))) stop("AUROC grid must be complete and unique: ", file)
  matrix_from <- function(column) {
    result <- matrix(NA_real_, length(ql), length(rl), dimnames = list(ql, rl))
    result[index] <- d[[column]]
    result
  }
  list(scores = list(forward = matrix_from("query_recovery_AUROC"),
                     reverse = matrix_from("reference_recovery_AUROC"),
                     mean = matrix_from("mean_bidirectional_AUROC")),
       level = level,
       reference = sub("_(corrected_broad_label|comb_subcluster)_all_AUROCs\\.csv$", "", basename(file)))
}

tables <- lapply(files, make_scores)
if (!dir.create(out, recursive = TRUE)) stop("Cannot create plot output directory.")
heatmap_dir <- file.path(out, "heatmaps")
if (!dir.create(heatmap_dir, recursive = TRUE)) stop("Cannot create heatmap output directory.")
for (table in tables) {
  title <- paste("MetaNeighborUS", table$reference, table$level, if (legacy_query_labels) "[LEGACY QUERY LABELS]" else "")
  mn_heatmaps(table$scores, file.path(heatmap_dir, paste0(table$reference, "_", table$level, "_AUROC_heatmaps.pdf")), title)
}
write.csv(data.frame(path = files, md5 = unname(tools::md5sum(files))), file.path(out, "plot_input_provenance.csv"), row.names = FALSE)
write.csv(data.frame(legacy_query_labels = legacy_query_labels), file.path(out, "comparison_mode.csv"), row.names = FALSE)
file.copy(here::here("scripts", "R", "metaneighbor_results.R"), file.path(out, "metaneighbor_plot_helpers_snapshot.R"))
file.copy(here::here("scripts", "R", "reference_similarity_plots.R"), file.path(out, "query_order_helpers_snapshot.R"))
file.copy(here::here("scripts", "color_palette.R"), file.path(out, "color_palette_snapshot.R"))
writeLines("Complete: MetaNeighbor heatmaps redrawn from completed AUROC tables only.", file.path(out, "COMPLETED.txt"))
cat("Redrawn:", out, "\n")
