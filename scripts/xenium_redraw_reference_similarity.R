#!/usr/bin/env Rscript
# Redraw completed, current comparison CSVs without loading expression objects.
# Expected runtime <1-2 minutes; low confidence from CSV/PDF sizes, no HPC RDS IO.
args <- commandArgs(trailingOnly = TRUE)
legacy_flag <- "--legacy-query-labels" %in% args
allowed <- c("input", "out")
opt <- list()
for (a in setdiff(args, c("--dry-run", "--legacy-query-labels"))) {
  if (!grepl("^--[^=]+=.+$", a)) stop("Use --input=DIR --out=NEW_DIR [--dry-run].")
  k <- sub("=.*$", "", substring(a, 3L))
  if (!k %in% allowed || !is.null(opt[[k]])) stop("Unknown/repeated option: ", a)
  opt[[k]] <- sub("^[^=]+=", "", a)
}
if (!all(allowed %in% names(opt))) stop("Supply --input=COMPLETED_RESULTS --out=NEW_PLOT_DIR.")
if (!requireNamespace("here", quietly = TRUE)) stop("Missing here in project renv.")
here::i_am("scripts/xenium_redraw_reference_similarity.R")
source(here::here("scripts", "R", "reference_similarity_plots.R"))
resolve <- function(p) if (grepl("^/", p)) p else here::here(p)
input <- resolve(opt$input); out <- resolve(opt$out)
files <- file.path(input, paste0(rep(c("aldinger", "sepp", "science"), each = 2), "_",
                                rep(c("corrected_broad_label", "comb_subcluster"), 3), "_all_correlations.csv"))
if (file.exists(out)) stop("Output exists; choose a new plot directory: ", out)
if (!file.exists(file.path(input, "COMPLETED.txt")) || any(!file.exists(files))) stop("Input must contain a complete all-reference comparison and all six correlation CSVs.")
mode_file <- file.path(input, "comparison_mode.csv")
legacy_query_labels <- FALSE
if (file.exists(mode_file)) {
  mode <- read.csv(mode_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(mode) != 1L || !all(c("comparison_mode", "legacy_query_labels") %in% names(mode))) stop("Invalid comparison_mode.csv.")
  legacy_query_labels <- isTRUE(mode$legacy_query_labels[[1L]])
}
if (legacy_flag && !legacy_query_labels) legacy_query_labels <- TRUE
cat("Redraw:", length(files), "CSV inputs; no query/reference objects loaded. Query-label mode:", if (legacy_query_labels) "legacy" else "current", "\n")
if ("--dry-run" %in% args) quit(save = "no")
# Validate every matrix and cohort before writing any plot.
tables <- lapply(files, function(file) {
  d <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("query_label", "reference_label", "spearman", "specificity_spearman", "reference_annotation_status")
  if (!all(needed %in% names(d)) || anyNA(d[, needed]) || any(!is.finite(d$spearman)) || any(!is.finite(d$specificity_spearman))) stop("Incomplete correlation CSV: ", file)
  if (any(grepl("^Unannotated ", d$reference_label)) || any(d$reference_annotation_status == "unannotated_original_cluster")) stop("Old Sepp results include unannotated groups. Rerun the comparison to recompute centering/rankings after their exclusion.")
  if (grepl("^science_", basename(file)) && any(d$reference_annotation_status != "published_cell_class")) stop("Old Science numeric-cluster scores cannot be merged into broad-class scores. Rerun with original CellClass grouping first.")
  level <- if (grepl("_comb_subcluster_", basename(file), fixed = TRUE)) "comb_subcluster" else "corrected_broad_label"
  ql <- reference_similarity_order(unique(d$query_label), level, legacy_query_labels = legacy_query_labels); rl <- unique(d$reference_label)
  idx <- cbind(match(d$query_label, ql), match(d$reference_label, rl))
  if (nrow(d) != length(ql) * length(rl) || anyDuplicated(data.frame(idx))) stop("Correlation grid must be complete and unique: ", file)
  raw <- specific <- matrix(NA_real_, length(ql), length(rl), dimnames = list(ql, rl))
  raw[idx] <- d$spearman; specific[idx] <- d$specificity_spearman
  list(raw = raw, specific = specific, level = level)
})
if (!dir.create(out, recursive = TRUE)) stop("Cannot create plot output.")
heatmap_dir <- file.path(out, "heatmaps")
if (!dir.create(heatmap_dir, recursive = TRUE)) stop("Cannot create heatmap output directory.")
for (i in seq_along(files)) {
  stem <- sub("_all_correlations.csv$", "", basename(files[i]))
  plot_order <- write_reference_heatmaps(tables[[i]]$raw, tables[[i]]$specific,
                           file.path(heatmap_dir, paste0(stem, "_heatmaps.pdf")),
                           paste(stem, if (legacy_query_labels) "[LEGACY QUERY LABELS]" else ""),
                           query_level = tables[[i]]$level, legacy_query_labels = legacy_query_labels)
  write.csv(data.frame(top_to_bottom = seq_along(plot_order), query_label = plot_order),
            file.path(out, paste0(stem, "_heatmap_row_order.csv")), row.names = FALSE)
}
write.csv(data.frame(path = files, md5 = unname(tools::md5sum(files))), file.path(out, "plot_input_provenance.csv"), row.names = FALSE)
write.csv(data.frame(legacy_query_labels = legacy_query_labels), file.path(out, "comparison_mode.csv"), row.names = FALSE)
file.copy(here::here("scripts", "R", "reference_similarity_plots.R"), file.path(out, "plot_helpers_snapshot.R"))
file.copy(here::here("scripts", "color_palette.R"), file.path(out, "color_palette_snapshot.R"))
cat("Redrawn:", out, "\n")
