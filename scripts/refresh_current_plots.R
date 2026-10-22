#!/usr/bin/env Rscript
# Regenerate current production figure families from saved objects into a new
# snapshot. Never writes RDS, reruns integration, or edits source annotations.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/refresh_current_plots.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "palette_refresh.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
args <- commandArgs(trailingOnly = TRUE)
valid <- args %in% c("--list", "--dry-run", "--validate") |
  grepl("^--group=(references|samples|pooled|combined)$", args) |
  grepl("^--tag=[A-Za-z0-9][A-Za-z0-9_-]*$", args) | grepl("^[0-9]+$", args)
if (any(!valid)) stop("Use --group=references|samples|pooled|combined, --tag=NAME, [--list|--dry-run|--validate], [TASK_ID].")
option <- function(prefix, default) {
  value <- args[startsWith(args, paste0(prefix, "="))]
  if (length(value) > 1L) stop("Repeated option: ", prefix)
  if (length(value)) sub("^[^=]*=", "", value) else default
}
group <- option("--group", "references")
tag <- option("--tag", "palette_20261019")
if (sum(c("--list", "--dry-run", "--validate") %in% args) > 1L) stop("Choose one inspection mode.")
root <- here(config$project$outputs_dir)
snapshot <- file.path(root, "plot_refresh", tag, group)
rows <- list()
add <- function(task, id, stage, path, column, kind, spatial = FALSE) {
  rows[[length(rows) + 1L]] <<- data.frame(task_id = task, id = id, stage = stage,
    input_path = path, label_column = column, kind = kind, spatial = spatial,
    output_dir = file.path(snapshot, id, stage), stringsAsFactors = FALSE)
}
if (group == "references") {
  for (reference in names(config$inputs$references)) {
    panel <- resolve_config_path(config$inputs$references[[reference]], config)
    full <- sub("newUMAPv2_5k[.]rds$", "newUMAPv1.rds", panel)
    if (identical(full, panel)) stop("Cannot derive the standardized full-reference path: ", panel)
    add(length(rows) + 1L, reference, "full", full, "clusters_refined", "broad")
    add(length(rows) + 1L, reference, "panel", panel, "clusters_refined", "broad")
  }
} else if (group == "samples") {
  production <- get_production_output_paths(config)
  corrected_broad_paths <- refresh_corrected_broad_paths(root, samples$sample_id)
  for (i in seq_len(nrow(samples))) {
    id <- samples$sample_id[[i]]
    add(i, id, "consensus", file.path(production$consensus_rds_dir, paste0(id, "_Consensus_annotated.rds")), "consensus_label", "broad", TRUE)
    add(i, id, "vz", file.path(root, "xenium", "vz", "07_mapped", "rds", paste0(id, "_VZ_QC_Subclusters.rds")), "VZ_subcluster", "vz", TRUE)
    add(i, id, "rl", file.path(root, "xenium", "rl", "07_mapped", "rds", paste0(id, "_RL_QC_Subclusters.rds")), "RL_subcluster", "rl", TRUE)
    add(i, id, "residual", file.path(root, "xenium", "residual", "06_mapped", "rds", paste0(id, "_Residual_Subclusters.rds")), "residual_subcluster", "residual", TRUE)
    add(i, id, "combined", file.path(root, "xenium", "vz_rl", "01_combined_labels", "rds", paste0(id, "_VZ_RL_QC_Subclusters.rds")), "comb_subcluster", "combined", TRUE)
    add(i, id, "corrected_broad", corrected_broad_paths[[id]], "corrected_broad_label", "broad", TRUE)
  }
} else if (group == "pooled") {
  add(1L, "vz", "subclusters", file.path(root, "xenium", "vz", "06_subclusters", "rds", "Xenium_VZ_subclusters.rds"), "VZ_subcluster", "vz")
  add(2L, "rl", "subclusters", file.path(root, "xenium", "rl", "06_subclusters", "rds", "Xenium_RL_subclusters.rds"), "RL_subcluster", "rl")
  add(3L, "residual", "subclusters", file.path(root, "xenium", "residual", "05_subclusters", "rds", "Xenium_Residual_Subclusters.rds"), "residual_subcluster", "residual")
} else {
  for (label in c("consensus_label", "comb_subcluster")) {
    add(1L, "nonspatial", label, file.path(root, "xenium", "vz_rl", "nonspatial", "04_processed", "rds", "Xenium_VZRL_processed.rds"), label, if (label == "consensus_label") "broad" else "combined")
    add(2L, "spatial", label, file.path(root, "xenium", "vz_rl", "spatial", "04_integrated", "rds", "Xenium_VZRL_spatial_integrated.rds"), label, if (label == "consensus_label") "broad" else "combined")
  }
}
task_map <- do.call(rbind, rows)
if ("--list" %in% args) {
  utils::write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
task_args <- args[grepl("^[0-9]+$", args)]
if (length(task_args) > 1L) stop("Specify at most one TASK_ID.")
task_value <- if (length(task_args)) task_args else Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (nzchar(task_value)) {
  task <- suppressWarnings(as.integer(task_value))
  if (is.na(task) || !task %in% task_map$task_id) stop("TASK_ID is outside --list mapping.")
  task_map <- task_map[task_map$task_id == task, , drop = FALSE]
} else if (!any(c("--dry-run", "--validate") %in% args)) stop("Specify TASK_ID; inspect --list first.")
review_path <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels, config)
code_paths <- here("scripts", c("color_palette.R", "refresh_current_plots.R", "R/palette_refresh.R", "R/dotplot_helpers.R", "R/residual_spatial_reports.R", "R/spatial_report_helpers.R"))
manifests <- file.path(task_map$output_dir, "Refresh_manifest.csv")
if ("--dry-run" %in% args) {
  inputs <- unique(c(task_map$input_path, code_paths, if (any(task_map$kind == "residual")) review_path))
  ready <- all(file.exists(inputs))
  protected <- vapply(task_map$output_dir, function(path) dir.exists(path) &&
    length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L, logical(1))
  compact_dry_run(paste("Plot snapshot", group, tag), inputs = inputs, outputs = manifests,
                 checks = c(snapshot_directories_empty = !any(protected)))
  cat("  Tasks: ", length(unique(task_map$task_id)), "; object/label views: ", nrow(task_map),
      "; path-only: no RDS loads; existing snapshots are never overwritten.\n", sep = "")
  quit(save = "no", status = if (ready && !any(protected)) 0L else 1L)
}
palette_hash <- unname(tools::md5sum(code_paths[[1]]))
code_hash <- paste(unname(tools::md5sum(code_paths)), collapse = ":")
if ("--validate" %in% args) {
  for (i in seq_len(nrow(task_map))) {
    if (!file.exists(manifests[[i]])) stop("Missing completion manifest: ", manifests[[i]])
    manifest <- utils::read.csv(manifests[[i]], stringsAsFactors = FALSE)
    required <- c("output_path", "input_rds", "input_size", "input_mtime", "palette_md5", "code_md5", "review_md5", "expected_output_count")
    if (!nrow(manifest) || !all(required %in% names(manifest)) || anyDuplicated(manifest$output_path)) stop("Invalid refresh manifest.")
    info <- file.info(manifest$output_path)
    source_info <- file.info(task_map$input_path[[i]])
    expected_review <- if (task_map$kind[[i]] == "residual") unname(tools::md5sum(review_path)) else "not_used"
    if (anyNA(manifest$expected_output_count) || any(manifest$expected_output_count != nrow(manifest)) ||
        anyNA(info$size) || any(info$size <= 0) ||
        any(manifest$input_rds != normalizePath(task_map$input_path[[i]])) ||
        any(manifest$input_size != source_info$size) || any(manifest$input_mtime != as.numeric(source_info$mtime)) ||
        any(manifest$palette_md5 != palette_hash) || any(manifest$code_md5 != code_hash) ||
        any(manifest$review_md5 != expected_review)) stop("Stale or partial snapshot: ", task_map$output_dir[[i]])
  }
  cat("PASS | ", group, " | ", nrow(task_map), " completed object/label views | source provenance and output files current\n", sep = "")
  quit(save = "no", status = 0L)
}
missing <- unique(task_map$input_path[!file.exists(task_map$input_path)])
if (length(missing)) stop("Missing inputs:\n- ", paste(missing, collapse = "\n- "))
if (any(vapply(task_map$output_dir, function(path) dir.exists(path) &&
    length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L, logical(1)))) stop("Snapshot outputs exist. Choose a new --tag.")
dependencies <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "scales", "scattermore")
missing_dependencies <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_dependencies) || !capabilities("cairo")) stop("Missing plotting dependencies/cairo: ",
  paste(missing_dependencies, collapse = ", "), ". Resolve through the project renv workflow; no packages are installed here.")
source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))
source(here("scripts", "R", "residual_spatial_reports.R"))
source(here("scripts", "R", "palette_refresh.R"))
set.seed(config$runtime$random_seed)
review <- if (any(task_map$kind == "residual")) load_reviewed_residual_subcluster_labels(config) else NULL
# Null device prevents layout measurement from creating an incidental Rplots.pdf.
grDevices::pdf(file = NULL)
for (path in unique(task_map$input_path)) {
  source_info <- file.info(path)
  object <- readRDS(path)
  for (i in which(task_map$input_path == path)) {
    row <- task_map[i, , drop = FALSE]
    manifest <- refresh_plot_object(object, path, row$output_dir, row$label_column,
                                    row$kind, paste(row$id, row$stage), row$spatial, review)
    manifest$palette_md5 <- palette_hash
    manifest$code_md5 <- code_hash
    manifest$review_md5 <- if (row$kind == "residual") unname(tools::md5sum(review_path)) else "not_used"
    manifest$expected_output_count <- nrow(manifest)
    if (!identical(source_info[c("size", "mtime")], file.info(path)[c("size", "mtime")])) {
      stop("Source RDS changed while plotting: ", path)
    }
    partial_manifest <- paste0(manifests[[i]], ".partial")
    utils::write.csv(manifest, partial_manifest, row.names = FALSE)
    if (!file.rename(partial_manifest, manifests[[i]])) stop("Could not install completion manifest.")
    message("Completed: ", row$output_dir)
  }
  rm(object)
  gc()
}
grDevices::dev.off()
