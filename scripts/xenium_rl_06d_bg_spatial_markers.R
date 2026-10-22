#!/usr/bin/env Rscript
# Read-only per-sample BG identity diagnostics, after RL Stage 7 mapping.
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "rl_bg_spatial_diagnostics.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)$sample_id
if (!length(samples) || anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples) ||
    any(!grepl("^[A-Za-z0-9_]+$", samples))) stop("Invalid or duplicate manifest sample IDs.")
root <- here(config$project$outputs_dir, "xenium", "rl")
task_map <- data.frame(task_id = seq_along(samples), sample_id = samples,
  input_path = file.path(root, "07_mapped", "rds", paste0(samples, "_RL_QC_Subclusters.rds")))
output_root <- file.path(root, "06_subclusters", "diagnostics", "bg_spatial_markers")
paths_for <- function(id) setNames(file.path(output_root, id,
  paste0(id, "_BG_", rl_bg_suffixes)), rl_bg_suffixes)
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
options_allowed <- c("--dry-run", "--overwrite", "--summarize")
unknown <- args[startsWith(args, "--") & !args %in% options_allowed]
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
summarize <- "--summarize" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
task_args <- args[!args %in% options_allowed]
if (length(task_args) > 1L || (summarize && length(task_args))) stop("Supply one TASK_ID or --summarize.")
if (summarize) {
  suffixes <- c("Summary.csv", "Marker_Availability.csv", "Marker_Detection.csv", "Coexpression.csv")
  expected <- unlist(lapply(samples, paths_for), use.names = FALSE)
  combined_paths <- file.path(output_root, paste0("AllSamples_BG_", suffixes))
  if (dry_run) {
    ok <- compact_dry_run("BG summary: all manifest samples required", expected, combined_paths)
    quit(save = "no", status = if (ok) 0L else 1L)
  }
  if (any(!file.exists(expected))) stop("Partial sample output set; run --summarize --dry-run for missing inputs.")
  if (any(file.exists(combined_paths)) && !overwrite) stop("Combined tables already exist; explicit --overwrite required.")
  # Check completion records and input/config/manifest provenance before combining.
  for (i in seq_along(samples)) {
    p <- read.csv(paths_for(samples[i])[["Provenance.csv"]], stringsAsFactors = FALSE)
    info <- file.info(task_map$input_path[i])
    if (nrow(p) != 1L || !identical(p$sample_id, samples[i]) ||
        !identical(as.integer(p$task_id), as.integer(i)) ||
        !file.exists(task_map$input_path[i]) ||
        !identical(p$input_path, normalizePath(task_map$input_path[i])) ||
        !isTRUE(all.equal(p$input_size, as.numeric(info$size))) ||
        !isTRUE(all.equal(p$input_mtime, as.numeric(info$mtime))) ||
        !identical(p$config_md5, unname(tools::md5sum(attr(config, "config_path")))) ||
        !identical(p$manifest_md5, unname(tools::md5sum(resolve_config_path(config$manifests$samples, config)))) ||
        !identical(p$helper_md5, unname(tools::md5sum(here("scripts", "R", "rl_bg_spatial_diagnostics.R"))))) {
      stop("Stale or invalid completion provenance: ", samples[i])
    }
  }
  tables <- lapply(suffixes, function(suffix) {
    pieces <- lapply(seq_along(samples), function(i) {
      z <- read.csv(paths_for(samples[i])[[suffix]], stringsAsFactors = FALSE)
      if (!nrow(z) || anyNA(z$sample_id) || !all(z$sample_id == samples[i]) ||
          anyNA(z$task_id) || !all(z$task_id == i)) stop("Invalid sample/task table: ", samples[i], " / ", suffix)
      z
    })
    do.call(rbind, pieces)
  })
  if (nrow(tables[[1]]) != length(samples) || anyDuplicated(tables[[1]]$sample_id) ||
      sum(tables[[1]]$BG_cells) == 0L) {
    stop("Invalid cohort coverage or no BG Progenitors cells in the entire cohort; verify the current labels.")
  }
  for (i in seq_along(tables)) write.csv(tables[[i]], combined_paths[i], row.names = FALSE)
  message("Combined ", length(samples), " samples; sample-level rows retained without pooling donor inference.")
  quit(save = "no", status = 0L)
}
task_value <- if (length(task_args)) task_args else Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (!nzchar(task_value) && dry_run) {
  ok <- compact_dry_run("All BG spatial/marker tasks (path-only)", task_map$input_path,
    unlist(lapply(samples, paths_for), use.names = FALSE))
  cat("Samples: ", length(samples), "; no RDS objects loaded. Use --list for task mapping.\n", sep = "")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!grepl("^[1-9][0-9]*$", task_value)) stop("Supply a positive integer TASK_ID; use --list.")
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id > nrow(task_map)) stop("TASK_ID exceeds manifest size: ", nrow(task_map))
sample_id <- samples[task_id]
input <- task_map$input_path[task_id]
output_paths <- paths_for(sample_id)
if (dry_run) {
  ok <- compact_dry_run(paste("BG diagnostic task", task_id, sample_id), input, output_paths)
  cat("Path-only; no RDS objects loaded.\n")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!file.exists(input)) stop("Mapped RL object missing: ", input)
if (any(file.exists(output_paths)) && !overwrite) {
  stop("Existing/partial diagnostic outputs protected: ", sample_id,
       ". Review before using --overwrite.")
}
# An old completion record must not survive an interrupted explicit overwrite.
if (overwrite && file.exists(output_paths[["Provenance.csv"]])) {
  if (!file.remove(output_paths[["Provenance.csv"]])) stop("Cannot invalidate old completion record.")
}
rl_bg_analyze(input, sample_id, task_id, output_paths,
  attr(config, "config_path"), resolve_config_path(config$manifests$samples, config))
