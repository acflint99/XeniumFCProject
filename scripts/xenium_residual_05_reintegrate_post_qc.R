#!/usr/bin/env Rscript
# Expected 15-120 min after removals, 1-10 min for reuse (~10k cells; low
# confidence). Path-only --dry-run: <1 min, no RDS loads. HPC production only.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_05_reintegrate_post_qc.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
source(here("scripts", "R", "residual_post_qc_processing.R"))
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(!args %in% c("--dry-run", "--overwrite", "--no-removals")) ||
      all(c("--dry-run", "--overwrite") %in% args)) stop("Usage: [--no-removals] [--dry-run|--overwrite]")
  config <- load_pipeline_config()
  p <- residual_post_qc_paths(config)
  q <- residual_qc_removal_paths(config)
  no_removals <- "--no-removals" %in% args
  input <- if (no_removals) q$source else q$filtered
  merge_manifest <- here(config$project$outputs_dir, "xenium", "residual", "02_merged", "tables", "Xenium_Residual_Merge_Manifest.csv")
  inputs <- c(input, merge_manifest, if (!no_removals) unlist(q[c("decision", "source", "removed", "manifest", "provenance")]))
  resolutions <- c(0.5, 0.8, 1)
  plots <- file.path(dirname(dirname(p$input)), "plots", paste0("Xenium_Residual_PostQC_UMAP_Res_", resolutions, ".tif"))
  outputs <- c(p$input, p$provenance, p$coverage, plots)
  if ("--dry-run" %in% args) {
    ok <- compact_dry_run("Residual post-QC processing (paths only; no RDS load)", inputs, outputs,
                         checks = c(explicit_qc_choice = !no_removals || !file.exists(q$decision)))
    return(invisible(if (ok) 0L else 1L))
  }
  if (any(!file.exists(inputs))) stop("Missing post-QC inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  if (any(file.exists(outputs)) && !"--overwrite" %in% args) stop("Protected post-QC outputs exist; review before --overwrite.")
  state <- load_residual_qc_state(config)
  if (identical(state$active, no_removals)) stop("QC choice conflicts with the approved decision.")
  if (!identical(normalizePath(input), normalizePath(state$input))) stop("Unexpected QC source.")
  ids <- as.character(load_sample_manifest(config)$sample_id)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid sample manifest.")
  merge <- read.csv(merge_manifest, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "residual_cells") %in% names(merge)) ||
      nrow(merge) != length(ids) || anyDuplicated(merge$sample_id) || !setequal(merge$sample_id, ids)) stop("Invalid residual merge sample map.")
  merge <- merge[match(ids, merge$sample_id), , drop = FALSE]
  if (state$active && !identical(as.integer(merge$residual_cells), as.integer(state$manifest$input_cells))) stop("QC input counts differ from the merge manifest.")
  expected <- if (state$active) state$manifest$remaining_cells else merge$residual_cells
  if (any(!is.finite(expected)) || any(expected < 0 | expected != floor(expected))) stop("Invalid expected per-sample counts.")
  suppressPackageStartupMessages({library(Seurat); library(harmony); library(Cairo); library(ggplot2); library(future)})
  future::plan(future::sequential)
  source_info <- file.info(inputs)
  sample_md5 <- unname(tools::md5sum(resolve_config_path(config$manifests$samples, config)))
  object <- readRDS(input)
  validate_residual_qc_object(object, state)
  original_ids <- colnames(object)
  samples <- as.character(object$orig.ident)
  if (anyNA(samples) || any(!samples %in% ids) ||
      !identical(as.integer(table(factor(samples, levels = ids))), as.integer(expected))) stop("Retained cell counts differ from expected sample counts.")
  object <- residual_postqc_representations(object, config, no_removals)
  validate_residual_qc_object(object, state)
  object@misc$residual_post_qc <- list(reintegrated = !no_removals, input = normalizePath(input),
    decision_md5 = state$decision_md5, seed = config$runtime$random_seed)
  dir.create(dirname(p$provenance), recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(dirname(p$provenance), ".postqc_lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Residual post-QC output is locked; inspect running/failed jobs.")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (any(file.exists(outputs)) && !"--overwrite" %in% args) stop("Post-QC output appeared during computation; refusing replacement.")
  for (path in outputs) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(resolutions)) {
    plot <- DimPlot(object, reduction = "umap_harmony", group.by = paste0("residual_postQC_Xenium_snn_res.", resolutions[i]),
      label = TRUE, repel = TRUE, raster = TRUE) + ggtitle(paste("Residual post-QC resolution", resolutions[i]))
    Cairo::CairoTIFF(plots[i], width = 10, height = 8, units = "in", res = 600)
    print(plot); grDevices::dev.off()
  }
  after <- file.info(inputs)
  if (!identical(after$size, source_info$size) || !identical(after$mtime, source_info$mtime) ||
      !identical(load_residual_qc_state(config)$decision_md5, state$decision_md5) ||
      !identical(unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))), sample_md5)) stop("QC inputs or sample manifest changed during processing.")
  saveRDS(object, p$input, compress = FALSE)
  write.csv(data.frame(sample_id = ids, remaining_cells = as.integer(expected)), p$coverage, row.names = FALSE)
  info <- file.info(p$input)
  write.csv(data.frame(input_rds = normalizePath(input), input_size = source_info$size[1], input_mtime = as.numeric(source_info$mtime[1]),
    output_rds = normalizePath(p$input), output_size = info$size, output_mtime = as.numeric(info$mtime),
    cells_after_qc = length(original_ids), resolutions = paste(resolutions, collapse = "|"), random_seed = config$runtime$random_seed,
    decision_md5 = state$decision_md5, reintegrated = !no_removals, cluster_prefix = "residual_postQC_Xenium_snn_res.",
    dimensions = config$regional_subsets$integration_dimensions_default, neighbor_k = config$regional_subsets$neighbor_k_post_qc,
    sample_manifest_md5 = sample_md5, coverage_md5 = unname(tools::md5sum(p$coverage)), stringsAsFactors = FALSE), p$provenance, row.names = FALSE)
  load_residual_analysis_state(config)
  message("Residual post-QC processing complete. Regenerate marker evidence and approve fresh labels: ", p$input)
  invisible(0L)
}
if (sys.nframe() == 0L) quit(save = "no", status = main())
