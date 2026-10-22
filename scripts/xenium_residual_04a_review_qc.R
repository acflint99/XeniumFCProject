#!/usr/bin/env Rscript
# Independent descriptive QC report; never writes an RDS or applies a decision.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_04a_review_qc.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_residual_qc.R"))
source(here("scripts", "R", "xenium_replicate_marker_evidence.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
ids <- as.character(samples$sample_id)
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid sample manifest IDs.")
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite"))) {
  stop("Usage: Rscript scripts/xenium_residual_04a_review_qc.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
root <- here(config$project$outputs_dir, "xenium", "residual")
input <- file.path(root, "03_integrated", "rds", "Xenium_Residual_Integrated.rds")
merge_manifest <- file.path(root, "02_merged", "tables", "Xenium_Residual_Merge_Manifest.csv")
out <- file.path(root, "04_qc_review")
outputs <- c(
  cells = file.path(out, "tables", "Xenium_Residual_QC_PerCell.csv"),
  pooled = file.path(out, "tables", "Xenium_Residual_QC_ByCluster.csv"),
  paired = file.path(out, "tables", "Xenium_Residual_QC_BySampleCluster.csv"),
  balanced = file.path(out, "tables", "Xenium_Residual_QC_SampleBalanced.csv"),
  segmentation = file.path(out, "tables", "Xenium_Residual_QC_Segmentation.csv"),
  coverage = file.path(out, "tables", "Xenium_Residual_QC_SampleCoverage.csv"),
  provenance = file.path(out, "tables", "Xenium_Residual_QC_Provenance.csv"),
  plot = file.path(out, "plots", "Xenium_Residual_QC_AllClusters.pdf")
)
if (!file.exists(merge_manifest)) {
  if (dry_run) {
    compact_dry_run("Residual QC (raw paths require merge manifest)",
                    inputs = c(input, merge_manifest), outputs = outputs)
    quit(save = "no", status = 1L)
  }
  stop("Missing residual merge manifest: ", merge_manifest)
}
m <- read.csv(merge_manifest, stringsAsFactors = FALSE, na.strings = character())
if (!all(c("sample_id", "residual_cells", "has_residual_cells") %in% names(m)) ||
    nrow(m) != length(ids) || anyDuplicated(m$sample_id) || !setequal(m$sample_id, ids)) {
  stop("Merge manifest must map one-to-one to the current sample manifest.")
}
m <- m[match(ids, m$sample_id), , drop = FALSE]
expected <- residual_qc_numeric(m$residual_cells, "manifest residual_cells")
if (any(expected != round(expected)) || anyNA(as.logical(m$has_residual_cells)) ||
    any(as.logical(m$has_residual_cells) != (expected > 0)) || !any(expected > 0)) {
  stop("Invalid residual counts/statuses in merge manifest.")
}
raw_paths <- setNames(vapply(which(expected > 0), function(i) {
  directory <- resolve_sample_paths(samples[i, , drop = FALSE], config)$input_dir
  candidates <- file.path(directory, c("cells.csv.gz", "cells.csv"))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 1L) stop("Ambiguous raw cell tables for ", ids[[i]])
  if (length(existing)) existing[[1]] else candidates[[1]]
}, character(1)), ids[expected > 0])
inputs <- c(input, merge_manifest, raw_paths)
if (dry_run) {
  compact_dry_run("Residual cluster QC (paths/CSV manifest only; no RDS/raw reads)",
                  inputs = inputs, outputs = outputs)
  quit(save = "no", status = if (all(file.exists(inputs))) 0L else 1L)
}
if (any(!file.exists(inputs))) stop("Missing QC inputs:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
if (any(file.exists(outputs)) && !overwrite) {
  stop("Refusing to overwrite residual QC outputs. Review them before using --overwrite.")
}
suppressPackageStartupMessages({library(SeuratObject); library(ggplot2)})
input_info <- file.info(inputs)
object <- readRDS(input)
metadata <- object[[]]
cluster_column <- "Xenium_snn_res.0.8"
required <- c("orig.ident", cluster_column, "consensus_label", "residual_reference_route",
              "nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control")
if (!all(required %in% names(metadata))) {
  stop("Residual metadata lacks: ", paste(setdiff(required, names(metadata)), collapse = ", "))
}
if (!identical(colnames(object), rownames(metadata)) || anyDuplicated(rownames(metadata)) ||
    nrow(metadata) != sum(expected)) stop("Residual cell IDs/counts differ from expected metadata.")
cells <- data.frame(cell_id = rownames(metadata), sample_id = as.character(metadata$orig.ident),
                    cluster = as.character(metadata[[cluster_column]]), stringsAsFactors = FALSE)
if (anyNA(cells$sample_id) || any(!cells$sample_id %in% ids) ||
    anyNA(cells$cluster) || any(!nzchar(trimws(cells$cluster)))) stop("Invalid sample/cluster identities.")
observed <- as.integer(table(factor(cells$sample_id, levels = ids)))
if (!identical(observed, as.integer(expected))) stop("Per-sample cell counts differ from merge manifest.")
if (anyNA(metadata$consensus_label) || any(!is_unknown_consensus_label(metadata$consensus_label)) ||
    anyNA(metadata$residual_reference_route) ||
    any(!metadata$residual_reference_route %in% c("route_other", "route_unresolved"))) {
  stop("QC input contains a non-Unknown or non-residual cell.")
}
for (metric in c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control")) {
  cells[[metric]] <- residual_qc_numeric(metadata[[metric]], metric)
}
if (any(cells$cell_area <= 0) || any(cells$percent_control > 100)) stop("Invalid area/control percentage.")
for (metric in c("nCount_Xenium", "nFeature_Xenium")) {
  if (any(cells[[metric]] != round(cells[[metric]]))) stop("Noninteger ", metric, ".")
}
cells$nucleus_area <- cells$nucleus_count <- NA_real_
max_difference <- numeric(length(raw_paths))
for (i in seq_along(raw_paths)) {
  id <- names(raw_paths)[[i]]
  target <- which(cells$sample_id == id)
  prefix <- paste0(id, "_")
  if (any(!startsWith(cells$cell_id[target], prefix))) stop("Invalid sample prefix for ", id)
  raw_ids <- substring(cells$cell_id[target], nchar(prefix) + 1L)
  message("Joining raw morphology for ", id, " (", length(target), " cells)")
  morphology <- residual_qc_read_morphology(raw_paths[[i]], raw_ids, cells$cell_area[target])
  max_difference[[i]] <- attr(morphology, "max_cell_area_difference")
  cells$cell_area[target] <- morphology$cell_area
  cells$nucleus_area[target] <- morphology$nucleus_area
  cells$nucleus_count[target] <- morphology$nucleus_count
}
cells$cell_to_nucleus_area_ratio <- ifelse(
  is.finite(cells$nucleus_area) & cells$nucleus_area > 0,
  cells$cell_area / cells$nucleus_area, NA_real_)
cells$transcripts_per_cell_area <- cells$nCount_Xenium / cells$cell_area
metrics <- c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control",
             "nucleus_area", "nucleus_count", "cell_to_nucleus_area_ratio",
             "transcripts_per_cell_area")
# Reuse the established marker-review eligibility rule, not a filtering threshold.
min_cells <- eval(formals(build_replicate_marker_evidence)$min_cells_per_group)
summary <- residual_qc_summarize(cells, metrics, min_cells)
coverage <- data.frame(task_id = seq_along(ids), sample_id = ids,
                       expected_cells = expected, observed_cells = observed,
                       status = ifelse(expected > 0, "included", "no_residual_cells"))
provenance <- data.frame(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_path = normalizePath(inputs, mustWork = TRUE),
  input_size = input_info$size, input_mtime = as.numeric(input_info$mtime),
  cluster_column = cluster_column, cells = nrow(cells),
  clusters = length(unique(cells$cluster)), manifest_samples = length(ids),
  samples_with_cells = sum(expected > 0), min_measured_target_and_background = min_cells,
  background = "other residual clusters in the same biological sample",
  heatmap_reference = "all residual cells in the same biological sample (matches RL/VZ convention)",
  inference = "descriptive; samples equally weighted in sample-balanced summaries",
  cell_area_source = "raw Xenium cells table; retained area validated before replacement",
  cell_area_rounding_tolerance = 0.00501,
  maximum_cell_area_difference = max(max_difference),
  missing_values = "nucleus_area may be missing; ratio undefined for missing/zero nucleus area",
  decisions_applied = FALSE, stringsAsFactors = FALSE)
provenance$sample_id <- c("", "", names(raw_paths))
provenance$config_path <- normalizePath(attr(config, "config_path"), mustWork = TRUE)
provenance$sample_manifest_path <- resolve_config_path(config$manifests$samples, config)
after <- file.info(inputs)
if (!identical(as.numeric(input_info$size), as.numeric(after$size)) ||
    !identical(as.numeric(input_info$mtime), as.numeric(after$mtime))) stop("QC inputs changed during analysis.")

# Build complete outputs in a private directory; publish only after plot generation succeeds.
dir.create(out, recursive = TRUE, showWarnings = FALSE)
temporary <- tempfile(".qc_build_", tmpdir = out)
dir.create(temporary)
staged <- file.path(temporary, basename(outputs)); names(staged) <- names(outputs)
write.csv(cells, staged[["cells"]], row.names = FALSE, na = "")
for (name in names(summary)) write.csv(summary[[name]], staged[[name]], row.names = FALSE, na = "")
write.csv(coverage, staged[["coverage"]], row.names = FALSE)
write.csv(provenance, staged[["provenance"]], row.names = FALSE)
cluster_levels <- unique(summary$pooled$cluster)
cells$cluster <- factor(cells$cluster, levels = cluster_levels)
grDevices::pdf(staged[["plot"]], width = 12, height = 9)
tryCatch({
  for (metric in metrics) {
    plot_data <- cells[is.finite(cells[[metric]]), , drop = FALSE]
    if (!nrow(plot_data)) {
      plot <- ggplot() + annotate("text", x = 0, y = 0, label = paste(metric, "has no measured values")) + theme_void()
      print(plot); next
    }
    plot_data$qc_value <- plot_data[[metric]]
    plot <- ggplot(plot_data, aes(x = cluster, y = qc_value)) +
      geom_boxplot(outlier.shape = NA) + theme_bw() +
      labs(x = "Residual cluster", y = metric,
           subtitle = "Descriptive distributions; boxes omit outlier points, not cells")
    print(plot + labs(title = paste("Residual QC:", metric),
                     caption = "Pooled cells can reflect sample composition; use within-sample panels for review."))
    print(plot + facet_wrap(~sample_id, ncol = 3, scales = "free_y") +
            labs(title = paste("Within-sample residual QC:", metric)))
  }
  # Match RL/VZ review: sample-relative QC and sample-composition heatmaps.
  for (metric in c("nCount_Xenium", "nFeature_Xenium")) {
    heat <- summary$paired[summary$paired$metric == metric, , drop = FALSE]
    heat$cluster <- factor(heat$cluster, levels = cluster_levels)
    heat$sample_id <- factor(heat$sample_id, levels = rev(ids))
    heat$log2_ratio <- ifelse(is.finite(heat$median_ratio_to_sample) &
                               heat$median_ratio_to_sample > 0,
                             log2(heat$median_ratio_to_sample), NA_real_)
    print(ggplot(heat, aes(x = cluster, y = sample_id, fill = log2_ratio)) +
      geom_tile(color = "white", linewidth = .15) +
      scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, na.value = "grey90", name = "log2 ratio") +
      theme_bw(base_size = 9) +
      labs(title = paste("Residual QC: within-sample", metric, "ratio"),
           subtitle = "Cluster median / median of all residual cells in the same sample",
           caption = "Blank tiles: no cells. Grey tiles: undefined ratio; not a QC decision.",
           x = "Residual cluster", y = "Manifest sample"))
  }
  heat <- summary$paired[summary$paired$metric == "nCount_Xenium", , drop = FALSE]
  heat$cluster <- factor(heat$cluster, levels = cluster_levels)
  heat$sample_id <- factor(heat$sample_id, levels = rev(ids))
  print(ggplot(heat, aes(x = cluster, y = sample_id, fill = fraction_of_cluster)) +
    geom_tile(color = "white", linewidth = .15) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    scale_fill_gradient(low = "white", high = "#2166AC", name = "Fraction") +
    theme_bw(base_size = 9) +
    labs(title = "Residual subclusters: sample composition",
         subtitle = "Fraction of each cluster's cells contributed by each sample",
         caption = "Blank tiles: no cells; all manifest samples are shown.",
         x = "Residual cluster", y = "Manifest sample"))
}, finally = grDevices::dev.off())
for (name in c(setdiff(names(outputs), "provenance"), "provenance")) {
  dir.create(dirname(outputs[[name]]), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(outputs[[name]]) && !overwrite) stop("QC output appeared during analysis: ", outputs[[name]])
  if (!file.rename(staged[[name]], outputs[[name]])) stop("Could not publish ", outputs[[name]])
}
unlink(temporary, recursive = TRUE)
message("Residual QC report complete: ", out, ". No cells removed, labels applied, or RDS files written.")
