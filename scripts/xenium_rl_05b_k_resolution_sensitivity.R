#!/usr/bin/env Rscript

# Non-destructive RL graph sensitivity analysis. The existing post-QC
# expression reduction and UMAP are held fixed while k.param and clustering
# resolution vary.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "subcluster_review_markers.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_rl_05b_k_resolution_sensitivity.R ",
    "[--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

k_values <- c(10L, 20L)
resolutions <- c(0.8, 1.0)
dimensions <- 30L
prune_snn <- 1 / 15
cluster_algorithm <- 1L
cluster_starts <- 10L
cluster_iterations <- 10L
focus_production_cluster <- "6"
production_cluster_column <- "Xenium_snn_res.0.8"

resolution_label <- function(value) {
  if (isTRUE(all.equal(value, 1))) "1.0" else format(value, trim = TRUE)
}
combo_tag <- function(k_value, resolution) {
  paste0("k", k_value, "_Res", resolution_label(resolution))
}
diagnostic_column <- function(k_value, resolution) {
  paste0("RL_diag_k", k_value, "_res", resolution_label(resolution))
}

output_root <- here(config$project$outputs_dir)
input_path <- file.path(
  output_root, "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC.rds"
)
input_manifest_path <- file.path(
  output_root, "xenium", "rl", "05_post_qc", "tables", "Xenium_RL_postQC_manifest.csv"
)
diagnostic_root <- file.path(
  output_root, "xenium", "rl", "05_post_qc", "diagnostics", "k_resolution_sensitivity"
)
plot_dir <- file.path(diagnostic_root, "plots")
table_dir <- file.path(diagnostic_root, "tables")

common_table_names <- c(
  "RL_k10_k20_Res0.8_1.0_Assignments.csv.gz",
  "RL_k10_k20_Res0.8_1.0_Cluster_Summary.csv",
  "RL_k10_k20_Res0.8_1.0_Sample_Composition.csv",
  "RL_k10_k20_Res0.8_1.0_InputLabel_Composition.csv",
  "RL_k10_k20_Res0.8_1.0_Contingencies.csv",
  "RL_ProductionCluster6_Splits.csv",
  "RL_k10_k20_Res0.8_1.0_Manifest.csv"
)
combo_grid <- expand.grid(
  k = k_values,
  resolution = resolutions,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
combo_plot_names <- unlist(lapply(seq_len(nrow(combo_grid)), function(i) {
  tag <- combo_tag(combo_grid$k[[i]], combo_grid$resolution[[i]])
  c(
    paste0("RL_", tag, "_UMAP.tif"),
    paste0("RL_ProductionCluster6_", tag, "_UMAP.tif"),
    paste0("RL_", tag, "_FocusedMarker_DotPlot.tif"),
    paste0("RL_", tag, "_FocusedMarker_DotPlot.pdf")
  )
}), use.names = FALSE)
combo_table_names <- unlist(lapply(seq_len(nrow(combo_grid)), function(i) {
  tag <- combo_tag(combo_grid$k[[i]], combo_grid$resolution[[i]])
  c(
    paste0("RL_", tag, "_top50_Markers.csv"),
    paste0("RL_", tag, "_FocusedMarker_Availability.csv")
  )
}), use.names = FALSE)
expected_outputs <- c(
  file.path(plot_dir, combo_plot_names),
  file.path(table_dir, c(common_table_names, combo_table_names))
)

if (dry_run) {
  compact_dry_run(
    "RL k.param/resolution sensitivity",
    inputs = c(input_path, input_manifest_path),
    outputs = expected_outputs,
    checks = c(
      sample_manifest_has_34_samples = length(sample_ids) == 34L,
      k_values_are_10_and_20 = identical(k_values, c(10L, 20L)),
      resolutions_are_0.8_and_1.0 = identical(resolutions, c(0.8, 1.0))
    )
  )
  quit(save = "no", status = 0L)
}

missing_inputs <- c(input_path, input_manifest_path)[
  !file.exists(c(input_path, input_manifest_path))
]
if (length(missing_inputs)) {
  stop("Missing sensitivity input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
if (length(sample_ids) != 34L || anyDuplicated(sample_ids) || any(!nzchar(sample_ids))) {
  stop("config/samples.csv must contain exactly 34 unique, nonblank sample IDs.")
}
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL sensitivity outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(ggplot2)
})
source(here("scripts", "color_palette.R"))

review_marker_sets <- add_palette_markers_to_review(
  focused_marker_sets = rl_focused_review_markers,
  palette_marker_sets = rl_markers,
  palette_object_name = "rl_markers"
)
set.seed(config$runtime$random_seed)
plan(sequential)
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)

check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message("\n[", Sys.time(), "] --- ", step_label, " ---")
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

save_tiff <- function(plot, path, width = 12, height = 9) {
  Cairo::CairoTIFF(
    filename = path, width = width, height = height,
    units = "in", res = 600
  )
  print(plot)
  grDevices::dev.off()
}

staging_root <- file.path(
  diagnostic_root,
  paste0(".staging_", Sys.getenv("SLURM_JOB_ID", unset = "local"), "_", Sys.getpid())
)
staging_plot_dir <- file.path(staging_root, "plots")
staging_table_dir <- file.path(staging_root, "tables")
if (dir.exists(staging_root)) stop("Sensitivity staging directory already exists: ", staging_root)
dir.create(staging_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_table_dir, recursive = TRUE, showWarnings = FALSE)
completed <- FALSE
on.exit({
  if (!completed && dir.exists(staging_root)) {
    message("Sensitivity run failed; partial files remain only in staging: ", staging_root)
  }
}, add = TRUE)

input_manifest <- read.csv(input_manifest_path, stringsAsFactors = FALSE)
required_manifest_columns <- c(
  "output_rds", "output_size", "output_mtime", "random_seed", "processing_mode",
  "pca_dimensions", "neighbor_k"
)
if (nrow(input_manifest) != 1L ||
    !all(required_manifest_columns %in% names(input_manifest))) {
  stop("RL Stage-5 manifest is invalid: ", input_manifest_path)
}
input_info <- file.info(input_path)
manifest_matches <-
  identical(
    normalizePath(as.character(input_manifest$output_rds), mustWork = TRUE),
    normalizePath(input_path, mustWork = TRUE)
  ) &&
  isTRUE(all.equal(as.numeric(input_manifest$output_size), as.numeric(input_info$size))) &&
  isTRUE(all.equal(as.numeric(input_manifest$output_mtime), as.numeric(input_info$mtime))) &&
  identical(as.integer(input_manifest$random_seed), as.integer(config$runtime$random_seed))
if (!manifest_matches) {
  stop("RL Stage-5 RDS does not match its provenance manifest.")
}
if (!identical(as.integer(input_manifest$pca_dimensions), dimensions) ||
    !identical(as.integer(input_manifest$neighbor_k), 30L)) {
  stop(
    "This controlled sensitivity analysis requires the production RL Stage-5 ",
    "object built with 30 dimensions and k.param=30."
  )
}

message(Sys.time(), ": loading RL Stage-5 object...")
obj <- readRDS(input_path)
check_mem("RL STAGE-5 OBJECT LOADED")
if (anyDuplicated(colnames(obj))) stop("RL Stage-5 object contains duplicate cell IDs.")
if (!"Xenium" %in% Assays(obj)) stop("RL Stage-5 object lacks the Xenium assay.")
required_metadata <- c(
  "orig.ident", "regional_consensus_label_pooled", production_cluster_column
)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) {
  stop("RL Stage-5 object lacks metadata: ", paste(missing_metadata, collapse = ", "))
}
unexpected_samples <- setdiff(unique(as.character(obj$orig.ident)), sample_ids)
if (length(unexpected_samples)) {
  stop("RL Stage-5 object contains samples absent from config/samples.csv: ",
       paste(unexpected_samples, collapse = ", "))
}
if (!"umap_clean" %in% Reductions(obj)) stop("RL Stage-5 object lacks umap_clean.")
clustering_reduction <- if ("harmony_clean" %in% Reductions(obj)) {
  "harmony_clean"
} else if ("harmony" %in% Reductions(obj)) {
  "harmony"
} else {
  stop("RL Stage-5 object lacks harmony_clean and harmony reductions.")
}
if (ncol(Embeddings(obj, clustering_reduction)) < dimensions) {
  stop(clustering_reduction, " has fewer than ", dimensions, " dimensions.")
}
if (!focus_production_cluster %in%
    unique(as.character(obj[[production_cluster_column, drop = TRUE]]))) {
  stop("Production focus cluster is absent: ", focus_production_cluster)
}

assignments <- data.frame(
  cell_id = colnames(obj),
  production_res0.8 = as.character(obj[[production_cluster_column, drop = TRUE]]),
  stringsAsFactors = FALSE
)
cluster_summaries <- list()
sample_compositions <- list()
label_compositions <- list()
focus_splits <- list()

composition_table <- function(cluster_values, category_values, k_value, resolution, category_name) {
  counts <- as.data.frame(
    table(cluster = cluster_values, category = category_values),
    stringsAsFactors = FALSE
  )
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  names(counts)[names(counts) == "Freq"] <- "cell_count"
  counts$fraction_of_cluster <- counts$cell_count /
    ave(counts$cell_count, counts$cluster, FUN = sum)
  counts$k_param <- k_value
  counts$resolution <- resolution
  counts$category_type <- category_name
  counts
}

for (k_value in k_values) {
  nn_name <- paste0("RL_sensitivity_nn_k", k_value)
  snn_name <- paste0("RL_sensitivity_snn_k", k_value)
  message(Sys.time(), ": building dimensions 1:", dimensions,
          " graph with k.param=", k_value, "...")
  obj <- FindNeighbors(
    obj,
    reduction = clustering_reduction,
    dims = seq_len(dimensions),
    k.param = k_value,
    prune.SNN = prune_snn,
    graph.name = c(nn_name, snn_name),
    verbose = TRUE
  )
  obj <- FindClusters(
    obj,
    graph.name = snn_name,
    resolution = resolutions,
    algorithm = cluster_algorithm,
    n.start = cluster_starts,
    n.iter = cluster_iterations,
    random.seed = config$runtime$random_seed,
    verbose = TRUE
  )

  for (resolution in resolutions) {
    source_resolution_label <- format(resolution, scientific = FALSE, trim = TRUE)
    source_column <- paste0(snn_name, "_res.", source_resolution_label)
    if (!source_column %in% colnames(obj[[]])) {
      stop("FindClusters did not create expected metadata: ", source_column)
    }
    output_column <- diagnostic_column(k_value, resolution)
    cluster_values <- as.character(obj[[source_column, drop = TRUE]])
    numeric_levels <- as.character(sort(unique(as.numeric(cluster_values))))
    obj[[output_column]] <- factor(cluster_values, levels = numeric_levels)
    assignments[[output_column]] <- cluster_values
    tag <- combo_tag(k_value, resolution)

    cluster_summary <- as.data.frame(table(cluster = cluster_values), stringsAsFactors = FALSE)
    names(cluster_summary)[names(cluster_summary) == "Freq"] <- "cell_count"
    cluster_summary$sample_count <- vapply(
      cluster_summary$cluster,
      function(cluster) length(unique(as.character(obj$orig.ident)[cluster_values == cluster])),
      integer(1)
    )
    cluster_summary$k_param <- k_value
    cluster_summary$resolution <- resolution
    cluster_summaries[[tag]] <- cluster_summary
    sample_compositions[[tag]] <- composition_table(
      cluster_values, as.character(obj$orig.ident),
      k_value, resolution, "orig.ident"
    )
    label_compositions[[tag]] <- composition_table(
      cluster_values, as.character(obj$regional_consensus_label_pooled),
      k_value, resolution, "regional_consensus_label_pooled"
    )

    focus_cells <- assignments$production_res0.8 == focus_production_cluster
    focus_split <- as.data.frame(
      table(
        diagnostic_cluster = cluster_values[focus_cells],
        sample_id = as.character(obj$orig.ident)[focus_cells]
      ),
      stringsAsFactors = FALSE
    )
    focus_split <- focus_split[focus_split$Freq > 0L, , drop = FALSE]
    names(focus_split)[names(focus_split) == "Freq"] <- "cell_count"
    focus_split$fraction_of_production_cluster <- focus_split$cell_count /
      sum(focus_split$cell_count)
    focus_split$k_param <- k_value
    focus_split$resolution <- resolution
    focus_splits[[tag]] <- focus_split

    umap_plot <- DimPlot(
      obj, reduction = "umap_clean", group.by = output_column,
      label = TRUE, label.size = 4, repel = TRUE, raster = TRUE
    ) +
      ggtitle(paste0("RL sensitivity: k=", k_value, ", resolution=", resolution_label(resolution))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_tiff(umap_plot, file.path(staging_plot_dir, paste0("RL_", tag, "_UMAP.tif")))

    focus_column <- paste0(output_column, "_production_cluster6")
    obj[[focus_column]] <- ifelse(
      focus_cells, paste0("Cluster ", cluster_values), "Other RL cells"
    )
    focus_levels <- c(
      "Other RL cells",
      paste0("Cluster ", numeric_levels[numeric_levels %in% unique(cluster_values[focus_cells])])
    )
    obj[[focus_column]] <- factor(obj[[focus_column, drop = TRUE]], levels = focus_levels)
    focus_colors <- c(
      "Other RL cells" = "grey90",
      stats::setNames(
        grDevices::hcl.colors(length(focus_levels) - 1L, "Dynamic"),
        focus_levels[-1L]
      )
    )
    focus_plot <- DimPlot(
      obj, reduction = "umap_clean", group.by = focus_column,
      cols = focus_colors, raster = TRUE, pt.size = 0.5
    ) +
      ggtitle(paste0(
        "Production cluster 6 split: k=", k_value,
        ", resolution=", resolution_label(resolution)
      )) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_tiff(
      focus_plot,
      file.path(staging_plot_dir, paste0("RL_ProductionCluster6_", tag, "_UMAP.tif"))
    )

    Idents(obj) <- output_column
    message(Sys.time(), ": identifying markers for ", tag, "...")
    all_markers <- FindAllMarkers(
      obj,
      assay = "Xenium",
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident = 1000,
      random.seed = config$runtime$random_seed
    )
    if (!nrow(all_markers)) stop("No markers returned for ", tag, ".")
    top50_markers <- all_markers %>%
      group_by(.data$cluster) %>%
      slice_max(n = 50, order_by = .data$avg_log2FC, with_ties = FALSE) %>%
      ungroup()
    write.csv(
      top50_markers,
      file.path(staging_table_dir, paste0("RL_", tag, "_top50_Markers.csv")),
      row.names = FALSE
    )

    write_focused_marker_dotplot(
      object = obj,
      marker_sets = review_marker_sets,
      assay = "Xenium",
      group_by = output_column,
      branch = paste0("RL k=", k_value),
      resolution_label = resolution_label(resolution),
      availability_path = file.path(
        staging_table_dir, paste0("RL_", tag, "_FocusedMarker_Availability.csv")
      ),
      tiff_path = file.path(
        staging_plot_dir, paste0("RL_", tag, "_FocusedMarker_DotPlot.tif")
      ),
      pdf_path = file.path(
        staging_plot_dir, paste0("RL_", tag, "_FocusedMarker_DotPlot.pdf")
      ),
      width = 21,
      height = 10
    )
    obj[[focus_column]] <- NULL
    rm(umap_plot, focus_plot, all_markers, top50_markers)
    gc()
  }

  obj[[nn_name]] <- NULL
  obj[[snn_name]] <- NULL
  check_mem(paste0("COMPLETED k.param=", k_value))
}

assignment_connection <- gzfile(
  file.path(staging_table_dir, common_table_names[[1]]), open = "wt"
)
write.csv(assignments, assignment_connection, row.names = FALSE)
close(assignment_connection)
write.csv(bind_rows(cluster_summaries), file.path(staging_table_dir, common_table_names[[2]]), row.names = FALSE)
write.csv(bind_rows(sample_compositions), file.path(staging_table_dir, common_table_names[[3]]), row.names = FALSE)
write.csv(bind_rows(label_compositions), file.path(staging_table_dir, common_table_names[[4]]), row.names = FALSE)

make_contingency <- function(from_column, to_column, comparison) {
  result <- as.data.frame(
    table(from_cluster = assignments[[from_column]], to_cluster = assignments[[to_column]]),
    stringsAsFactors = FALSE
  )
  result <- result[result$Freq > 0L, , drop = FALSE]
  names(result)[names(result) == "Freq"] <- "cell_count"
  result$fraction_of_from_cluster <- result$cell_count /
    ave(result$cell_count, result$from_cluster, FUN = sum)
  result$comparison <- comparison
  result$from_column <- from_column
  result$to_column <- to_column
  result
}

contingencies <- list()
for (resolution in resolutions) {
  contingencies[[paste0("k_at_res_", resolution_label(resolution))]] <- make_contingency(
    diagnostic_column(10L, resolution),
    diagnostic_column(20L, resolution),
    paste0("k10_vs_k20_at_resolution_", resolution_label(resolution))
  )
}
for (k_value in k_values) {
  contingencies[[paste0("resolution_at_k_", k_value)]] <- make_contingency(
    diagnostic_column(k_value, 0.8),
    diagnostic_column(k_value, 1.0),
    paste0("resolution_0.8_vs_1.0_at_k", k_value)
  )
}
for (i in seq_len(nrow(combo_grid))) {
  k_value <- combo_grid$k[[i]]
  resolution <- combo_grid$resolution[[i]]
  contingencies[[paste0("production_vs_", combo_tag(k_value, resolution))]] <- make_contingency(
    "production_res0.8",
    diagnostic_column(k_value, resolution),
    paste0("production_k30_res0.8_vs_", combo_tag(k_value, resolution))
  )
}
write.csv(
  bind_rows(contingencies),
  file.path(staging_table_dir, common_table_names[[5]]),
  row.names = FALSE
)
write.csv(
  bind_rows(focus_splits),
  file.path(staging_table_dir, common_table_names[[6]]),
  row.names = FALSE
)

manifest <- data.frame(
  input_rds = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime),
  input_stage5_manifest = normalizePath(input_manifest_path, mustWork = TRUE),
  input_processing_mode = as.character(input_manifest$processing_mode),
  cells = ncol(obj),
  samples = length(unique(as.character(obj$orig.ident))),
  clustering_reduction = clustering_reduction,
  dimensions = dimensions,
  k_values = paste(k_values, collapse = "|"),
  resolutions = paste(vapply(resolutions, resolution_label, character(1)), collapse = "|"),
  prune_snn = prune_snn,
  cluster_algorithm = cluster_algorithm,
  cluster_starts = cluster_starts,
  cluster_iterations = cluster_iterations,
  random_seed = config$runtime$random_seed,
  fixed_umap = "umap_clean",
  production_cluster_column = production_cluster_column,
  focus_production_cluster = focus_production_cluster,
  marker_min_pct = 0.25,
  marker_logfc_threshold = 0.25,
  marker_max_cells_per_ident = 1000,
  focused_marker_panel = "rl_focused_review_markers + color_palette.R::rl_markers",
  focused_marker_signature = focused_marker_set_signature(review_marker_sets),
  writes_production_rds = FALSE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(staging_table_dir, common_table_names[[7]]),
  row.names = FALSE
)

staged_outputs <- c(
  file.path(staging_plot_dir, combo_plot_names),
  file.path(staging_table_dir, c(common_table_names, combo_table_names))
)
missing_staged_outputs <- staged_outputs[!file.exists(staged_outputs)]
if (length(missing_staged_outputs)) {
  stop("Sensitivity run did not create expected output(s):\n- ",
       paste(missing_staged_outputs, collapse = "\n- "))
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
copy_ok <- file.copy(staged_outputs, expected_outputs, overwrite = overwrite)
if (!all(copy_ok)) {
  stop("Failed to publish sensitivity output(s):\n- ",
       paste(expected_outputs[!copy_ok], collapse = "\n- "))
}
unlink(staging_root, recursive = TRUE)
completed <- TRUE
plan(sequential)
check_mem("RL k.param/RESOLUTION SENSITIVITY COMPLETE")
message("RL k.param/resolution sensitivity outputs written under: ", diagnostic_root)
