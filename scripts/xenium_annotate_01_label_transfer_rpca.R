# Manifest-driven annotation for every biological sample and reference.
# Use --list or --dry-run before launching heavy Seurat work.

rm(list = ls())

library(here)

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
sample_list <- sample_manifest$sample_id

pilot_manifest <- load_resolution2_pilot_manifest(config)
required_pilot_columns <- c("task_id", "sample_id", "PCW", "age_group")
if (!all(required_pilot_columns %in% names(pilot_manifest))) {
  stop(
    "Resolution-2.0 pilot manifest must contain: ",
    paste(required_pilot_columns, collapse = ", ")
  )
}
pilot_manifest_counts <- vapply(
  pilot_manifest$sample_id,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(pilot_manifest_counts != 1L)) {
  stop("Every clustering-pilot sample must map to exactly one config/samples.csv row.")
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

args <- commandArgs(trailingOnly = TRUE)

valid_options <- c(
  "--dry-run", "--overwrite", "--pilot-res2", "--pilot-res3",
  "--pilot-res4", "--pilot-res5", "--all-samples-res4", "--list"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

pilot_res2 <- "--pilot-res2" %in% args
pilot_res3 <- "--pilot-res3" %in% args
pilot_res4 <- "--pilot-res4" %in% args
pilot_res5 <- "--pilot-res5" %in% args
all_samples_res4 <- "--all-samples-res4" %in% args
pilot_flags <- c(pilot_res2, pilot_res3, pilot_res4, pilot_res5)
if (sum(pilot_flags) > 1L) {
  stop(
    "Choose only one of --pilot-res2, --pilot-res3, --pilot-res4, or --pilot-res5."
  )
}
pilot_mode <- any(pilot_flags)
if (all_samples_res4 && pilot_mode) {
  stop("--all-samples-res4 cannot be combined with a pilot-resolution option.")
}
resolution_mode <- pilot_mode || all_samples_res4
pilot_resolution <- if (all_samples_res4) 4.0 else if (pilot_res5) 5.0 else if (pilot_res4) 4.0 else if (pilot_res3) 3.0 else 2.0
pilot_resolution_tag <- sprintf("%.1f", pilot_resolution)
pilot_stage <- if (all_samples_res4) {
  "03g_resolution4_all_samples"
} else if (pilot_res5) {
  "03e_resolution5_pilot"
} else if (pilot_res4) {
  "03d_resolution4_pilot"
} else if (pilot_res3) {
  "03c_resolution3_pilot"
} else {
  "03b_resolution2_pilot"
}
pilot_cluster_column <- paste0("whole_tissue_cluster_res", pilot_resolution_tag)
pilot_graph_column <- paste0("Xenium_snn_res.", format(pilot_resolution, trim = TRUE))
facet_point_size <- if (all_samples_res4) {
  0.02
} else if (pilot_res3 || pilot_res4 || pilot_res5) {
  0.03
} else {
  0.1
}
list_requested <- "--list" %in% args
sample_list <- if (pilot_mode) pilot_manifest$sample_id else sample_manifest$sample_id
mode_description <- if (all_samples_res4) {
  "resolution-4.0 all-sample analysis"
} else if (pilot_mode) {
  paste0("resolution-", pilot_resolution_tag, " pilot")
} else {
  "production"
}

if (list_requested) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept REFERENCE or TASK_ID arguments.")
  cat(
    "Mode:", mode_description, "\n"
  )
  cat("References:", paste(names(reference_paths), collapse = ", "), "\n")
  write.table(
    data.frame(task_id = seq_along(sample_list), sample_id = sample_list),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) {
  stop("--dry-run and --overwrite cannot be used together.")
}

job_args <- args[!args %in% valid_options]
if (length(job_args) != 2L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01_label_transfer_rpca.R ",
    "[--pilot-res2|--pilot-res3|--pilot-res4|--pilot-res5|--all-samples-res4] ",
    "[--dry-run|--overwrite] REFERENCE TASK_ID"
  )
}

reference_match <- match(tolower(job_args[[1]]), tolower(names(reference_paths)))
if (is.na(reference_match)) {
  stop("REFERENCE must be Aldinger, Sepp, or Science.")
}
reference_name <- names(reference_paths)[reference_match]
reference_path <- unname(reference_paths[[reference_name]])
reference_key <- tolower(reference_name)

task_id <- suppressWarnings(as.integer(job_args[[2]]))
if (is.na(task_id) || task_id < 1L || task_id > length(sample_list)) {
  stop("TASK_ID must be between 1 and ", length(sample_list), ".")
}

current_sample <- sample_list[[task_id]]
if (resolution_mode) {
  pilot_root <- here(
    "outputs", "xenium", "preprocess", pilot_stage
  )
  input_file <- file.path(
    pilot_root, "rds",
    paste0(
      current_sample, "_whole_tissue_Res", pilot_resolution_tag,
      if (all_samples_res4) "" else "_pilot", ".rds"
    )
  )
  annotation_root <- if (all_samples_res4) {
    here("outputs", "xenium", "annotation", "resolution4_all_samples")
  } else {
    file.path(pilot_root, "annotation")
  }
} else {
  input_file <- here(
    "outputs", "xenium", "preprocess", "03_clustered", "rds",
    paste0(current_sample, "_CB_QC_cluster.rds")
  )
  annotation_root <- here("outputs", "xenium", "annotation")
}
annotation_dir <- file.path(annotation_root, "01_label_transfer", reference_key)
plots_dir <- file.path(annotation_dir, "plots")
tables_dir <- file.path(annotation_dir, "tables")
rds_dir <- file.path(annotation_dir, "rds")

expected_outputs <- c(
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_majority_vs_weighted.csv")),
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_prediction_cellcounts.csv")),
  file.path(plots_dir, paste0(current_sample, "_", reference_name, "_Broad_ClusterWeighted_UMAP.tif")),
  file.path(plots_dir, paste0(current_sample, "_", reference_name, "_Broad_ClusterMajority_UMAP.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_PredictionScores_Hist.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_PredictionScores_Hist.pdf")),
  file.path(plots_dir, paste0(current_sample, "_Broad_GlobalSpatial_ClusterWeighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_GlobalSpatial_ClusterMajority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_FacetSpatial_ClusterWeighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_FacetSpatial_ClusterMajority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Weighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Weighted.pdf")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Majority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Majority.pdf")),
  file.path(rds_dir, paste0(current_sample, "_", reference_name, "_annotated.rds"))
)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]

if (dry_run) {
  cat(
    "Mode:", mode_description, "\n"
  )
  cat("Reference:", reference_name, "\n")
  cat("Task:", task_id, "of", length(sample_list), "\n")
  cat("Sample:", current_sample, "\n")
  cat("Input:", input_file, "\n")
  cat("Input exists:", file.exists(input_file), "\n")
  cat("Reference RDS:", reference_path, "\n")
  cat("Reference exists:", file.exists(reference_path), "\n")
  cat("Dry-run scope: path inspection only; no RDS is loaded.\n")
  write.table(
    data.frame(output = expected_outputs, exists = file.exists(expected_outputs)),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  inputs_ready <- file.exists(input_file) && file.exists(reference_path)
  quit(save = "no", status = if (inputs_ready) 0L else 1L)
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (!file.exists(reference_path)) stop("Reference file not found: ", reference_path)

if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Refusing to overwrite existing ", reference_name,
    " annotation outputs for ", current_sample, ":\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}
if (length(existing_outputs) > 0L) {
  warning(
    "Overwriting ", length(existing_outputs), " existing ", reference_name,
    " annotation outputs for ", current_sample, "."
  )
}

library(future)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Cairo)

options(future.globals.maxSize = 400 * 1024^3)

message(
  "### Processing ", reference_name, " annotation for sample [",
  task_id, "/", length(sample_list), "]: ", current_sample, " ###"
)

# =========================================================
# Annotation Function
# =========================================================
annotate_xenium_from_ref <- function(xenium_obj,
                                     sample_name,
                                     reference_name,
                                     reference_path,
                                     annotation_dir) {
  
  ## ----------------------------
  ## 0. Setup & Paths
  ## ----------------------------
  set.seed(42)
  plan("sequential") 
  
  # Load custom palette and ordering from your specific script
  source(here("scripts", "color_palette.R")) 
  
  pred_score_thresh <- 0.4
  plots_dir <- file.path(annotation_dir, "plots")
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- file.path(annotation_dir, "tables")
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  ## ----------------------------
  ## 1. Load Reference
  ## ----------------------------
  ref_path <- reference_path
  if (!file.exists(ref_path)) stop("Reference file not found: ", ref_path)
  reference <- readRDS(ref_path)
  
  ## ----------------------------
  ## 2. Preparation & Anchor Finding
  ## ----------------------------
  # Set reference assay safely
  if ("RNA" %in% Assays(reference)) {
    DefaultAssay(reference) <- "RNA"
  } else {
    # Fallback to whatever assay exists (e.g., "originalexp")
    DefaultAssay(reference) <- Assays(reference)[1] 
  }
  
  # Set xenium assay safely
  if ("Xenium" %in% Assays(xenium_obj)) {
    DefaultAssay(xenium_obj) <- "Xenium"
  } else {
    DefaultAssay(xenium_obj) <- Assays(xenium_obj)[1]
  }
  
  reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  shared_genes <- intersect(rownames(reference), rownames(xenium_obj))
  
  # Downsample reference for speed and balance
  # Set the active identity to your cell types BEFORE downsampling
  Idents(reference) <- "clusters_refined"
  reference_balanced <- subset(reference, downsample = 1000)
  transfer_features  <- intersect(VariableFeatures(reference_balanced), shared_genes)
  
  reference_balanced <- ScaleData(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium_obj         <- ScaleData(xenium_obj, features = transfer_features, verbose = FALSE)
  
  reference_balanced <- RunPCA(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium_obj         <- RunPCA(xenium_obj, features = transfer_features, verbose = FALSE)
  
  cat("Shared features for transfer:", length(transfer_features), "\n")
  
  anchors <- FindTransferAnchors(
    reference = reference_balanced,
    query = xenium_obj,
    normalization.method = "LogNormalize", 
    reduction = "rpca", 
    features = transfer_features,
    dims = 1:30,
    k.anchor = 20,
    k.score = 30,
    approx.pca = TRUE
  )
  
  ## ----------------------------
  ## 3. Transfer & Voting Logic
  ## ----------------------------
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_balanced$clusters_refined,
    dims = 1:30,
    store.weights = TRUE
  )
  
  xenium_obj <- AddMetaData(xenium_obj, predictions)
  xenium_obj$high_conf <- xenium_obj$prediction.score.max > pred_score_thresh
  
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:30, reduction = "pca")
  
  # -------------------
  # Thresholded Majority Voting
  # -------------------
  majority_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1],
      mean_max_score = mean(prediction.score.max), # Calculate cluster average confidence
      .groups = "drop"
    ) %>%
    mutate(
      cluster_majority = ifelse(mean_max_score >= pred_score_thresh, cluster_majority, "Unknown")
    )
  
  xenium_obj$cluster_majority <- majority_labels$cluster_majority[match(xenium_obj$seurat_clusters, majority_labels$seurat_clusters)]
  
  # -------------------
  # Thresholded Rigorous Weighted Voting 
  # -------------------
  score_cols <- grep("^prediction\\.score\\.", colnames(xenium_obj@meta.data), value = TRUE)
  score_cols <- score_cols[!score_cols %in% c("prediction.score.max", "prediction.score.id")]
  
  # Get cluster sizes to normalize the total probability sums
  cluster_counts <- xenium_obj@meta.data %>% count(seurat_clusters, name = "n_cells")
  
  weighted_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), sum), .groups = "drop") %>%
    pivot_longer(cols = -seurat_clusters, names_to = "cluster_weighted", values_to = "total_score") %>%
    group_by(seurat_clusters) %>%
    slice_max(total_score, n = 1, with_ties = FALSE) %>%
    left_join(cluster_counts, by = "seurat_clusters") %>%
    mutate(
      cluster_weighted = sub("^prediction\\.score\\.", "", cluster_weighted),
      avg_prob = total_score / n_cells, # Calculate the mean probability per cell for the winning class
      cluster_weighted = ifelse(avg_prob >= pred_score_thresh, cluster_weighted, "Unknown")
    )
  
  xenium_obj$cluster_weighted <- weighted_labels$cluster_weighted[match(xenium_obj$seurat_clusters, weighted_labels$seurat_clusters)]
  
  ## ----------------------------
  ## 4. Export Comparison Table
  ## ----------------------------
  comparison <- xenium_obj@meta.data %>%
    select(seurat_clusters, cluster_majority, cluster_weighted) %>%
    distinct() %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(comparison, file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted.csv")), row.names = FALSE)
  
  ## ----------------------------
  ## 4c. Export Cluster-by-CellType Contingency Table
  ## ----------------------------
  # Note: Retains `predicted.id` so you can see the raw background noise/variance that caused a cluster to become "Unknown"
  count_matrix <- xenium_obj@meta.data %>%
    count(seurat_clusters, predicted.id) %>%
    pivot_wider(names_from = predicted.id, values_from = n, values_fill = list(n = 0)) %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(count_matrix, 
            file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")), 
            row.names = FALSE)
  
  message("Contingency matrix exported for: ", sample_name)
  
  ## ----------------------------
  ## 5. Visualizations
  ## ----------------------------
  
  # Safely add "Unknown" to your palette if it isn't already in color_palette.R
  if (!"Unknown" %in% celltype_order) {
    celltype_order <- c(celltype_order, "Unknown")
    if (!"Unknown" %in% names(cluster_colors)) {
      cluster_colors["Unknown"] <- "grey80" 
    }
  }
  
  # Validate palette (exclude "Unknown" from the strict check so it doesn't crash)
  clusters_to_validate <- setdiff(unique(c(xenium_obj$cluster_weighted, xenium_obj$cluster_majority)), "Unknown")
  validate_palette(clusters_to_validate)
  
  xenium_obj$cluster_weighted <- factor(xenium_obj$cluster_weighted, levels = celltype_order)
  xenium_obj$cluster_majority <- factor(xenium_obj$cluster_majority, levels = celltype_order)
  
  # --- UMAPs ---
  p_umap_wei <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterWeighted_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_wei)
  dev.off()
  
  p_umap_maj <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_majority", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterMajority_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_maj)
  dev.off()
  
  # --- Histogram of scores ---
  p_hist <- ggplot(xenium_obj@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
    theme_minimal() + labs(title = paste(sample_name, "Scores"), x = "Max Prediction Score")
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.tif")), width = 6, height = 4, units = "in", res = 600)
  print(p_hist)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.pdf")), p_hist, device = grDevices::cairo_pdf, width = 6, height = 4)
  
  # --- Spatial Plots ---
  p_wei <- ImageDimPlot(xenium_obj, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Weighted"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterWeighted.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_wei)
  dev.off() # Fixed missing parentheses
  
  p_maj <- ImageDimPlot(xenium_obj, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Majority")) # Fixed title
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterMajority.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_maj)
  dev.off()
  
  # --- Spatial Facet Plots ---
  coords <- GetTissueCoordinates(xenium_obj) 
  
  # Weighted Facet
  plot_data_wei <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_weighted), levels = celltype_order))
  p_facet_wei <- ggplot(plot_data_wei, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = facet_point_size) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterWeighted.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_wei)
  dev.off()
  
  # Majority Facet
  plot_data_maj <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_majority), levels = celltype_order))
  p_facet_maj <- ggplot(plot_data_maj, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = facet_point_size) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterMajority.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_maj)
  dev.off()
  
  # --- DotPlots ---
  existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium_obj)))
  
  # Weighted DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_weighted, levels = rev(celltype_order))
  p_dot_wei <- DotPlot(
    xenium_obj,
    features = existing_markers,
    assay = "Xenium",
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot_wei <- standardize_broad_dotplot(p_dot_wei) +
    ggtitle(paste(sample_name, "Markers (Weighted)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_wei)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.pdf")), p_dot_wei, device = grDevices::cairo_pdf, width = 10, height = 6)
  
  # Majority DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_majority, levels = rev(celltype_order))
  p_dot_maj <- DotPlot(
    xenium_obj,
    features = existing_markers,
    assay = "Xenium",
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot_maj <- standardize_broad_dotplot(p_dot_maj) +
    ggtitle(paste(sample_name, "Markers (Majority)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_maj)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.pdf")), p_dot_maj, device = grDevices::cairo_pdf, width = 10, height = 6)
  
  ## ----------------------------
  ## 6. Save & Return
  ## ----------------------------
  output_dir <- file.path(annotation_dir, "rds")
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(xenium_obj, file.path(output_dir, paste0(sample_name, "_", reference_name, "_annotated.rds")))
  message("Successfully annotated and saved: ", sample_name)
  
  return(xenium_obj)
}

# =========================================================
# Execution
# =========================================================
seu <- readRDS(input_file)
if (resolution_mode) {
  required_pilot_metadata <- c(
    "whole_tissue_cluster_res1.5",
    pilot_cluster_column,
    pilot_graph_column,
    "seurat_clusters"
  )
  missing_pilot_metadata <- setdiff(required_pilot_metadata, colnames(seu[[]]))
  if (length(missing_pilot_metadata)) {
    stop(
      "Resolution-", pilot_resolution_tag, " input lacks metadata: ",
      paste(missing_pilot_metadata, collapse = ", ")
    )
  }
  candidate_clusters <- as.character(seu[[pilot_cluster_column]][, 1])
  seurat_clusters <- as.character(seu$seurat_clusters)
  if (anyNA(candidate_clusters) || !identical(candidate_clusters, seurat_clusters)) {
    stop(
      "Resolution-specific input seurat_clusters does not exactly match ",
      pilot_cluster_column, " for every cell."
    )
  }
  if (anyDuplicated(Cells(seu))) {
    stop("Resolution-specific input contains duplicate cell IDs.")
  }
}
annotate_xenium_from_ref(
  xenium_obj = seu,
  sample_name = current_sample,
  reference_name = reference_name,
  reference_path = reference_path,
  annotation_dir = annotation_dir
)
