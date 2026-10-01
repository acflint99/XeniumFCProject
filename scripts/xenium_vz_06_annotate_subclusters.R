#!/usr/bin/env Rscript

# Review post-QC VZ clusters or explicitly apply the reviewed labels.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--review-only", "--apply-reviewed-labels", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_vz_06_annotate_subclusters.R ",
    "(--review-only|--apply-reviewed-labels) [--dry-run|--overwrite]"
  )
}
review_only <- "--review-only" %in% args
apply_reviewed_labels <- "--apply-reviewed-labels" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (sum(c(review_only, apply_reviewed_labels)) != 1L) {
  stop("Choose exactly one mode: --review-only or --apply-reviewed-labels.")
}

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "vz", "06_subclusters", "plots")
table_path <- file.path(output_root, "xenium", "vz", "06_subclusters", "tables")
rds_dir <- file.path(output_root, "xenium", "vz", "06_subclusters", "rds")
merged_path <- file.path(output_root, "xenium", "vz", "05_post_qc", "rds", "Xenium_VZ_postQC.rds")
output_path <- file.path(rds_dir, "Xenium_VZ_subclusters.rds")
label_manifest_path <- file.path(table_path, "Xenium_VZ_subcluster_labels.csv")
review_marker_path <- file.path(
  table_path, "VZ_RawSubcluster_postQC_top10_Markers_Res0.5.csv"
)
review_manifest_path <- file.path(table_path, "Xenium_VZ_subcluster_review_manifest.csv")
review_umap_path <- file.path(
  plot_path, "Xenium_VZ_RawSubcluster_postQC_UMAP_Res0.5.tif"
)
review_outputs <- c(review_marker_path, review_manifest_path, review_umap_path)
expected_outputs <- c(
  output_path,
  label_manifest_path,
  file.path(table_path, "Xenium_VZ_Subcluster_QC_Summary.csv"),
  file.path(plot_path, "Xenium_VZ_Subcluster_QC_Violins.pdf"),
  file.path(plot_path, "Xenium_VZ_Subcluster_UMAP.tif"),
  file.path(plot_path, "Xenium_VZ_Subcluster_rmGCP,eCN,Cyc_UMAP.tif"),
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_VZ_Subcluster_Top5_Heatmap.tif"),
  file.path(plot_path, "Xenium_VZ_Subcluster_Top5_Heatmap.pdf")
)

file_is_nonempty <- function(path) {
  info <- suppressWarnings(file.info(path))
  file.exists(path) && !is.na(info$size) && info$size > 0
}

review_manifest_matches_input_file <- function() {
  if (!file.exists(merged_path) || !file.exists(review_manifest_path)) return(FALSE)
  review <- tryCatch(
    read.csv(review_manifest_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(review) || nrow(review) != 1L ||
      !all(c("input_size", "input_mtime") %in% names(review))) return(FALSE)
  input_info <- file.info(merged_path)
  isTRUE(all.equal(as.numeric(review$input_size), as.numeric(input_info$size))) &&
    isTRUE(all.equal(as.numeric(review$input_mtime), as.numeric(input_info$mtime)))
}

if (dry_run) {
  compact_dry_run(
    if (review_only) "VZ post-QC subcluster review" else "VZ reviewed-label application",
    inputs = if (review_only) merged_path else c(merged_path, review_outputs),
    outputs = if (review_only) review_outputs else expected_outputs,
    checks = if (apply_reviewed_labels) {
      c(
        review_outputs_nonempty = all(vapply(review_outputs, file_is_nonempty, logical(1))),
        review_manifest_matches_input = review_manifest_matches_input_file()
      )
    } else logical()
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(merged_path)) stop("VZ post-QC input not found: ", merged_path)
if (apply_reviewed_labels) {
  missing_review_outputs <- review_outputs[!file.exists(review_outputs)]
  if (length(missing_review_outputs)) {
    stop(
      "Run --review-only and inspect its outputs before applying VZ labels. Missing:\n- ",
      paste(missing_review_outputs, collapse = "\n- ")
    )
  }
  empty_review_outputs <- review_outputs[
    !vapply(review_outputs, file_is_nonempty, logical(1))
  ]
  if (length(empty_review_outputs)) {
    stop("VZ review outputs are empty:\n- ", paste(empty_review_outputs, collapse = "\n- "))
  }
}
mode_outputs <- if (review_only) review_outputs else expected_outputs
existing_outputs <- mode_outputs[file.exists(mode_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ ",
    if (review_only) "review" else "subcluster", " outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(patchwork)
})
source(here("scripts", "color_palette.R"))

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(config$runtime$random_seed)

with_future_workers <- function(expr) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  if (workers > 1L) plan(multisession, workers = workers)
  force(expr)
}

check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(merged_path)
required_metadata <- c(
  "orig.ident", "consensus_label", "regional_consensus_label",
  "regional_consensus_label_pooled", "regional_unknown_label",
  "PCW", "Xenium_snn_res.0.5"
)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("VZ post-QC input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
  stop("VZ post-QC input has invalid sample-qualified Unknown labels.")
}
if (!"umap_clean" %in% Reductions(obj)) stop("VZ post-QC input lacks the umap_clean reduction.")
if (!"Xenium" %in% Assays(obj)) stop("VZ post-QC input lacks the Xenium assay.")
if (ncol(obj) == 0L) stop("VZ post-QC input contains zero cells.")

cluster_column <- "Xenium_snn_res.0.5"
Idents(obj) <- cluster_column
current_clusters <- levels(Idents(obj))

if (review_only) {
  DefaultAssay(obj) <- "Xenium"
  check_mem("STARTING POST-QC SUBCLUSTER REVIEW")
  message(Sys.time(), ": identifying raw VZ cluster markers for review...")
  all_markers <- with_future_workers(
    FindAllMarkers(
      obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident = 1000,
      random.seed = config$runtime$random_seed
    )
  )
  if (!nrow(all_markers)) stop("No VZ post-QC markers were returned for review.")
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  write.csv(top10_markers, review_marker_path, row.names = FALSE)

  p_review <- DimPlot(
    obj,
    reduction = "umap_clean",
    group.by = cluster_column,
    label = TRUE,
    label.size = 4,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("VZ post-QC raw clusters - resolution 0.5") +
    theme_classic()
  Cairo::CairoTIFF(
    filename = review_umap_path,
    width = 12,
    height = 9,
    units = "in",
    res = 600
  )
  print(p_review)
  grDevices::dev.off()

  input_info <- file.info(merged_path)
  write.csv(
    data.frame(
      input_rds = normalizePath(merged_path, mustWork = TRUE),
      input_size = as.numeric(input_info$size),
      input_mtime = as.numeric(input_info$mtime),
      cells = ncol(obj),
      cluster_column = cluster_column,
      clusters = paste(current_clusters, collapse = "|"),
      random_seed = config$runtime$random_seed,
      min_pct = 0.25,
      logfc_threshold = 0.25,
      max_cells_per_ident = 1000,
      stringsAsFactors = FALSE
    ),
    review_manifest_path,
    row.names = FALSE
  )
  message("VZ post-QC review outputs are ready. No subcluster labels were applied.")
  quit(save = "no", status = 0L)
}

review <- read.csv(review_manifest_path, stringsAsFactors = FALSE)
required_review_columns <- c(
  "input_rds", "input_size", "input_mtime", "cells", "cluster_column", "clusters",
  "random_seed"
)
if (nrow(review) != 1L || !all(required_review_columns %in% names(review))) {
  stop("VZ review manifest is invalid; rerun --review-only.")
}
input_info <- file.info(merged_path)
review_is_current <-
  identical(as.character(review$input_rds), normalizePath(merged_path, mustWork = TRUE)) &&
  isTRUE(all.equal(as.numeric(review$input_size), as.numeric(input_info$size))) &&
  isTRUE(all.equal(as.numeric(review$input_mtime), as.numeric(input_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj))) &&
  identical(as.character(review$cluster_column), cluster_column) &&
  identical(as.character(review$clusters), paste(current_clusters, collapse = "|")) &&
  identical(as.integer(review$random_seed), as.integer(config$runtime$random_seed))
if (!review_is_current) {
  stop("VZ post-QC input or clusters changed after review; rerun --review-only.")
}

# 1. Create a named vector for the mapping
# The names (0, 1, 2...) must match your current cluster IDs exactly
new_labels <- c(
  "0"  = "Early-born PCs",
  "1"  = "Maturing PCs",
  "2"  = "GABA Progenitors",
  "3"  = "BG",
  "4"  = "Cycling Cells",
  "5"  = "RG Progenitors",
  "6"  = "Astrocytes/Ependyma",
  "7"  = "eCN",
  "8"  = "Patterning PCs",
  "9"  = "iCN",
  "10" = "VZPs",
  "11" = "Golgi Cells",
  "12" = "OPCs",
  "13" = "Late-born PCs",
  "14" = "Late-born PCs",
  "15" = "GCPs",
  "16" = "MLIs",
  "17" = "eCN",
  "18" = "Maturing PCs"

)

missing_label_mappings <- setdiff(current_clusters, names(new_labels))
unused_label_mappings <- setdiff(names(new_labels), current_clusters)
if (length(missing_label_mappings) || length(unused_label_mappings)) {
  stop(
    "VZ label mapping does not exactly match Xenium_snn_res.0.5. Missing: ",
    if (length(missing_label_mappings)) paste(missing_label_mappings, collapse = ", ") else "none",
    "; unused: ",
    if (length(unused_label_mappings)) paste(unused_label_mappings, collapse = ", ") else "none"
  )
}
if (anyNA(new_labels) || any(!nzchar(new_labels))) stop("VZ label mapping contains missing or blank labels.")

obj <- RenameIdents(obj, new_labels)

# Final metadata assignment
obj$VZ_subcluster <- Idents(obj)

write.csv(
  data.frame(
    cluster_id = names(new_labels),
    VZ_subcluster = unname(new_labels),
    source_column = "Xenium_snn_res.0.5",
    input_rds = merged_path,
    output_rds = output_path,
    stringsAsFactors = FALSE
  ),
  label_manifest_path,
  row.names = FALSE
)
saveRDS(obj, output_path, compress = FALSE)

# 3. APPLY STANDARDIZED ORDERING
# We use 'subcluster_order' from your color_palette.R script
obj$VZ_subcluster <- factor(Idents(obj), levels = vz_subcluster_order)
Idents(obj) <- "VZ_subcluster"

DefaultAssay(obj) <- "Xenium"

# Export Merged Cluster QC Statistics
merged_qc_stats <- obj@meta.data %>%
  group_by(VZ_subcluster) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(merged_qc_stats,
          file.path(table_path, "Xenium_VZ_Subcluster_QC_Summary.csv"),
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "Xenium_VZ_Subcluster_QC_Violins.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)

  p <- VlnPlot(
    obj,
    features = feat,
    group.by = "VZ_subcluster",
    cols = vz_palette,
    pt.size = 0      # Remove points for speed and smaller file size
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"), # Explicitly style Y
      legend.position = "none"
    ) +
    ylab(feat) +
    ggtitle(paste("Xenium VZ Subcluster QC-", feat))

  # Printing the plot inside the loop sends it to the PDF device
  print(p)
}

# 5. Close the device to finalize the file
dev.off()

message("Individual QC plots saved to: ", pdf_path)

# 4. UMAP WITH CONSISTENT COLORS
p <- DimPlot(obj,
             reduction = "umap_clean",
             label = TRUE,
             label.size = 4,
             repel = TRUE,
             cols = vz_palette) + # Uses your defined palette
  ggtitle("Xenium Merged VZ Subclusters") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_Subcluster_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()

# 4. UMAP WITH "Other" Clusters removed
# Define the clusters you want to REMOVE
clusters_to_remove <- c("GCPs", "eCN", "Cycling Cells")

# Create a subset of the object
# ! means 'not', %in% checks if the identity is in your removal list
obj_subset <- subset(obj, idents = clusters_to_remove, invert = TRUE)

# Now run your plotting code on 'obj_subset' instead of 'obj'
p2 <- DimPlot(obj_subset,
              reduction = "umap_clean",
              label = TRUE,
              label.size = 4,
              repel = TRUE,
              cols = vz_palette) +
  ggtitle("Xenium Merged VZ Subclusters (Filtered)") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_Subcluster_rmGCP,eCN,Cyc_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p2)
grDevices::dev.off()

rm(obj_subset)

# 5. DOTPLOT WITH REVERSED ORDER
# Create a temporary reversed factor for the Y-axis
# (Reverse the levels of the factor to get VZP at the top)
obj$VZ_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(vz_subcluster_order))
Idents(obj) <- "VZ_subcluster_rev"

# 3. Create the Plot
# Adding 'assay = "Xenium"' is the safest way to avoid multi-assay conflicts
p1 <- DotPlot(obj,
              features = vz_markers,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged VZ Subcluster Markers")

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_SubclusterMarker_DotPlot.tif")),
  width = 12,
  height = 8,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.pdf"),
  p1, device = grDevices::cairo_pdf, width = 12, height = 8
)

# ----------------------------------------------------------------
# 6. GENERATE TOP 5 MARKERS HEATMAP
# ----------------------------------------------------------------
check_mem("STARTING HEATMAP GENERATION")

# 1. Identify Markers (if not already in environment from earlier)
# We use a lower max.cells.per.ident to speed up the heatmap calculation
message("Finding top 5 markers for heatmap...")
heatmap_markers <- with_future_workers(
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    random.seed = config$runtime$random_seed
  )
)

# 2. Extract top 5 per cluster
top5_markers <- heatmap_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

top5_markers <- rev(top5_markers)

# 3. Scale the data for the specific markers only
# (Necessary for DoHeatmap to show relative expression)
# Just run the plot on a temporary scaled version:
temp_obj <- ScaleData(obj, features = top5_markers, verbose = FALSE)

# 4. Create Heatmap
# We downsample the plot to 100 cells per group so the labels are readable
# 1. Prepare the plot without the default lines
p4 <- DoHeatmap(
  subset(temp_obj, downsample = 100),
  features = top5_markers,
  group.by = "VZ_subcluster",
  group.colors = vz_palette,
  size = 4,
  angle = 45,
  draw.lines = TRUE,
  raster = FALSE
) +
  scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
  guides(color = "none") +
  theme(
    axis.text.y = element_text(size = 6, face = "italic"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Top 5 Markers per VZ Subcluster")

# 4. SAVE AS PDF
ggplot2::ggsave(
  file.path(plot_path, paste0("Xenium_VZ_Subcluster_Top5_Heatmap.pdf")),
  p4, width = 14, height = 12, device = grDevices::cairo_pdf, useDingbats = FALSE
)

# 5. Save the Heatmap
Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_Subcluster_Top5_Heatmap.tif")),
  width = 14,
  height = 12,
  units = "in",
  res = 600
)
print(p4)
grDevices::dev.off()

check_mem("HEATMAP COMPLETE")

plan("sequential")

check_mem("PIPELINE COMPLETE - FILE SAVED")
