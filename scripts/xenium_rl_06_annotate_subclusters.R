#!/usr/bin/env Rscript

# Apply the reviewed RL subcluster labels and create summary figures.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_rl_06_annotate_subclusters.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "rl", "06_subclusters", "plots")
table_path <- file.path(output_root, "xenium", "rl", "06_subclusters", "tables")
rds_dir <- file.path(output_root, "xenium", "rl", "06_subclusters", "rds")
merged_path <- file.path(output_root, "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC_Res1.5.rds")
output_path <- file.path(rds_dir, "Xenium_RL_subclusters_Res1.5.rds")
label_manifest_path <- file.path(table_path, "Xenium_RL_subcluster_labels.csv")
expected_outputs <- c(
  output_path,
  label_manifest_path,
  file.path(table_path, "XenAld_RL_Subcluster_QC_Summary.csv"),
  file.path(table_path, "XenAld_RL_Subcluster_Markers.csv"),
  file.path(plot_path, "XenAld_RL_Subcluster_QC_Violins.pdf"),
  file.path(plot_path, "XenAld_RL_Subcluster_UMAP.tif"),
  file.path(plot_path, "XenAld_RL_Subcluster_rmBG,Ep,PC_UMAP.tif"),
  file.path(plot_path, "XenAld_RL_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "XenAld_RL_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "XenAld_RL_Subcluster_Top5_Heatmap.tif"),
  file.path(plot_path, "XenAld_RL_Subcluster_Top5_Heatmap.pdf")
)

if (dry_run) {
  cat("RL post-QC input:", merged_path, "\n")
  cat("RL post-QC input exists:", file.exists(merged_path), "\n")
  cat("Annotation cluster column: Xenium_snn_res.0.5\n")
  cat("Configured workers:", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"), "\n")
  write.table(
    data.frame(output = expected_outputs, exists = file.exists(expected_outputs)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(merged_path)) stop("RL post-QC input not found: ", merged_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL subcluster outputs:\n- ",
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
required_metadata <- c("orig.ident", "consensus_label", "PCW", "Xenium_snn_res.0.5")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("RL post-QC input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!"umap_clean" %in% Reductions(obj)) stop("RL post-QC input lacks the umap_clean reduction.")
if (!"Xenium" %in% Assays(obj)) stop("RL post-QC input lacks the Xenium assay.")
if (ncol(obj) == 0L) stop("RL post-QC input contains zero cells.")

# # 2. SET THE IDENTITY
# # Ensure we are looking at the resolution you liked best (e.g., 0.3)
# Idents(obj) <- "Xenium_snn_res.0.3"
# 
# check_mem("STARTING FINDALLMARKERS")
# 
# # 3. RUN FINDALLMARKERS (The Optimized Way)
# message(Sys.time(), ": Starting Marker Identification...")
# 
# all_markers <- FindAllMarkers(
#   obj,
#   only.pos = TRUE,          # Only look for upregulated genes (standard for cell types)
#   min.pct = 0.25,           # Gene must be in 25% of the cluster
#   logfc.threshold = 0.25,   # Minimum 1.28x fold change
#   test.use = "wilcox",      # Standard fast test
#   
#   # --- THE SPEED TRICK ---
#   # Downsampling to 1000 cells per cluster gives 99% the same results
#   # but runs 5x faster on massive Xenium datasets.
#   max.cells.per.ident = 1000
# )
# 
# check_mem("MARKERS COMPLETE")
# 
# # 4. FILTER & SAVE TOP 10
# top10_markers <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# write.csv(top10_markers,
#           file.path(table_path, "RL_RawSubluster_postQC_top10_Markers_Res0.3.csv"),
#           row.names = FALSE)
# 
# message("Marker analysis complete! Results saved to CSV.")


####input markers into Gemini####



Idents(obj) <- "Xenium_snn_res.0.5"

# 1. Create a named vector for the mapping
# The names (0, 1, 2...) must match your current cluster IDs exactly
new_labels <- c(
  "0" = "Prolif GCPs",
  "1" = "Mature UBCs",
  "2" = "Immature UBCs",
  "3" = "Migrating GCs",
  "4" = "Migrating GCs",
  "5" = "Maturing GCPs",
  "6" = "Differentiating GCs",
  "7" = "Prolif GCPs",
  "8" = "Maturing PCs",
  "9" = "RL SVZ",
  "10" = "BG",
  "11" = "Mature GCs",
  "12" = "Maturing GCPs",
  "13" = "RL VZ",
  "14" = "Prolif GCPs",
  "15" = "Ependymal Cells",
  "16" = "Intermediate Progenitors"
)

current_clusters <- levels(Idents(obj))
missing_label_mappings <- setdiff(current_clusters, names(new_labels))
unused_label_mappings <- setdiff(names(new_labels), current_clusters)
if (length(missing_label_mappings) || length(unused_label_mappings)) {
  stop(
    "RL label mapping does not exactly match Xenium_snn_res.0.5. Missing: ",
    if (length(missing_label_mappings)) paste(missing_label_mappings, collapse = ", ") else "none",
    "; unused: ",
    if (length(unused_label_mappings)) paste(unused_label_mappings, collapse = ", ") else "none"
  )
}
if (anyNA(new_labels) || any(!nzchar(new_labels))) stop("RL label mapping contains missing or blank labels.")

obj <- RenameIdents(obj, new_labels)

# Final metadata assignment
obj$RL_subcluster <- Idents(obj)

write.csv(
  data.frame(
    cluster_id = names(new_labels),
    RL_subcluster = unname(new_labels),
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
obj$RL_subcluster <- factor(Idents(obj), levels = rl_subcluster_order)
Idents(obj) <- "RL_subcluster"

DefaultAssay(obj) <- "Xenium"

# Export Merged Cluster QC Statistics
merged_qc_stats <- obj@meta.data %>%
  group_by(RL_subcluster) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(merged_qc_stats,
          file.path(table_path, "XenAld_RL_Subcluster_QC_Summary.csv"),
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "XenAld_RL_Subcluster_QC_Violins.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)
  
  p <- VlnPlot(
    obj,
    features = feat,
    group.by = "RL_subcluster",
    cols = rl_palette,
    pt.size = 0      # Remove points for speed and smaller file size
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"), # Explicitly style Y
      legend.position = "none"
    ) +
    ylab(feat) +
    ggtitle(paste("Xenium RL Subcluster QC-", feat))
  
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
             cols = rl_palette) + # Uses your defined palette
  ggtitle("Xenium Merged RL Subclusters") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_RL_Subcluster_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()

# 4. UMAP WITH "Other" Clusters removed
# Define the clusters you want to REMOVE
clusters_to_remove <- c("BG", "Ependymal Cells", "Maturing PCs")  #"RG Progenitors",

# Create a subset of the object
# ! means 'not', %in% checks if the identity is in your removal list
obj_subset <- subset(obj, idents = clusters_to_remove, invert = TRUE)

# Now run your plotting code on 'obj_subset' instead of 'obj'
p2 <- DimPlot(obj_subset,
              reduction = "umap_clean",
              label = TRUE,
              label.size = 4,
              repel = TRUE,
              cols = rl_palette) +
  ggtitle("Xenium Merged RL Subclusters (Filtered)") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_RL_Subcluster_rmBG,Ep,PC_UMAP.tif")),
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
obj$RL_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(rl_subcluster_order))
Idents(obj) <- "RL_subcluster_rev"

# 3. Create the Plot
# Adding 'assay = "Xenium"' is the safest way to avoid multi-assay conflicts
p1 <- DotPlot(obj,
              features = rl_markers,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged RL Subcluster Markers")

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_RL_SubclusterMarker_DotPlot.tif")),
  width = 12,
  height = 8,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "XenAld_RL_SubclusterMarker_DotPlot.pdf"),
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

write.csv(heatmap_markers,
          file.path(table_path, "XenAld_RL_Subcluster_Markers.csv"),
          row.names = FALSE)

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
  group.by = "RL_subcluster",
  group.colors = rl_palette,
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
  ggtitle("Top 5 Markers per RL Subcluster")

# 4. SAVE AS PDF
ggplot2::ggsave(
  file.path(plot_path, paste0("XenAld_RL_Subcluster_Top5_Heatmap.pdf")),
  p4, width = 14, height = 12, device = grDevices::cairo_pdf, useDingbats = FALSE
)

# 5. Save the Heatmap
Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_RL_Subcluster_Top5_Heatmap.tif")),
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
