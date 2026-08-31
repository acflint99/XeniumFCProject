#!/usr/bin/env Rscript

# Restrict the standardized Science reference to the Xenium gene panel.

# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# Load only the path dependency before the path-only dry-run exits.
library(here)

source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/science_02_subset_gene_panel.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

# Define output directories
plot_path <- here("outputs", "references", "science", "plots")
RDS_path <- here("outputs", "references", "science", "rds")

input_path <- file.path(RDS_path, "Science_newClusters_newUMAPv1.rds")
panel_path <- here("inputs", "xenium_5k_genes.rds")
output_paths <- c(
  file.path(plot_path, "ScienceUMAP_FC_newClusters_5k.tiff"),
  file.path(plot_path, "ScienceUMAP_newClusters_5k_newUMAP50v2.tiff"),
  file.path(plot_path, "ScienceDotPlot_newClusters_5k_newUMAP50v2_markers.tiff"),
  file.path(plot_path, "ScienceDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"),
  file.path(RDS_path, "Science_newClusters_UMAPv1_5k.rds"),
  file.path(RDS_path, "Science_newClusters_newUMAPv2_5k.rds")
)

if (dry_run) {
  cat("Standardized Science input:", input_path, "\n")
  cat("Input exists:", file.exists(input_path), "\n")
  cat("Xenium panel:", panel_path, "\n")
  cat("Panel exists:", file.exists(panel_path), "\n")
  write.table(
    data.frame(output = output_paths, exists = file.exists(output_paths)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = 0L)
}

missing_inputs <- c(input_path, panel_path)[!file.exists(c(input_path, panel_path))]
if (length(missing_inputs)) stop("Science panel stage is missing input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing Science panel outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the standardized Science reference."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(Cairo)
})
source(here("scripts", "color_palette.R"))

# Ensure directories exist before saving any files
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
if (!dir.exists(RDS_path)) dir.create(RDS_path, recursive = TRUE)

# Load the FC dataset ----
Science <- readRDS(input_path)

if (!"clusters_refined" %in% colnames(Science[[]])) {
  stop("Standardized Science input lacks 'clusters_refined'.")
}
if (anyDuplicated(colnames(Science))) stop("Standardized Science input contains duplicate cell IDs.")
if (!all(c("Glia", "RL") %in% unique(as.character(Science$clusters_refined)))) {
  stop("Standardized Science input must contain both Glia and RL labels.")
}

xenium_genes <- readRDS(panel_path)

genes_present <- intersect(xenium_genes, rownames(Science))
if (!length(genes_present)) stop("Science reference and Xenium panel have no genes in common.")

ScienceSubset <- subset(
  Science,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

# Apply celltype order from color_palette.R as factor levels
ScienceSubset$clusters_refined <- factor(
  ScienceSubset$clusters_refined,
  levels = intersect(celltype_order, unique(ScienceSubset$clusters_refined))
)

# check UMAP
p <- DimPlot(ScienceSubset, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_FC_newClusters_5k.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p)
dev.off()

# export new Seurat object as .rds ----
saveRDS(ScienceSubset, file = file.path(RDS_path, "Science_newClusters_UMAPv1_5k.rds"))


## redo PCA & UMAP ----
# 1️⃣ Normalize data
ScienceSubset_newUMAP <- NormalizeData(ScienceSubset, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
ScienceSubset_newUMAP <- FindVariableFeatures(ScienceSubset_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
ScienceSubset_newUMAP <- ScaleData(ScienceSubset_newUMAP, features = rownames(ScienceSubset_newUMAP))

# 4️⃣ Run PCA
ScienceSubset_newUMAP <- RunPCA(ScienceSubset_newUMAP, features = VariableFeatures(ScienceSubset_newUMAP))

# 5️⃣ Find neighbors
ScienceSubset_newUMAP <- FindNeighbors(ScienceSubset_newUMAP, dims = 1:50)

# retain previous cluster identities & order
ScienceSubset_newUMAP$clusters_refined <- factor(
  ScienceSubset_newUMAP$clusters_refined,
  levels = intersect(celltype_order, unique(ScienceSubset_newUMAP$clusters_refined))
)
Idents(ScienceSubset_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
ScienceSubset_newUMAP <- RunUMAP(ScienceSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP with color_palette.R colors
p2 <- DimPlot(ScienceSubset_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_newClusters_5k_newUMAP50v2.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

# remove scale.data for all assays to reduce file size ----
# Get assay names
assay_names <- names(ScienceSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  ScienceSubset_newUMAP[[assay]]@scale.data <- matrix()
}

# export new Seurat object as .rds ----
saveRDS(ScienceSubset_newUMAP, file = file.path(RDS_path, "Science_newClusters_newUMAPv2_5k.rds"))

# Reverse the celltype order for the DotPlot so the y-axis renders top-to-bottom
Idents(ScienceSubset_newUMAP) <- factor(
  ScienceSubset_newUMAP$clusters_refined, 
  levels = rev(intersect(celltype_order, unique(ScienceSubset_newUMAP$clusters_refined)))
)

dotplot_assay <- DefaultAssay(ScienceSubset_newUMAP)
existing_markers <- lapply(
  markers,
  function(features) intersect(features, rownames(ScienceSubset_newUMAP[[dotplot_assay]]))
)
existing_markers <- existing_markers[lengths(existing_markers) > 0L]
if (!length(existing_markers)) {
  stop("None of the configured broad-cell markers are present in the Science panel object.")
}

p3 <- DotPlot(
  ScienceSubset_newUMAP,
  features = existing_markers,
  assay = dotplot_assay,
  col.min = broad_dotplot_col_min,
  col.max = broad_dotplot_col_max,
  dot.min = broad_dotplot_dot_min / 100,
  dot.scale = broad_dotplot_dot_scale,
  scale.min = broad_dotplot_dot_min,
  scale.max = broad_dotplot_dot_max
)
p3 <- standardize_broad_dotplot(p3) +
  ggtitle("High-Specificity Marker Expression by Cluster")

CairoTIFF(filename = file.path(plot_path, "ScienceDotPlot_newClusters_5k_newUMAP50v2_markers.tiff"), 
          width = 10, height = 6, units = "in", res = 600)
print(p3)
dev.off()

# 2. Save Dotplot as PDF
ggplot2::ggsave(
  filename = file.path(plot_path, "ScienceDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"), 
  plot = p3, 
  width = 10, 
  height = 6,
  device = grDevices::cairo_pdf
) 
