#!/usr/bin/env Rscript

# Standardize the published Science reference labels. Existing outputs are
# protected because every downstream Science transfer depends on this object.

# Clear the environment
rm(list = ls())

# Load only the path dependency before the path-only dry-run exits.
suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/science_01_standardize_cell_types.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

# Define output directories
plot_path <- here("outputs", "references", "science", "plots")
RDS_path <- here("outputs", "references", "science", "rds")

input_path <- file.path(
  config$inputs$published_data_root,
  "ScienceBraunFC", "Luo", "Science_9-15pcw_res1.2_broad_label_clusters.rds"
)
output_paths <- c(
  file.path(plot_path, "ScienceUMAP_newClusters.tiff"),
  file.path(plot_path, "ScienceUMAP_newClusters_newUMAP50v1.tiff"),
  file.path(RDS_path, "Science_newClusters.rds"),
  file.path(RDS_path, "Science_newClusters_newUMAPv1.rds")
)

if (dry_run) {
  compact_dry_run(
    "Science reference standardization",
    inputs = input_path,
    outputs = output_paths
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("Science source object not found: ", input_path)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing Science reference outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the Science reference changes."
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

# Load the FC dataset
Science <- readRDS(input_path)

if (!"science_broad_label" %in% colnames(Science[[]])) {
  stop("Science source object lacks required metadata column 'science_broad_label'.")
}
if (anyDuplicated(colnames(Science))) stop("Science source object contains duplicate cell IDs.")
if (!"umap" %in% Reductions(Science)) stop("Science source object lacks the expected 'umap' reduction.")

Idents(Science) <- "science_broad_label"

# remove clusters----
Science <- subset(Science, idents = c("RBC", "Brainstem", "Cycling", "Other"), invert = TRUE)

expected_labels <- c(
  "Purkinje", "GABA", "Granule", "RL", "UBC", "Glia",
  "OPC", "Immune", "Endothelial", "Meninges", "VZ"
)

unexpected_labels <- setdiff(
  unique(as.character(Science$science_broad_label)),
  expected_labels
)

if (length(unexpected_labels) > 0) {
  stop(
    "Unexpected science_broad_label values: ",
    paste(unexpected_labels, collapse = ", ")
  )
}

Science$clusters_refined <- as.character(Science$science_broad_label)

if (anyNA(Science$clusters_refined) || any(!nzchar(as.character(Science$clusters_refined)))) {
  stop("Science standardized labels contain missing or blank values.")
}


# Apply celltype order from color_palette.R as factor levels
Science$clusters_refined <- factor(
  Science$clusters_refined, 
  levels = intersect(celltype_order, unique(Science$clusters_refined))
)

Idents(Science) <- "clusters_refined"

# plot UMAP again with new cluster labels and color_palette.R colors ----
p1 <- DimPlot(Science, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_newClusters.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p1)
dev.off()

# save Seurat object w/ new cluster labels----
saveRDS(Science, file = file.path(RDS_path, "Science_newClusters.rds"))

# redo PCA & UMAP----
# 1️⃣ Normalize data
Science_newUMAP <- NormalizeData(Science, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
Science_newUMAP <- FindVariableFeatures(Science_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
Science_newUMAP <- ScaleData(Science_newUMAP, features = rownames(Science_newUMAP))

# 4️⃣ Run PCA
Science_newUMAP <- RunPCA(Science_newUMAP, features = VariableFeatures(Science_newUMAP))

# 5️⃣ Find neighbors
Science_newUMAP <- FindNeighbors(Science_newUMAP, dims = 1:50)

# retain previous cluster identities & order
Science_newUMAP$clusters_refined <- factor(
  Science_newUMAP$clusters_refined,
  levels = intersect(celltype_order, unique(Science_newUMAP$clusters_refined))
)
Idents(Science_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
Science_newUMAP <- RunUMAP(
  Science_newUMAP,
  dims = 1:50,
  seed.use = config$runtime$random_seed
)

# 8️⃣ Plot UMAP with color_palette.R colors ----
p2 <- DimPlot(Science_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_newClusters_newUMAP50v1.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

# remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(Science_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  Science_newUMAP[[assay]]@scale.data <- matrix()
}

# save Seurat object w/ new cluster labels & new UMAP reductions----
saveRDS(Science_newUMAP, file = file.path(RDS_path, "Science_newClusters_newUMAPv1.rds"))
