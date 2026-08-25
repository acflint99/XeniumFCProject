# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)

plot_dir <- here("outputs", "Aldinger_VZ_Plots")
rds_dir  <- here("outputs", "Aldinger_VZ_RDS")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(rds_dir))  dir.create(rds_dir, recursive = TRUE)


obj <- readRDS(here("outputs", "SingleCellRDS", "Aldinger_newClusters_newUMAPv2_5k.rds"))

target_clusters <- c("Glia", "GABA", "Purkinje", "OPC") 

# Subset the object
vz_obj <- subset(obj, idents = target_clusters)

# 4. STANDARD WORKFLOW
message("Normalizing and Finding Variable Features...")
vz_obj <- NormalizeData(vz_obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message("Scaling data...")
vz_obj <- ScaleData(vz_obj, features = VariableFeatures(vz_obj), verbose = TRUE)

message("Running PCA...")
vz_obj <- RunPCA(vz_obj, npcs = 30, verbose = FALSE)

# 5. BATCH CHECK (UNCORRECTED UMAP)
message("Generating uncorrected UMAP...")
vz_obj <- RunUMAP(vz_obj, reduction = "pca", dims = 1:30, 
               reduction.name = "umap_uncorrected", verbose = TRUE)

p1 <- DimPlot(vz_obj, reduction = "umap_uncorrected", group.by = "clusters_refined")

Cairo::CairoTIFF(
  filename = file.path(plot_dir, paste0("Aldinger_VZ_VZClusterUMAP.tif")),
  width = 10,
  height = 10,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()

vz_obj <- FindClusters(vz_obj, resolution = 0.8, verbose = TRUE)

p2 <- DimPlot(vz_obj, reduction = "umap_uncorrected", group.by = "RNA_snn_res.0.8")

Cairo::CairoTIFF(
  filename = file.path(plot_dir, paste0("Aldinger_VZ_RawCluster_UMAP_Res_0.8.tif")),
  width = 10,
  height = 10,
  units = "in",
  res = 600
)
print(p2)
grDevices::dev.off()

output_path <- file.path(rds_dir, "Aldinger_VZ_4126.rds")
saveRDS(vz_obj, output_path, compress = FALSE)


