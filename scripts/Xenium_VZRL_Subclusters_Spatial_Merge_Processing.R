# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

library(Seurat)
library(here)
library(dplyr)
library(ggplot2)
library(future)

source(here("scripts", "color_palette.R"))

# --- 1. Parallelization Setup ---
# Using 8 cores with a high memory limit for workers
plan("multisession", workers = 8)
options(future.globals.maxSize = 100 * 1024^3) 

# Ensure output directory exists
plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- 2. Load the CLEAN Merged Object ---
message("Loading the cleanly merged object...")
merged_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS", "XenAld_VZ_RL_QC_Subclusters_Spatial_Merged_4-22-26.rds")
obj <- readRDS(merged_path)

# Ensure the full object is "joined" before starting to prevent layer fragmentation
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])

# --- 3. Correct Pre-Processing Order (Full Object) ---
message("Normalizing and Finding Variable Features (Full Object)...")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

# --- 4. Sketching 100k Cells ---
message("Sketching 100k cells...")
obj <- SketchData(
  object = obj,
  ncells = 100000, 
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# --- 5. Reconstruct Clean Sketch Assay ---
DefaultAssay(obj) <- "Xenium"
sketch_cells <- colnames(obj[["sketch"]])

# Extract raw counts safely using Seurat v5 accessor
sketch_counts_matrix <- LayerData(obj, assay = "Xenium", layer = "counts")[, sketch_cells]

# Create a fresh Assay object to wipe the fragmented layer history
new_sketch_assay <- CreateAssay5Object(counts = sketch_counts_matrix)
obj[["sketch"]] <- new_sketch_assay

# --- 6. Pre-process the Clean Sketch ---
DefaultAssay(obj) <- "sketch"

# Explicitly grab sample IDs just for sketched cells to ensure perfect splitting
sketch_sample_ids <- obj@meta.data[sketch_cells, "sample_id"]

# Split by sample (should work perfectly on the clean assay)
obj[["sketch"]] <- split(obj[["sketch"]], f = sketch_sample_ids)

message("Re-processing clean sketch for Harmony...")
n_features_to_use <- min(2000, nrow(obj[["sketch"]]))

# Enforce correct normalization -> variable features -> scale order
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_features_to_use)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)


# --- 7. Harmony Integration ---
message("Running Harmony Integration...")
obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  group.by = "sample_id", 
  verbose = TRUE
)


# --- 8. Project to Full Dataset ---
message("Projecting integration and UMAP to full dataset...")
obj[["sketch"]] <- JoinLayers(obj[["sketch"]])

# Project the Harmony dimensions
obj <- ProjectIntegration(
  object = obj,
  sketched.assay = "sketch",
  assay = "Xenium",
  reduction = "integrated.harmony"
)

# --- FIX: Add return.model = TRUE so ProjectData can use it ---
obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = TRUE)

# Project the UMAP coordinates to the remaining unsketched cells
obj <- ProjectData(
  object = obj,
  sketched.assay = "sketch",
  assay = "Xenium",
  sketched.reduction = "integrated.harmony", 
  full.reduction = "integrated.harmony",     
  dims = 1:30                                
)


# --- FIX: Reset Default Assay back to the full dataset for plotting ---
DefaultAssay(obj) <- "Xenium"


# --- 9. Visualization ---
message("Generating Plots...")

# Enforce cluster ordering safely
if (exists("celltype_order") && "cluster_weighted" %in% colnames(obj@meta.data)) {
  obj$cluster_weighted <- factor(
    obj$cluster_weighted, 
    levels = intersect(celltype_order, unique(obj$cluster_weighted))
  )
}

if (exists("master_subcluster_order") && "comb_subcluster" %in% colnames(obj@meta.data)) {
  obj$comb_subcluster <- factor(
    obj$comb_subcluster, 
    levels = intersect(master_subcluster_order, unique(obj$comb_subcluster))
  )
}

# Generate the Plots
p1 <- DimPlot(obj, reduction = "umap.harmony", group.by = "sample_id", raster = TRUE) + 
  ggtitle("Integrated by Sample")
ggsave(file.path(plot_dir, "XenAld_VZ_RL_Subclusters_Spatial_Merged_SampleID_UMAP.png"), p1, width = 14, height = 14)

if ("cluster_weighted" %in% colnames(obj@meta.data) && exists("cluster_colors")) {
  p2 <- DimPlot(obj, reduction = "umap.harmony", group.by = "cluster_weighted", label = TRUE, raster = TRUE) + 
    scale_color_manual(values = cluster_colors) + 
    ggtitle("Integrated: Broad Clusters")
  ggsave(file.path(plot_dir, "XenAld_VZ_RL_Subclusters_Spatial_Merged_BroadClusters_UMAPs2.png"), p2, width = 14, height = 14)
}

if ("comb_subcluster" %in% colnames(obj@meta.data) && exists("subcluster_palette")) {
  p3 <- DimPlot(obj, 
                reduction = "umap.harmony", 
                group.by = "comb_subcluster", 
                label = TRUE, 
                label.size = 3.5,      
                repel = TRUE,          
                label.box = TRUE,      
                raster = TRUE, 
                alpha = 0.6,
                cols = subcluster_palette) + 
    ggtitle("Integrated: Subclusters") +
    theme(
      legend.text = element_text(size = 9), 
      legend.key.size = unit(0.5, "cm")     
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))
  ggsave(file.path(plot_dir, "XenAld_VZ_RL_Subclusters_Spatial_Merged_Subclusters_UMAP2.png"), p3, width = 14, height = 14)
}

# --- 10. Save ---
save_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS", "XenAld_VZ_RL_QC_Subclusters_Spatial_Merged_Integrated2_4-23-26.rds")
message("Saving final integrated object to: ", save_path)

# Return to sequential plan for safer disk writing
plan("sequential")
saveRDS(obj, save_path, compress = FALSE)

message("Done! Ready for Trajectory Analysis.")