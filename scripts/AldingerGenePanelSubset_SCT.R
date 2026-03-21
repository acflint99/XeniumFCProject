# Clear the environment
rm(list = ls())

# Load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(future) # New: for parallelization

source(here("scripts", "color_palette.R")) # Adjust path as necessary

# --- HPC OPTIMIZATION SETUP ---
# With 440GB RAM, we can be aggressive. 
# Using 12-16 workers is the "sweet spot" to avoid overhead while maximizing speed.
plan("multisession", workers = 12) 

# Set a 50GB memory limit for passing objects between cores
options(future.globals.maxSize = 50 * 1024^3) 

# --- Load the FC dataset ---
Aldinger <- readRDS(here("outputs", "SingleCellRDS", "Aldinger_newClusters_SCT_newUMAPv1.rds"))
xenium_genes <- readRDS(here("inputs", "xenium_5k_genes.rds"))
genes_present <- intersect(xenium_genes, rownames(Aldinger))

AldingerSubset <- subset(Aldinger, features = genes_present)

## redo PCA & UMAP using SCTransform (Parallelized) ----
AldingerSubset_newUMAP <- SCTransform(
  AldingerSubset, 
  assay = "RNA", 
  new.assay.name = "SCT", 
  verbose = TRUE
)

# Return to sequential processing to free up resources for the rest of the script
plan("sequential")

# Run Dimensionality Reduction
AldingerSubset_newUMAP <- RunPCA(AldingerSubset_newUMAP, verbose = FALSE)
AldingerSubset_newUMAP <- FindNeighbors(AldingerSubset_newUMAP, dims = 1:50)
AldingerSubset_newUMAP <- RunUMAP(AldingerSubset_newUMAP, dims = 1:50)

# Lock the cell type order based on your master palette
AldingerSubset_newUMAP$clusters_refined <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = celltype_order
)

# Plot UMAP
# Validate colors first (using your custom function)
validate_palette(unique(AldingerSubset_newUMAP$clusters_refined))

p2 <- DimPlot(AldingerSubset_newUMAP, reduction = "umap", group.by = "clusters_refined", cols = cluster_colors)
p2
ggsave(here("outputs", "SingleCellPlots", "AldingerUMAP_newClusters_5k_SCT_newUMAP50v2.pdf"), 
       plot = p2, width = 7, height = 6)


# --- check key cell type markers ---
markers <- list(
  "RL" = c("MKI67", "LTBP1", "OTX2"),
  "UBC" = c("EOMES"),
  "Granule" = c("ATOH1", "PAX6", "NEUROD1", "RELN"),
  #"VZ" = c("PRDM13", "DLL1", "ASCL1"),
  "Purkinje" = c("FOXP2", "CALB1", "DAB1"),
  "GABA" = c("PAX2", "GAD1", "GAD2"),
  "Glia" = c("SOX9", "ADCY2"),
  "OPC" = c("PDGFRA", "OLIG1"),
  "Meninges" = c("FOXC1", "SLC7A11"),
  "Endothelial" = c("CLDN5", "PECAM1"),
  "Immune" = c("P2RY12")
)


# Set Idents and Order
Idents(AldingerSubset_newUMAP) <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = rev(celltype_order)
)

p3 <- DotPlot(
  AldingerSubset_newUMAP,
  features = markers,
  assay = "SCT" # Ensure we use the SCT assay for the plot
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

p3
ggsave(here("outputs", "SingleCellPlots", "AldingerDotPlot_newclusters_5k_SCT_newUMAP50v2_markers.pdf"), 
       plot = p3, width = 14, height = 6)

#### --- Proportion Plot: clusters_refined by PCW --- ####

# 1. Calculate proportions
prop_data <- AldingerSubset_newUMAP@meta.data %>%
  group_by(age, clusters_refined) %>%
  tally() %>%
  group_by(age) %>%
  mutate(proportion = n / sum(n))

# 2. Generate the Stacked Bar Plot
p4 <- ggplot(prop_data, aes(x = age, y = proportion, fill = clusters_refined)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(
    title = "Cluster Proportions Across Fetal Development",
    x = "Post-Conceptional Week (PCW)",
    y = "Proportion of Cells",
    fill = "Cell Type"
  ) +
  # Apply your master palette here:
  scale_fill_manual(values = cluster_colors) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Show and Save
p4
ggsave(here("outputs", "SingleCellPlots", "AldingerBarPlot_newClusters_5k_SCT_newUMAP50v2_prop.pdf"), 
       plot = p4, width = 8, height = 6)

# Clean up scale.data to save space
assay_names <- names(AldingerSubset_newUMAP@assays)
for (assay in assay_names) {
  AldingerSubset_newUMAP[[assay]]@scale.data <- matrix()
}

# export final Seurat object
saveRDS(AldingerSubset_newUMAP, file = here("outputs", "SingleCellRDS", "Aldinger_newClusters_5k_SCT_newUMAPv2.rds"))
