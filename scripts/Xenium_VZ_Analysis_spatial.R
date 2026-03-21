# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# Load your new palette and order
source(here("scripts", "color_palette.R"))

check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 1. PARALLELIZATION (adjust to session cores and memory)
plan("multisession", workers = 7) 
options(future.globals.maxSize = 50 * 1024^3)

# UPDATED: Load the new Banksy Integrated object
merged_path <- here("outputs", "XenAld_VZ_Banksy_RDS", "Xenium_VZ_Banksy_Integrated_lambda0.1_31826.rds")
obj <- readRDS(merged_path)

# Ensure layers are joined for FindAllMarkers
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])

# 2. SET THE IDENTITY
# UPDATED: Use the correct Banksy k50 resolution metadata column
Idents(obj) <- "clust_M0_lam0.1_k50_res1"

check_mem("STARTING FINDALLMARKERS")

# # 3. RUN FINDALLMARKERS
# message(Sys.time(), ": Starting Marker Identification...")
# 
# all_markers <- FindAllMarkers(
#   obj,
#   only.pos = TRUE,          
#   min.pct = 0.25,           
#   logfc.threshold = 0.25,   
#   test.use = "wilcox",      
#   max.cells.per.ident = 1000 # Critical speed/memory trick
# )
# 
# check_mem("MARKERS COMPLETE")
# 
# # 4. FILTER & SAVE TOP 10
# top10_markers <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# # UPDATED: Save to the Banksy output folder
# write.csv(top10_markers,
#           here("outputs", "XenAld_VZ_Banksy_Tables", "Banksy_VZ_top10_Markers_Res1.csv"),
#           row.names = FALSE)
# 
# message("Marker analysis complete! Results saved to CSV.")

# ==============================================================================
# STOP HERE!
# Open 'Banksy_VZ_top10_Markers_Res0.8.csv'. 
# You must map the new Banksy cluster numbers to their biological identities
# before running the code below.
# ==============================================================================

# 5. RENAME IDENTITIES
# (Update these mappings based on your new CSV results!)
new_labels <- c(
  "1" = "GABAergic Progenitors", "2" = "Bergmann Glia", "3" = "Glial Progenitors", 
  "4" = "Maturing PCs", "5" = "Astrocytes", "6" = "Migratory PCs", 
  "7" = "Migratory PCs", "8" = "Glutamatergic DCN", "9" = "PC Progenitors", 
  "10" = "GABAergic DCN", "11" = "Cycling Cells", "12" = "Intermediate PCs", 
  "13" = "Rostral PCs", "14" = "Mature PCs", "15" = "Ependymal Cells", 
  "16" = "OPCs", "17" = "Ependymal Cells"
)
obj <- RenameIdents(obj, new_labels)

# 6. APPLY STANDARDIZED ORDERING
# Ensure vz_subcluster_order covers all your new Banksy cell types
obj$VZ_B_subcluster <- factor(Idents(obj), levels = vz_B_subcluster_order)
Idents(obj) <- "VZ_B_subcluster"

DefaultAssay(obj) <- "Xenium"

# 7. UMAP WITH CONSISTENT COLORS
# UPDATED: reduction changed to umap_banksy
p <- DimPlot(obj,
             reduction = "umap_banksy",
             label = TRUE,
             label.size = 4,
             repel = TRUE,
             cols = vz_B_palette) +
  ggtitle("Banksy Spatial VZ Subclusters") +
  theme(legend.text = element_text(size = 8))

ggsave(here("outputs", "XenAld_VZ_Banksy_Plots", "Banksy_VZ_Subcluster_UMAP.png"),
       p, width = 12, height = 9, dpi = 300)

# 8. DOTPLOT
obj$VZ_B_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(vz_B_subcluster_order))
Idents(obj) <- "VZ_B_subcluster_rev"

plot_features <- unique(as.vector(unlist(vz_B_markers)))
available_genes <- rownames(obj[["Xenium"]])
plot_features <- plot_features[plot_features %in% available_genes]

p1 <- DotPlot(obj,
              features = plot_features,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Banksy VZ Subcluster Spatial Markers")

ggsave(here("outputs", "XenAld_VZ_Banksy_Plots", "Banksy_VZ_SubclusterMarker_DotPlot.png"),
       p1, width = 12, height = 8, dpi = 300)

# 9. SAVE FINAL RESULT
Idents(obj) <- "VZ_B_subcluster"
output_path <- here("outputs", "XenAld_VZ_Banksy_RDS", "Xenium_VZ_Banksy_Annotated.rds")
saveRDS(obj, output_path, compress = FALSE)

plan("sequential")
check_mem("PIPELINE COMPLETE")