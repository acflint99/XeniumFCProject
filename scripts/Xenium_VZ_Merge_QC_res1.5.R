# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(harmony)
library(dplyr)
library(future)
library(future.apply)
library(ggplot2)
library(patchwork)

# Load your new palette and order
source(here("scripts", "color_palette.R"))


# 1. PARALLELIZATION 
plan("multisession", workers = 8) 
options(future.globals.maxSize = 200 * 1024^3)

plot_path <- here("outputs", "XenAld_VZ_QC_Res1.5_Plots")
if(!dir.exists(plot_path)) dir.create(plot_path)

table_path <- here("outputs", "XenAld_VZ_QC_Res1.5_Tables")
if(!dir.exists(table_path)) dir.create(table_path, recursive = TRUE)

merged_path <- here("outputs", "XenAld_VZ_Res1.5_RDS", "Xenium_VZ_Res1.5_32726.rds")
obj <- readRDS(merged_path)

# This merges the 15 separate sample layers into one unified matrix
obj <- JoinLayers(obj)

# 2. SET THE IDENTITY
# Ensure we are looking at the resolution you liked best (e.g., 0.3)
Idents(obj) <- "Xenium_snn_res.0.8"

# Export Merged Cluster QC Statistics
raw_merged_qc_stats <- obj@meta.data %>%
  group_by(Xenium_snn_res.0.8) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(raw_merged_qc_stats, 
          file.path(table_path, "XenAld_VZ_RawSubcluster_QC_Summary.csv"), 
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "XenAld_VZ_RawSubcluster_QC_Violins.pdf")
pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)
  
  p5 <- VlnPlot(
    obj, 
    features = feat, 
    group.by = "Xenium_snn_res.0.8",
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
  print(p5)
}

# 5. Close the device to finalize the file
dev.off()

message("Individual QC plots saved to: ", pdf_path)



# 3. RUN FINDALLMARKERS (The Optimized Way)
message(Sys.time(), ": Starting Marker Identification...")

all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,          # Only look for upregulated genes (standard for cell types)
  min.pct = 0.25,           # Gene must be in 25% of the cluster
  logfc.threshold = 0.25,   # Minimum 1.28x fold change
  test.use = "wilcox",      # Standard fast test

  # --- THE SPEED TRICK ---
  # Downsampling to 1000 cells per cluster gives 99% the same results
  # but runs 5x faster on massive Xenium datasets.
  max.cells.per.ident = 1000
)


# 4. FILTER & SAVE TOP 10
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(top10_markers,
          file.path(table_path, "VZ_RawSubcluster_top10_Markers_Res0.8.csv"),
          row.names = FALSE)

message("Marker analysis complete! Results saved to CSV.")


# 2. Define Paths and "Master Key"
input_dir  <- here("outputs", "Xenium_ConsensusABT_Res1.5_RDS")
output_dir <- here("outputs", "Xenium_AldingerABT_VZ_QC_Res1.5_RDS")
if(!dir.exists(output_dir)) dir.create(output_dir)

# Update pattern to match your original whole objects
sample_files <- list.files(input_dir, pattern = "_Consensus_annotated\\.rds$", full.names = TRUE)

# Identify barcodes in Cluster 7 from the integrated 'obj'
# We use names() because these barcodes include the 'SampleName_' prefix
cells_to_remove <- colnames(obj)[obj$Xenium_snn_res.0.8 == "7"]

# 3. Run Parallel Processing (Modified to remove Cluster 7)
message("Starting filtering of Cluster 7 from whole objects...")

updated_status <- future_lapply(sample_files, function(f) {
  
  # Strip suffix to get sample name
  s_name <- gsub("_Consensus_annotated\\.rds", "", basename(f))
  
  # Load the WHOLE object
  temp_obj <- readRDS(f)
  
  # Convert local barcodes to integrated-style barcodes for matching
  # (e.g., "ATGC..." -> "SampleName_ATGC...")
  integrated_style_barcodes <- paste0(s_name, "_", colnames(temp_obj))
  
  # Identify which cells in this specific file are part of the 'cells_to_remove' list
  cells_is_cluster_7 <- integrated_style_barcodes %in% cells_to_remove
  
  # SUBSET: Keep only cells that are NOT in cluster 7
  # The '!' operator negates the logical vector
  temp_obj <- subset(temp_obj, cells = colnames(temp_obj)[!cells_is_cluster_7])
  
  # 4. Success Check
  removed_count <- sum(cells_is_cluster_7)
  remaining_count <- ncol(temp_obj)
  
  # Save the cleaned whole object
  out_path <- file.path(output_dir, paste0(s_name, "_Ald_VZ_QC.rds"))
  saveRDS(temp_obj, file = out_path, compress = FALSE)
  
  # Clean up memory
  rm(temp_obj, integrated_style_barcodes, cells_is_cluster_7)
  gc()
  
  return(paste0(s_name, ": Removed ", removed_count, " cells. Remaining: ", remaining_count))
}, future.seed = TRUE)

plan(sequential)
message("Done! Cluster 7 cells removed from all tissue objects.")
