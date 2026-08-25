# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# Load required libraries
library(Seurat)
library(ggplot2)
library(here)

source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()

# Define your directory and sample list
input_dir <- here("outputs", "xenium", "vz_rl", "01_combined_labels", "rds")

output_dir <- here("outputs", "xenium", "vz_rl", "01_combined_labels", "plots", "consensus")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sample_list <- load_sample_manifest(config)$sample_id

# Loop through each sample
for (sample_name in sample_list) {
  
  # FIX: Use file.path() instead of paste0() to ensure slashes are correct
  file_path <- file.path(input_dir, paste0(sample_name, "_Ald_VZ_RL_QC_Subclusters.rds"))
  
  # Check if the file actually exists
  if (file.exists(file_path)) {
    cat("Processing:", sample_name, "\n")
    
    # Load the Seurat object
    xenium_obj <- readRDS(file_path)
    
    # Generate the plot
    # Note: Ensure 'cluster_colors' is defined in your environment beforehand
    p_wei <- ImageDimPlot(
      xenium_obj, 
      group.by = "consensus_label", 
      size = 1.5, 
      cols = cluster_colors
    ) + 
      ggtitle(paste(sample_name, "-Post QC (Consensus)"))
    
    # Save with a specific background color to ensure ggsave doesn't add white
    Cairo::CairoTIFF(
      filename = file.path(output_dir, paste0(sample_name, "_postQC_Consensus_GlobalSpatial_plot1.5.tif")),
      width = 12,
      height = 10,
      units = "in",
      res = 600,
      bg = "black"
    )
    print(p_wei)
    grDevices::dev.off()
    
    # Clean up memory (optional but recommended for large Xenium objects)
    rm(xenium_obj)
    gc()
    
  } else {
    warning(paste("File not found:", file_path))
  }
}
