# Clear the environment
rm(list = ls())

# Load required libraries
library(Seurat)
library(ggplot2)
library(here)

source(here("scripts", "color_palette.R"))

# Define your directory and sample list
input_dir <- here("outputs", "Xenium_AldingerABT_Res1.5_RDS")

output_dir <- here("outputs", "Xenium_AldingerABT_Res1.5_GlobalPlots")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sample_list <- c(
  "GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G", 
  "FB330_1_X_G", "FB78_X_G", "GZFB5_X_G", "GZFB_12_X_G_1", 
  "GZFB_12_X_G_2", "GZFB_12_X_G_3", "GZFB_12_X_G_4", "GZFB_12_X_G_5", 
  "GZFB_1_X_G", "GZFB_9_X_G_1", "GZFB_9_X_G_2", "GZFB_9_X_G_3"
)

# Loop through each sample
for (sample_name in sample_list) {
  
  # FIX: Use file.path() instead of paste0() to ensure slashes are correct
  file_path <- file.path(input_dir, paste0(sample_name, "_Aldinger_annotated.rds"))
  
  # Check if the file actually exists
  if (file.exists(file_path)) {
    cat("Processing:", sample_name, "\n")
    
    # Load the Seurat object
    xenium_obj <- readRDS(file_path)
    
    # Generate the plot
    # Note: Ensure 'cluster_colors' is defined in your environment beforehand
    p_wei <- ImageDimPlot(
      xenium_obj, 
      group.by = "cluster_weighted", 
      size = 0.75, 
      cols = cluster_colors
    ) + 
      ggtitle(paste(sample_name, "(Weighted)"))
    
    # Save with a specific background color to ensure ggsave doesn't add white
    Cairo::CairoTIFF(
      filename = file.path(output_dir, paste0(sample_name, "_GlobalSpatial_plot0.75.tif")),
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
