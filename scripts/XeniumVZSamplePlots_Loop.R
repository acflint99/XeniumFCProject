# Clear the environment
rm(list = ls())

# 0. Load Libraries
source("renv/activate.R")
library(Seurat)
library(here)

# 1. Load the function and your palette/config
# This defines generate_spatial_reports() and your color variables
source(here("scripts", "Xenium_VZ_SamplePlots.R"))
source(here("scripts", "color_palette.R"))

# 3. Define Samples (run one line (3-4 samples) at a time)
samples_to_run <- c(
  "GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G"
  #"FB330_1_X_G", "FB78_X_G", "GZFB5_X_G", "GZFB_12_X_G_1"
  #"GZFB_12_X_G_2", "GZFB_12_X_G_3", "GZFB_12_X_G_4", "GZFB_12_X_G_5"
  #"GZFB_1_X_G", "GZFB_9_X_G_1", "GZFB_9_X_G_2", "GZFB_9_X_G_3"
)

# 3. Define Paths
input_base <- here("outputs", "XeniumAldingerABT_VZsubcluster_RDS")
output_base <- here("outputs", "XeniumAldingerABT_VZsubcluster_Results")

# 4. Run Sequential Loop
message("Starting sequential processing...")

for (s in samples_to_run) {
  
  input_path <- file.path(input_base, paste0(s, "_Ald_VZ.rds"))
  
  if (!file.exists(input_path)) {
    warning(paste("File not found, skipping:", input_path))
    next
  }
  
  message(paste(">>> Processing Sample:", s))
  
  tryCatch({
    generate_spatial_reports(
      sample_id = s, 
      input_rds_path = input_path, 
      base_plot_dir = output_base
    )
  }, error = function(e) {
    message(paste("Error in sample", s, ":", e$message))
    
    # Create error log
    err_path <- file.path(output_base, s)
    if(!dir.exists(err_path)) dir.create(err_path, recursive = TRUE)
    writeLines(as.character(e), file.path(err_path, "error_log.txt"))
  })
}

message("All tasks complete.")