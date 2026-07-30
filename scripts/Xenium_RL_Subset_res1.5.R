library(here)
library(Seurat)
library(dplyr)

# Capture Slurm Array ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No Task ID provided. This script should be run via Slurm array.")
}
task_id <- as.numeric(args[1])

# Define your 15 samples explicitly to match the array index
samples <- c(
  #"GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G", "FB330_1_X_G", "GZFB5_X_G", 
  "GZFB_12_X_G_1", "GZFB_12_X_G_2", "GZFB_12_X_G_3",
  #"FB78_X_G",  "GZFB_1_X_G", 
  "GZFB_12_X_G_4", "GZFB_12_X_G_5", "GZFB_9_X_G_3"
  #"GZFB_9_X_G_1", "GZFB_9_X_G_2", 
)

current_sample <- samples[task_id]

# Use here() to build the input path
# Structure: project_root/outputs/XeniumAldingerABT_RDS/...
input_path <- here("outputs", "Xenium_AldingerABT_Res1.5_PCW_RDS", paste0(current_sample, "_Aldinger_annotated.rds"))

message("Loading: ", input_path)
if (!file.exists(input_path)) stop("File not found: ", input_path)

obj <- readRDS(input_path)

# Identify the clusters you want to pull out for re-analysis
# (Update these strings to match your exact 'annotated' cluster names)
target_clusters <- c("RL", "Granule", "UBC") 

# 2. Find which of those actually exist in the object's active identities
existing_clusters <- intersect(target_clusters, levels(Idents(obj)))

# 3. Subset using only the clusters that were found
if (length(existing_clusters) > 0) {
  obj_subset <- subset(obj, idents = existing_clusters)
  message("Subset successful using: ", paste(existing_clusters, collapse = ", "))
} else {
  stop("None of the target clusters were found in the object.")
}

# Define and create output directory using here()
output_dir <- here("outputs", "XenAld_RL_Subsets_Res1.5_RDS")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_path <- file.path(output_dir, paste0(current_sample, "_RLsubset.rds"))
saveRDS(obj_subset, file = output_path)

message("Successfully saved subset for ", current_sample, " to ", output_path)