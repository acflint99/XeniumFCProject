# Clear the environment
rm(list = ls())

library(here)
library(Seurat)

# Source the function
source(here("scripts", "XeniumCropCerebellum.R"))

# Load and crop a sample
xenium_cereb <- XeniumCropCerebellum("GZFB5_X_G")

# Source the function
source(here("scripts", "XeniumQC.R"))

# Filter out low quality cells
xenium_cereb_QC <- qc_xenium(xenium_cereb, sample_name = "GZFB5_X_G")

# Source the function
source(here("scripts", "XeniumNormCluster.R"))

xenium_cereb_proc <- process_xenium_clusters(xenium_cereb_QC, sample_name = "GZFB5_X_G")



