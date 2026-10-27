
library(Seurat)

# 1. Define input and output directories
input_dir <- "/home/acflint/R/Projects/XeniumFCProject/outputs/xenium/vz_rl/nonspatial/04_processed/rds"

output_dir <- file.path(input_dir, "celltype_subsets")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load processed Xenium object
xenium <- readRDS(
  file.path(input_dir, "Xenium_VZRL_processed.rds")
)

# 3. Check cell type labels
table(xenium$corrected_broad_label, useNA = "ifany")

# 4. Subset Purkinje cells
purkinje <- subset(
  xenium,
  subset = corrected_broad_label == "Purkinje"
)

# 5. Subset Granule cells
granule <- subset(
  xenium,
  subset = corrected_broad_label == "Granule"
)

# 6. Print cell counts
cat("Purkinje cells:", ncol(purkinje), "\n")
cat("Granule cells:", ncol(granule), "\n")

# 7. Save each subset separately
saveRDS(
  purkinje,
  file = file.path(output_dir, "Xenium_VZRL_Purkinje.rds"),
  compress = FALSE
)

saveRDS(
  granule,
  file = file.path(output_dir, "Xenium_VZRL_Granule.rds"),
  compress = FALSE
)

cat("Both subsets saved successfully!\n")
