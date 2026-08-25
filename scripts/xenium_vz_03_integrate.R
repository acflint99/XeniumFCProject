# Clear the environment
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- samples$sample_id

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_03_integrate.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

plot_dir <- here(config$project$outputs_dir, "xenium", "vz", "03_integrated", "plots")
rds_dir <- here(config$project$outputs_dir, "xenium", "vz", "03_integrated", "rds")
merged_path <- here(
  config$project$outputs_dir, "xenium", "vz", "02_merged", "rds",
  "Xenium_Merged_VZSubsets.rds"
)
merge_manifest_path <- here(
  config$project$outputs_dir, "xenium", "vz", "02_merged", "tables",
  "Xenium_Merged_VZSubsets_manifest.csv"
)
res_list <- c(0.5, 0.8, 1)
output_path <- file.path(rds_dir, "Xenium_VZ_Res1.5.rds")
plot_paths <- c(
  file.path(plot_dir, paste0("XenAld_VZ_UMAP_Res_", res_list, ".tif")),
  file.path(plot_dir, "XenAld_VZ_Batch_Comp_UMAP.tif"),
  file.path(plot_dir, "Xenium_VZ_ConsensusLabel_UMAP.tif")
)
expected_outputs <- c(output_path, plot_paths)

manifest_ready <- FALSE
manifest_rows <- NA_integer_
manifest_ids_match <- FALSE
manifest_columns_match <- FALSE
manifest_values_valid <- FALSE
if (file.exists(merge_manifest_path)) {
  merge_manifest_preview <- read.csv(merge_manifest_path, stringsAsFactors = FALSE)
  manifest_rows <- nrow(merge_manifest_preview)
  manifest_columns_match <- all(c("sample_id", "input_path", "cells", "PCW") %in% names(merge_manifest_preview))
  manifest_ids_match <- "sample_id" %in% names(merge_manifest_preview) &&
    setequal(as.character(merge_manifest_preview$sample_id), sample_ids)
  if (manifest_columns_match) {
    manifest_cells <- suppressWarnings(as.numeric(merge_manifest_preview$cells))
    manifest_pcw <- as.character(merge_manifest_preview$PCW)
    manifest_values_valid <- all(is.finite(manifest_cells) & manifest_cells > 0) &&
      all(!is.na(manifest_pcw) & nzchar(manifest_pcw))
  }
  manifest_ready <- manifest_rows == length(sample_ids) && manifest_columns_match &&
    manifest_ids_match && manifest_values_valid
}

if (dry_run) {
  cat("Merged VZ input:", merged_path, "\n")
  cat("Merged VZ input exists:", file.exists(merged_path), "\n")
  cat("Merge manifest:", merge_manifest_path, "\n")
  cat("Manifest rows:", manifest_rows, "\n")
  cat("Manifest columns valid:", manifest_columns_match, "\n")
  cat("Manifest cell counts and PCW valid:", manifest_values_valid, "\n")
  cat("Manifest sample IDs match config:", manifest_ids_match, "\n")
  write.table(
    data.frame(output = expected_outputs, exists = file.exists(expected_outputs)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(merged_path)) stop("Merged VZ input not found: ", merged_path)
if (!manifest_ready) stop("VZ merge manifest is missing, incomplete, or inconsistent: ", merge_manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ processing outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(patchwork)
})

# --- HELPER: MEMORY MONITOR ---
# This prints the current RAM usage to your console at every step
check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 2. PARALLELIZATION SETUP
plan("multisession", workers = 8) 
options(future.globals.maxSize = 150 * 1024^3) # 150GB limit for globals

check_mem("PIPELINE START")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# 3. LOAD DATA
obj <- readRDS(merged_path)
required_metadata <- c("orig.ident", "consensus_label", "PCW")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("Merged VZ input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!setequal(unique(as.character(obj$orig.ident)), sample_ids)) {
  stop("Merged VZ input does not contain exactly the 34 configured sample IDs.")
}
manifest_order <- match(sample_ids, as.character(merge_manifest_preview$sample_id))
expected_cells <- as.integer(merge_manifest_preview$cells[manifest_order])
observed_cells <- as.integer(table(factor(as.character(obj$orig.ident), levels = sample_ids)))
if (!identical(observed_cells, expected_cells)) {
  stop("Merged VZ per-sample cell counts do not match the merge manifest.")
}
expected_pcw <- as.character(merge_manifest_preview$PCW[manifest_order])
observed_pcw <- vapply(sample_ids, function(sample_id) {
  values <- unique(as.character(obj$PCW[as.character(obj$orig.ident) == sample_id]))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 1L) values[[1]] else NA_character_
}, character(1))
if (anyNA(observed_pcw) || !identical(unname(observed_pcw), unname(expected_pcw))) {
  stop("Merged VZ per-sample PCW values do not match the merge manifest.")
}
if (ncol(obj) != sum(expected_cells)) {
  stop("Merged VZ cell count does not match its merge manifest.")
}
Idents(obj) <- "consensus_label"
check_mem("DATA LOADED")

# 4. STANDARD WORKFLOW
message("Normalizing and Finding Variable Features...")
obj <- NormalizeData(obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message("Scaling data...")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
check_mem("POST-SCALE")

message("Running PCA...")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# 5. BATCH CHECK (UNCORRECTED UMAP)
message("Generating uncorrected UMAP...")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, 
               reduction.name = "umap_uncorrected", verbose = TRUE)
check_mem("POST-UNCORRECTED UMAP")

# 6. RUN HARMONY
message("Running Harmony integration...")
obj <- RunHarmony(obj, group.by.vars = "orig.ident", 
                  dims.use = 1:30, reduction.save = "harmony", verbose = TRUE)
check_mem("POST-HARMONY")

# 7. POST-HARMONY REDUCTION & CLUSTERING
message("Running Integrated UMAP and Neighbors...")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, 
               reduction.name = "umap_harmony", verbose = TRUE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = TRUE)

# 8. GENERATE CLUSTERS (RUN IN BULK FOR SPEED)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters in parallel...")

# Running all resolutions at once is faster with future/parallelization
obj <- FindClusters(obj, resolution = res_list, verbose = FALSE)
check_mem("POST-BATCH-CLUSTERING")

# Now loop only for sorting and plotting
for(res in res_list) {
  # Dynamically build the column name so it ALWAYS matches
  res_col <- paste0(assay_prefix, "_snn_res.", res)
  
  # 1 & 2. Sort numerically
  clusters <- unique(na.omit(obj@meta.data[[res_col]]))
  numeric_order <- sort(as.numeric(as.character(clusters)))
  
  # 3. Re-assign (using res_col, not the hardcoded Xenium string)
  obj@meta.data[[res_col]] <- factor(as.character(obj@meta.data[[res_col]]), 
                                     levels = as.character(numeric_order))
  
  message(Sys.time(), ": Generating & Saving UMAP for Resolution: ", res)
  
  p <- DimPlot(obj, 
               reduction = "umap_harmony", 
               group.by = res_col, 
               label = TRUE, 
               label.size = 5,
               label.box = TRUE,
               raster = TRUE, 
               pt.size = 0.6,
               alpha = 0.8) + 
    ggtitle(paste0("Refined VZ UMAP - Resolution ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  Cairo::CairoTIFF(
    filename = file.path(plot_dir, paste0("XenAld_VZ_UMAP_Res_", res, ".tif")),
    width = 10,
    height = 8,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
  
  rm(p)
  gc()
}

# 8. VISUALIZATION (Using Raster to prevent RStudio lag)
message("Saving comparison plots...")
# raster = TRUE converts points to pixels; essential for >100k cells
p1 <- DimPlot(obj, reduction = "umap_uncorrected", group.by = "orig.ident", raster = TRUE) + 
  NoLegend() + ggtitle("Pre-Harmony")

p2 <- DimPlot(obj, reduction = "umap_harmony", group.by = "orig.ident", raster = TRUE) + 
  ggtitle("Post-Harmony")

Cairo::CairoTIFF(
  filename = file.path(plot_dir, "XenAld_VZ_Batch_Comp_UMAP.tif"),
  width = 16,
  height = 7,
  units = "in",
  res = 600
)
print(p1 + p2)
grDevices::dev.off()

# 2. Generate the Plot
p_orig <- DimPlot(obj, 
                  reduction = "umap_harmony", 
                  group.by = "consensus_label",
                  label = TRUE, 
                  label.size = 4,
                  label.box = TRUE,       # Makes original labels easier to see
                  raster = TRUE, 
                  pt.size = 0.5, 
                  alpha = 0.8) +
  ggtitle("Xenium VZ Integrated UMAP: Consensus Labels") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# 3. Save the Plot
Cairo::CairoTIFF(
  filename = file.path(plot_dir, "Xenium_VZ_ConsensusLabel_UMAP.tif"),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p_orig)
grDevices::dev.off()

saveRDS(obj, output_path, compress = FALSE)

plan("sequential")

check_mem("PIPELINE COMPLETE - FILE SAVED")
