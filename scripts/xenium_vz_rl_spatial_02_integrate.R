#!/usr/bin/env Rscript

# Sketch-based Harmony integration of the validated 34-sample spatial merge.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_spatial_02_integrate.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_path <- file.path(output_root, "xenium", "vz_rl", "spatial", "01_merged", "rds", "XenAld_VZRL_spatial_merged.rds")
input_manifest_path <- file.path(output_root, "xenium", "vz_rl", "spatial", "01_merged", "tables", "XenAld_VZRL_spatial_merged_manifest.csv")
rds_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "02_integrated", "rds")
plot_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "02_integrated", "plots")
table_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "02_integrated", "tables")
output_path <- file.path(rds_dir, "XenAld_VZRL_spatial_integrated.rds")
output_manifest_path <- file.path(table_dir, "XenAld_VZRL_spatial_integrated_manifest.csv")
plot_paths <- file.path(
  plot_dir,
  c(
    "XenAld_VZRL_Spatial_Integrated_SampleID_UMAP.tif",
    "XenAld_VZRL_Spatial_Integrated_BroadClusters_UMAP.tif",
    "XenAld_VZRL_Spatial_Integrated_Subclusters_UMAP.tif"
  )
)
expected_outputs <- c(output_path, output_manifest_path, plot_paths)

if (dry_run) {
  cat("Spatial merged input:", input_path, "\n")
  cat("Spatial merged input exists:", file.exists(input_path), "\n")
  cat("Spatial merge manifest exists:", file.exists(input_manifest_path), "\n")
  write.table(data.frame(output = expected_outputs, exists = file.exists(expected_outputs)),
              row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("Spatial merged input not found: ", input_path)
if (!file.exists(input_manifest_path)) stop("Spatial merge manifest not found: ", input_manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop("Refusing to overwrite spatial integration outputs:\n- ",
       paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
}

input_manifest <- read.csv(input_manifest_path, stringsAsFactors = FALSE)
if (nrow(input_manifest) != length(sample_ids) ||
    !setequal(input_manifest$sample_id, sample_ids)) {
  stop("Spatial merge manifest must contain exactly the configured 34 samples.")
}

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
options(future.globals.maxSize = 100 * 1024^3)
plan(multisession, workers = workers)
on.exit(plan(sequential), add = TRUE)
set.seed(config$runtime$random_seed)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading validated spatial merge...")
obj <- readRDS(input_path)
required_metadata <- c("sample_id", "orig.ident", "consensus_label", "comb_subcluster", "PCW")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("Spatial merge lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (ncol(obj) != sum(input_manifest$cells)) stop("Spatial merge cell count does not match its manifest.")
if (!setequal(unique(as.character(obj$sample_id)), sample_ids)) {
  stop("Spatial merge does not contain exactly the configured sample IDs.")
}

obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])
DefaultAssay(obj) <- "Xenium"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30, seed.use = config$runtime$random_seed)

sketch_size <- min(100000L, ncol(obj))
message("Sketching ", sketch_size, " cells...")
obj <- SketchData(
  object = obj, ncells = sketch_size, method = "LeverageScore",
  sketched.assay = "sketch"
)
DefaultAssay(obj) <- "Xenium"
sketch_cells <- colnames(obj[["sketch"]])
sketch_counts <- LayerData(obj, assay = "Xenium", layer = "counts")[, sketch_cells, drop = FALSE]
obj[["sketch"]] <- CreateAssay5Object(counts = sketch_counts)
DefaultAssay(obj) <- "sketch"
sketch_sample_ids <- obj@meta.data[sketch_cells, "sample_id"]
if (anyNA(sketch_sample_ids)) stop("Sketched cells contain missing sample IDs.")
obj[["sketch"]] <- split(obj[["sketch"]], f = sketch_sample_ids)

n_features <- min(2000L, nrow(obj[["sketch"]]))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_features)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE, seed.use = config$runtime$random_seed)
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration, orig.reduction = "pca",
  new.reduction = "integrated.harmony", group.by = "sample_id", verbose = TRUE
)

obj[["sketch"]] <- JoinLayers(obj[["sketch"]])
obj <- ProjectIntegration(
  object = obj, sketched.assay = "sketch", assay = "Xenium",
  reduction = "integrated.harmony"
)
obj <- RunUMAP(
  obj, reduction = "integrated.harmony", dims = 1:30,
  reduction.name = "umap.harmony", return.model = TRUE,
  seed.use = config$runtime$random_seed
)
obj <- ProjectData(
  object = obj, sketched.assay = "sketch", assay = "Xenium",
  sketched.reduction = "integrated.harmony",
  full.reduction = "integrated.harmony", dims = 1:30
)
DefaultAssay(obj) <- "Xenium"

observed_consensus <- unique(as.character(obj$consensus_label))
unexpected_consensus <- setdiff(observed_consensus, celltype_order)
if (length(unexpected_consensus)) stop("Unexpected consensus labels: ", paste(unexpected_consensus, collapse = ", "))
consensus_levels <- intersect(celltype_order, observed_consensus)
obj$consensus_label <- factor(obj$consensus_label, levels = consensus_levels)
observed_subclusters <- unique(as.character(obj$comb_subcluster))
unexpected_subclusters <- setdiff(observed_subclusters, master_subcluster_order)
if (length(unexpected_subclusters)) stop("Unexpected combined subclusters: ", paste(unexpected_subclusters, collapse = ", "))
subcluster_levels <- intersect(master_subcluster_order, observed_subclusters)
obj$comb_subcluster <- factor(obj$comb_subcluster, levels = subcluster_levels)

p_sample <- DimPlot(obj, reduction = "umap.harmony", group.by = "sample_id", raster = TRUE) +
  ggtitle("Integrated by Sample")
Cairo::CairoTIFF(plot_paths[[1]], width = 14, height = 14, units = "in", res = 600)
print(p_sample)
grDevices::dev.off()

p_broad <- DimPlot(
  obj, reduction = "umap.harmony", group.by = "consensus_label",
  label = TRUE, raster = TRUE, cols = cluster_colors[consensus_levels]
) + ggtitle("Integrated: Broad Clusters")
Cairo::CairoTIFF(plot_paths[[2]], width = 14, height = 14, units = "in", res = 600)
print(p_broad)
grDevices::dev.off()

p_subcluster <- DimPlot(
  obj, reduction = "umap.harmony", group.by = "comb_subcluster",
  label = TRUE, label.size = 3.5, repel = TRUE, label.box = TRUE,
  raster = TRUE, alpha = 0.6, cols = subcluster_palette[subcluster_levels]
) + ggtitle("Integrated: Subclusters") +
  theme(legend.text = element_text(size = 9), legend.key.size = grid::unit(0.5, "cm")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))
Cairo::CairoTIFF(plot_paths[[3]], width = 14, height = 14, units = "in", res = 600)
print(p_subcluster)
grDevices::dev.off()

plan(sequential)
saveRDS(obj, output_path, compress = FALSE)
input_info <- file.info(input_path)
output_manifest <- data.frame(
  input_path = input_path, input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime), cells = ncol(obj),
  samples = length(unique(obj$sample_id)), sketch_cells = sketch_size,
  random_seed = config$runtime$random_seed, stringsAsFactors = FALSE
)
write.csv(output_manifest, output_manifest_path, row.names = FALSE)
message("Saved integrated spatial object: ", output_path)
