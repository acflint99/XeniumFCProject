#!/usr/bin/env Rscript

# Process the complete combined VZ/RL merge, or create plots from the saved
# processed object without repeating integration.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(future)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")

# Seurat's uwot implementation derives its thread count from the active future
# plan. Restrict multisession use to UMAP itself so full Seurat objects are not
# exported to workers during normalization, scaling, PCA, or Harmony.
run_umap_with_workers <- function(object, ...) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (workers > 1L) future::plan(future::multisession, workers = workers)
  Seurat::RunUMAP(object, ...)
}

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--process-only", "--plots-only", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) stop("All arguments must be named options.")
process_only <- "--process-only" %in% args
plots_only <- "--plots-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (process_only && plots_only) stop("Choose --process-only or --plots-only, not both.")
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
rds_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS")
plot_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_Plots")
input_path <- file.path(rds_dir, "XenAld_VZRL_clean_merge.rds")
input_manifest_path <- file.path(rds_dir, "XenAld_VZRL_clean_merge_manifest.csv")
processed_path <- file.path(rds_dir, "XenAld_VZRL_clean_merge_processed.rds")
processed_manifest_path <- file.path(rds_dir, "XenAld_VZRL_clean_merge_processed_manifest.csv")

plot_names <- c(
  "XenAld_Batch_Comp_UMAP.tif",
  "Xenium_Merged_ConsensusLabel_UMAP.tif",
  "XenAld_Merged_consensus_label_DotPlot_markers.tif",
  "XenAld_Merged_consensus_label_DotPlot_markers.pdf",
  "XenAld_Purkinje_Specific_DotPlot.tif",
  "XenAld_Purkinje_Specific_DotPlot.pdf",
  "XenAld_Merged_Cluster_Top5_Heatmap.tif",
  "XenAld_Merged_Cluster_Top5_Heatmap.pdf",
  "XenAld_Signaling_Genes_dual_VlnPlot.tif",
  "XenAld_Signaling_Genes_dual_VlnPlot.pdf"
)
plot_paths <- file.path(plot_dir, plot_names)

if (dry_run) {
  cat("Combined merge:", input_path, "\n")
  cat("Combined merge exists:", file.exists(input_path), "\n")
  cat("Input manifest exists:", file.exists(input_manifest_path), "\n")
  cat("Processed object:", processed_path, "\n")
  cat("Processed object exists:", file.exists(processed_path), "\n")
  cat("Plot outputs existing:", sum(file.exists(plot_paths)), "of", length(plot_paths), "\n")
  cat("Requested mode:", if (process_only) "processing" else if (plots_only) "plotting" else "inspection only", "\n")
  quit(save = "no", status = 0L)
}
if (!process_only && !plots_only) {
  stop("Choose --process-only for integration or --plots-only after processing.")
}

if (process_only) {
  if (!file.exists(input_path)) stop("Combined merge not found: ", input_path)
  if (!file.exists(input_manifest_path)) stop("Combined merge manifest not found: ", input_manifest_path)
  existing_outputs <- c(processed_path, processed_manifest_path)[
    file.exists(c(processed_path, processed_manifest_path))
  ]
  if (length(existing_outputs) && !overwrite) {
    stop("Refusing to overwrite combined processed outputs:\n- ",
         paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
  }

  input_manifest <- read.csv(input_manifest_path, stringsAsFactors = FALSE)
  if (nrow(input_manifest) != length(sample_ids) ||
      !setequal(input_manifest$sample_id, sample_ids)) {
    stop("Combined input manifest must contain exactly the configured 34 samples.")
  }
  merged_obj <- readRDS(input_path)
  required_metadata <- c("orig.ident", "comb_subcluster", "consensus_label", "VZ_subcluster", "RL_subcluster", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(merged_obj[[]]))
  if (length(missing_metadata)) stop("Combined merge lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (ncol(merged_obj) != sum(input_manifest$cells)) {
    stop("Combined merge cell count does not match its manifest.")
  }
  if (!setequal(unique(as.character(merged_obj$orig.ident)), sample_ids)) {
    stop("Combined merge does not contain exactly the configured sample IDs.")
  }

  set.seed(config$runtime$random_seed)
  plan(sequential)
  DefaultAssay(merged_obj) <- "Xenium"
  merged_obj <- NormalizeData(merged_obj)
  merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000)
  merged_obj <- JoinLayers(merged_obj)
  merged_obj <- ScaleData(merged_obj)
  merged_obj <- RunPCA(merged_obj, npcs = 30, seed.use = config$runtime$random_seed)
  merged_obj <- run_umap_with_workers(
    merged_obj, reduction = "pca", dims = 1:30,
    reduction.name = "umap.unintegrated", seed.use = config$runtime$random_seed
  )
  merged_obj <- RunHarmony(
    merged_obj, group.by.vars = "orig.ident", reduction.use = "pca",
    reduction.save = "harmony"
  )
  merged_obj <- run_umap_with_workers(
    merged_obj, reduction = "harmony", dims = 1:30,
    reduction.name = "umap.harmony", seed.use = config$runtime$random_seed
  )
  saveRDS(merged_obj, processed_path, compress = FALSE)
  input_info <- file.info(input_path)
  processed_manifest <- data.frame(
    input_path = input_path, input_size = as.numeric(input_info$size),
    input_mtime = as.numeric(input_info$mtime), cells = ncol(merged_obj),
    samples = length(unique(merged_obj$orig.ident)),
    random_seed = config$runtime$random_seed, stringsAsFactors = FALSE
  )
  write.csv(processed_manifest, processed_manifest_path, row.names = FALSE)
  message("Saved processed combined object: ", processed_path)
  quit(save = "no", status = 0L)
}

if (!file.exists(processed_path)) stop("Processed combined object not found: ", processed_path)
if (!file.exists(processed_manifest_path)) stop("Processed manifest not found: ", processed_manifest_path)
existing_plots <- plot_paths[file.exists(plot_paths)]
if (length(existing_plots) && !overwrite) {
  stop("Refusing to overwrite combined plot outputs:\n- ",
       paste(existing_plots, collapse = "\n- "), "\nUse --overwrite only after review.")
}
processed_manifest <- read.csv(processed_manifest_path, stringsAsFactors = FALSE)
input_info <- file.info(input_path)
input_is_current <- file.exists(input_path) && nrow(processed_manifest) == 1L &&
  isTRUE(all.equal(as.numeric(processed_manifest$input_size), as.numeric(input_info$size))) &&
  isTRUE(all.equal(as.numeric(processed_manifest$input_mtime), as.numeric(input_info$mtime)))
if (!input_is_current) stop("Combined merge changed after processing. Rerun --process-only.")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
merged_obj <- readRDS(processed_path)
if (ncol(merged_obj) != processed_manifest$cells[[1]]) stop("Processed object cell count does not match its manifest.")
required_reductions <- c("umap.unintegrated", "umap.harmony")
if (!all(required_reductions %in% Reductions(merged_obj))) {
  stop("Processed object lacks reduction(s): ", paste(setdiff(required_reductions, Reductions(merged_obj)), collapse = ", "))
}
set.seed(config$runtime$random_seed)

p_pre <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "orig.ident", raster = TRUE) +
  ggtitle("Pre-Harmony (PCA UMAP)")
p_post <- DimPlot(merged_obj, reduction = "umap.harmony", group.by = "orig.ident", raster = TRUE) +
  ggtitle("Post-Harmony (Harmony UMAP)")
Cairo::CairoTIFF(plot_paths[[1]], width = 16, height = 7, units = "in", res = 600)
print(p_pre + p_post)
grDevices::dev.off()

  observed_consensus <- unique(as.character(merged_obj$consensus_label))
  unexpected_consensus <- setdiff(observed_consensus, celltype_order)
  if (length(unexpected_consensus)) {
    stop("Processed object has consensus labels absent from celltype_order: ",
         paste(unexpected_consensus, collapse = ", "))
  }
  consensus_levels <- intersect(celltype_order, observed_consensus)
validate_palette(setdiff(consensus_levels, "Unknown"))
p_consensus <- DimPlot(
  merged_obj, reduction = "umap.harmony", group.by = "consensus_label",
  label = TRUE, label.size = 4, label.box = TRUE, raster = TRUE,
  pt.size = 0.5, alpha = 0.8, cols = cluster_colors[consensus_levels]
) + ggtitle("Xenium Merged UMAP: Consensus Labels") + theme_classic()
Cairo::CairoTIFF(plot_paths[[2]], width = 12, height = 9, units = "in", res = 600)
print(p_consensus)
grDevices::dev.off()

existing_markers <- lapply(markers, function(x) intersect(x, rownames(merged_obj)))
existing_markers <- existing_markers[lengths(existing_markers) > 0L]
Idents(merged_obj) <- factor(merged_obj$consensus_label, levels = rev(consensus_levels))
p_dot <- DotPlot(merged_obj, features = existing_markers, assay = "Xenium") +
  RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("Xenium Merged Marker Expression (LogNorm)")
Cairo::CairoTIFF(plot_paths[[3]], width = 14, height = 6, units = "in", res = 600)
print(p_dot)
grDevices::dev.off()
ggsave(plot_paths[[4]], p_dot, device = grDevices::cairo_pdf, width = 14, height = 6)

pc_clusters <- c("Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs")
pc_markers <- c("FOXP1", "ITPR1", "COL5A1", "VSTM2L", "NDNF", "EBF2", "CALB1",
                "TRPC3", "NEFL", "ETV1", "PCDH10", "EBF1", "BCL11A", "RORB", "EN1")
pc_subset <- subset(merged_obj, subset = VZ_subcluster %in% pc_clusters)
pc_subset$VZ_subcluster <- factor(pc_subset$VZ_subcluster, levels = rev(pc_clusters))
Idents(pc_subset) <- "VZ_subcluster"
existing_pc_markers <- intersect(pc_markers, rownames(pc_subset))
if (!length(existing_pc_markers)) stop("No Purkinje markers are present in the Xenium assay.")
p_pc <- DotPlot(pc_subset, features = existing_pc_markers, assay = "Xenium",
                cols = c("lightgrey", "red"), dot.scale = 8) +
  RotatedAxis() + ggtitle("Purkinje Cell Lineage Marker Expression")
Cairo::CairoTIFF(plot_paths[[5]], width = 8, height = 5, units = "in", res = 600)
print(p_pc)
grDevices::dev.off()
ggsave(plot_paths[[6]], p_pc, device = grDevices::cairo_pdf, width = 8, height = 5)

Idents(merged_obj) <- factor(merged_obj$consensus_label, levels = consensus_levels)
heatmap_markers <- FindAllMarkers(
  merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
  max.cells.per.ident = 500, random.seed = config$runtime$random_seed
)
top_markers <- heatmap_markers %>% group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC, with_ties = FALSE) %>% pull(gene) %>% unique()
temp_obj <- ScaleData(merged_obj, features = top_markers, verbose = FALSE)
p_heatmap <- DoHeatmap(
  subset(temp_obj, downsample = 100), features = rev(top_markers),
  group.by = "consensus_label", group.colors = cluster_colors,
  size = 4, angle = 45, draw.lines = TRUE, raster = FALSE
) + scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
  ggtitle("Top 5 Markers per Cluster")
Cairo::CairoTIFF(plot_paths[[7]], width = 14, height = 12, units = "in", res = 600)
print(p_heatmap)
grDevices::dev.off()
ggsave(plot_paths[[8]], p_heatmap, device = grDevices::cairo_pdf, width = 14, height = 12)

kit_genes <- unique(c("APP", "TNFRSF21", "CADM3", "CADM4", "NECTIN3", "CNTN2", "L1CAM",
                      "CXCL12", "CXCR4", "EFNB2", "EPHA4", "GJA1", "NCAM1", "NRXN2",
                      "CLSTN1", "ADGRL1", "DAG1", "LRRTM1", "RELN", "VLDLR", "CD99", "NTF3", "NTRK2"))
kit_genes <- intersect(kit_genes, rownames(merged_obj))
kit_idents <- intersect(c("Granule", "Purkinje", "Glia", "GABA", "OPC", "UBC"), levels(Idents(merged_obj)))
if (!length(kit_genes) || !length(kit_idents)) stop("Signaling-gene violin inputs are absent.")
p_vln <- VlnPlot(
  merged_obj, features = kit_genes, pt.size = 0, ncol = 3,
  cols = cluster_colors[kit_idents], idents = kit_idents
)
Cairo::CairoTIFF(plot_paths[[9]], width = 20, height = 20, units = "in", res = 600)
print(p_vln)
grDevices::dev.off()
ggsave(plot_paths[[10]], p_vln, device = grDevices::cairo_pdf, width = 20, height = 20)

message("Saved combined plots without rerunning integration.")
