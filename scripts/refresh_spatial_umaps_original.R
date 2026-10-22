#!/usr/bin/env Rscript
# Redraw only the saved spatial integration's four original UMAP TIFFs.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/refresh_spatial_umaps_original.R")
source(here("scripts","R","config.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly=TRUE)
if(any(!args %in% c("--dry-run","--overwrite"))||all(c("--dry-run","--overwrite") %in% args)) stop("Use --dry-run or --overwrite.")
root <- here(config$project$outputs_dir,"xenium","vz_rl","spatial","04_integrated")
input <- file.path(root,"rds","Xenium_VZRL_spatial_integrated.rds")
plot_paths <- file.path(root,"plots",c("Xenium_VZRL_Spatial_Integrated_SampleID_UMAP.tif",
  "Xenium_VZRL_Spatial_Integrated_BroadClusters_UMAP.tif","Xenium_VZRL_Spatial_Integrated_Subclusters_UMAP.tif",
  "Xenium_VZRL_Spatial_Integrated_CorrectedBroad_UMAP.tif"))
if("--dry-run" %in% args) {ok<-compact_dry_run("Saved spatial UMAP redraw",input,plot_paths);quit(save="no",status=if(ok)0L else 1L)}
if(!file.exists(input)) stop("Saved spatially integrated RDS missing: ",input)
if(any(file.exists(plot_paths))&&!"--overwrite" %in% args) stop("Existing plots require --overwrite.")
suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
if(!requireNamespace("Cairo",quietly=TRUE)) stop("Missing project Cairo dependency.")
source(here("scripts","color_palette.R"))
source(here("scripts","R","corrected_broad_labels.R"))
input_info <- file.info(input)
obj <- readRDS(input)
metadata <- obj[[]];cells<-colnames(obj)
required <- c("sample_id","consensus_label","comb_subcluster","corrected_broad_label")
if(!length(cells)||anyDuplicated(cells)||!identical(rownames(metadata),cells)||!all(required %in% names(metadata))||!"umap.harmony" %in% Reductions(obj)) stop("Invalid saved integrated object.")
if(!setequal(as.character(obj$sample_id),load_sample_manifest(config)$sample_id)) stop("Saved spatial object does not match the current sample manifest.")
if(!identical(as.character(obj$corrected_broad_label),map_corrected_broad_labels(obj$comb_subcluster))) stop("Saved corrected broad labels differ from the approved map; refresh annotations before plotting.")
set.seed(config$runtime$random_seed)
dir.create(dirname(plot_paths[[1]]),recursive=TRUE,showWarnings=FALSE)

# Plotting statements below match the native spatial integration exporter.
observed_consensus <- unique(as.character(obj$consensus_label))
consensus_levels <- order_consensus_levels(observed_consensus, celltype_order)
obj$consensus_label <- factor(obj$consensus_label, levels = consensus_levels)
observed_subclusters <- unique(as.character(obj$comb_subcluster))
unexpected_subclusters <- setdiff(
  observed_subclusters[!is_unknown_consensus_label(observed_subclusters)],
  master_subcluster_order
)
if (length(unexpected_subclusters)) stop("Unexpected combined subclusters: ", paste(unexpected_subclusters, collapse = ", "))
subcluster_levels <- order_consensus_levels(observed_subclusters, master_subcluster_order)
obj$comb_subcluster <- factor(obj$comb_subcluster, levels = subcluster_levels)
subcluster_colors <- subcluster_palette[subcluster_levels]
unknown_subclusters <- subcluster_levels[is_unknown_consensus_label(subcluster_levels)]
subcluster_colors[unknown_subclusters] <- resolve_cluster_colors(unknown_subclusters)
if (anyNA(subcluster_colors)) {
  stop(
    "Combined subcluster labels lack colors: ",
    paste(names(subcluster_colors)[is.na(subcluster_colors)], collapse = ", ")
  )
}

p_sample <- DimPlot(obj, reduction = "umap.harmony", group.by = "sample_id", raster = TRUE) +
  ggtitle("Integrated by Sample")
Cairo::CairoTIFF(plot_paths[[1]], width = 14, height = 14, units = "in", res = 600)
print(p_sample)
grDevices::dev.off()

p_broad <- DimPlot(
  obj, reduction = "umap.harmony", group.by = "consensus_label",
  label = TRUE, raster = TRUE, cols = resolve_cluster_colors(consensus_levels)
) + ggtitle("Integrated: Broad Clusters")
Cairo::CairoTIFF(plot_paths[[2]], width = 14, height = 14, units = "in", res = 600)
print(p_broad)
grDevices::dev.off()

p_subcluster <- DimPlot(
  obj, reduction = "umap.harmony", group.by = "comb_subcluster",
  label = TRUE, label.size = 3.5, repel = TRUE, label.box = TRUE,
  raster = TRUE, alpha = 0.6, cols = subcluster_colors
) + ggtitle("Integrated: Subclusters") +
  theme(legend.text = element_text(size = 9), legend.key.size = grid::unit(0.5, "cm")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))
Cairo::CairoTIFF(plot_paths[[3]], width = 14, height = 14, units = "in", res = 600)
print(p_subcluster)
grDevices::dev.off()

corrected_levels <- intersect(celltype_order, unique(as.character(obj$corrected_broad_label)))
obj$corrected_broad_label <- factor(
  as.character(obj$corrected_broad_label), levels = corrected_levels
)
corrected_colors <- resolve_cluster_colors(corrected_levels)
p_corrected <- DimPlot(
  obj, reduction = "umap.harmony", group.by = "corrected_broad_label",
  label = TRUE, raster = TRUE, cols = corrected_colors
) + ggtitle("Integrated: Corrected Broad Labels")
Cairo::CairoTIFF(plot_paths[[4]], width = 14, height = 14, units = "in", res = 600)
print(p_corrected)
grDevices::dev.off()

if(!identical(colnames(obj),cells)) stop("Cell IDs changed during plotting.")
for(column in required) if(!identical(as.character(obj[[column]][,1]),as.character(metadata[[column]]))) stop("Source labels changed during plotting.")
if(!identical(input_info[c("size","mtime")],file.info(input)[c("size","mtime")])) stop("Saved RDS changed during plotting.")
if(anyNA(file.info(plot_paths)$size)||any(file.info(plot_paths)$size<=0)) stop("Partial UMAP output set.")
message("PASS: four original spatial UMAPs redrawn; saved reductions and RDS retained.")
