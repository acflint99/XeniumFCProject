#!/usr/bin/env Rscript

# Integrate retained Other/unresolved populations as an independent residual
# branch using the established regional dimensionality and clustering settings.

rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid <- c("--dry-run", "--overwrite")
if (any(!args %in% valid)) stop("Usage: Rscript scripts/xenium_residual_03_integrate.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

root <- here(config$project$outputs_dir, "xenium", "residual")
input_path <- file.path(root, "02_merged", "rds", "Xenium_Residual_Merged.rds")
merge_manifest_path <- file.path(root, "02_merged", "tables", "Xenium_Residual_Merge_Manifest.csv")
output_path <- file.path(root, "03_integrated", "rds", "Xenium_Residual_Integrated.rds")
manifest_path <- file.path(root, "03_integrated", "tables", "Xenium_Residual_Integration_Manifest.csv")
plot_paths <- file.path(
  root, "03_integrated", "plots",
  paste0("Xenium_Residual_UMAP_Res_", c(0.5, 0.8, 1), ".tif")
)
outputs <- c(output_path, manifest_path, plot_paths)
if (dry_run) {
  compact_dry_run("Residual integration", inputs=c(input_path, merge_manifest_path), outputs=outputs)
  quit(save="no", status=if(all(file.exists(c(input_path, merge_manifest_path)))) 0L else 1L)
}
missing <- c(input_path, merge_manifest_path)[!file.exists(c(input_path, merge_manifest_path))]
if (length(missing)) stop("Missing residual integration input(s):\n- ", paste(missing, collapse="\n- "))
existing <- outputs[file.exists(outputs)]
if (length(existing) && !overwrite) stop("Refusing to overwrite residual integration output(s):\n- ", paste(existing, collapse="\n- "))

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(Cairo)
})
object <- readRDS(input_path)
required <- c(
  "orig.ident", "consensus_label", "regional_consensus_label",
  "residual_reference_route", "residual_reference_support_n", "PCW"
)
missing_metadata <- setdiff(required, names(object[[]]))
if (length(missing_metadata)) stop("Residual merge lacks: ", paste(missing_metadata, collapse=", "))
if (ncol(object) == 0L || anyDuplicated(Cells(object))) stop("Residual merge is empty or has duplicate cells.")
if (any(!object$residual_reference_route %in% c("route_other", "route_unresolved"))) {
  stop("Residual merge contains a VZ/RL or invalid route.")
}
dimensions <- as.integer(config$regional_subsets$integration_dimensions_default)
neighbors <- as.integer(config$regional_subsets$neighbor_k_post_qc)
resolutions <- c(0.5, 0.8, 1)
DefaultAssay(object) <- if ("Xenium" %in% Assays(object)) "Xenium" else DefaultAssay(object)
set.seed(config$runtime$random_seed)
object <- NormalizeData(object, verbose=FALSE)
object <- FindVariableFeatures(object, selection.method="vst", nfeatures=2000, verbose=FALSE)
if (length(VariableFeatures(object)) <= dimensions) stop("Residual branch has too few variable genes for PCA.")
object <- ScaleData(object, features=VariableFeatures(object), verbose=FALSE)
object <- RunPCA(object, features=VariableFeatures(object), npcs=dimensions, verbose=FALSE, seed.use=config$runtime$random_seed)
object <- RunHarmony(
  object, group.by.vars="orig.ident", dims.use=seq_len(dimensions),
  reduction.save="harmony", verbose=TRUE
)
object <- RunUMAP(
  object, reduction="harmony", dims=seq_len(dimensions),
  reduction.name="umap_harmony", n.neighbors=neighbors,
  seed.use=config$runtime$random_seed, verbose=FALSE
)
object <- FindNeighbors(
  object, reduction="harmony", dims=seq_len(dimensions),
  k.param=neighbors, verbose=FALSE
)
object <- FindClusters(object, resolution=resolutions, verbose=FALSE)

assay_prefix <- DefaultAssay(object)
plot_dir <- dirname(plot_paths[[1]])
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  column <- paste0(assay_prefix, "_snn_res.", resolution)
  if (!column %in% names(object[[]])) stop("Residual clustering lacks ", column, ".")
  plot <- DimPlot(
    object, reduction="umap_harmony", group.by=column, label=TRUE,
    repel=TRUE, raster=TRUE
  ) + ggtitle(paste("Residual branch - resolution", resolution)) + theme_classic()
  Cairo::CairoTIFF(plot_paths[[i]], width=10, height=8, units="in", res=600)
  print(plot); grDevices::dev.off()
}
input_info <- file.info(input_path)
manifest <- data.frame(
  created_at_utc=format(Sys.time(), tz="UTC", usetz=TRUE),
  input_rds=normalizePath(input_path, mustWork=TRUE),
  input_size=as.numeric(input_info$size), input_mtime=as.numeric(input_info$mtime),
  cells=ncol(object), samples=length(unique(as.character(object$orig.ident))),
  variable_features=2000L, dimensions=dimensions, neighbor_k=neighbors,
  resolutions=paste(resolutions, collapse="|"),
  integration_variable="orig.ident", random_seed=config$runtime$random_seed,
  stringsAsFactors=FALSE
)
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(manifest_path), recursive=TRUE, showWarnings=FALSE)
saveRDS(object, output_path, compress=FALSE)
write.csv(manifest, manifest_path, row.names=FALSE)
message("Saved residual integration: ", output_path)
