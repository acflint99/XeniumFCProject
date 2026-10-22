#!/usr/bin/env Rscript
# Same stage sequence as VZ/RL: reprocess approved retained cells, or explicitly
# reuse the original integration when QC removed none. No QC/identity changes.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_05_reintegrate_post_qc.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite", "--no-removals")) ||
    all(c("--dry-run", "--overwrite") %in% args)) {
  stop("Usage: [--no-removals] [--dry-run|--overwrite]")
}
no_removals <- "--no-removals" %in% args
p <- residual_post_qc_paths(config)
q <- residual_qc_removal_paths(config)
# Permit path-only inspection before upstream jobs have produced QC outputs.
input <- if (no_removals) q$source else q$filtered
merge_manifest <- here(config$project$outputs_dir, "xenium", "residual", "02_merged",
                       "tables", "Xenium_Residual_Merge_Manifest.csv")
inputs <- c(input, merge_manifest, if (!no_removals) unlist(q[c("decision", "removed", "manifest", "provenance")]))
resolutions <- c(0.5, 0.8, 1) # Preserve established residual clustering settings.
plots <- file.path(dirname(dirname(p$input)), "plots",
                   paste0("Xenium_Residual_PostQC_UMAP_Res_", resolutions, ".tif"))
outputs <- c(p$input, p$provenance, plots)
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("Residual post-QC processing (paths only; no RDS load)", inputs, outputs,
                       checks = c(explicit_qc_choice = !no_removals || !file.exists(q$decision)))
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (any(!file.exists(inputs))) stop("Missing post-QC input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))
if (any(file.exists(outputs)) && !"--overwrite" %in% args) stop("Protected post-QC outputs exist; review before --overwrite.")
state <- load_residual_qc_state(config)
if (identical(state$active, no_removals)) stop("QC choice conflicts with approved decision; use --no-removals only after explicit no-removal review.")
if (!identical(normalizePath(input), normalizePath(state$input))) stop("Unexpected residual QC input.")
suppressPackageStartupMessages({library(Seurat); library(harmony); library(Cairo); library(ggplot2); library(future)})
future::plan(future::sequential)
input_info <- file.info(input)
object <- readRDS(input)
validate_residual_qc_object(object, state)
ids <- colnames(object)
metadata <- object[[]]
if (!length(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)) ||
    !identical(rownames(metadata), ids) || !"orig.ident" %in% names(metadata)) stop("Invalid residual cell IDs/metadata.")
sample_ids <- as.character(load_sample_manifest(config)$sample_id)
m <- read.csv(merge_manifest, stringsAsFactors = FALSE)
if (!all(c("sample_id", "residual_cells") %in% names(m)) ||
    anyDuplicated(m$sample_id) || !setequal(m$sample_id, sample_ids)) stop("Merge manifest must cover all active samples exactly once.")
m <- m[match(sample_ids, m$sample_id), , drop = FALSE]
expected <- if (state$active) state$manifest$remaining_cells else m$residual_cells
observed <- as.integer(table(factor(as.character(metadata$orig.ident), levels = sample_ids)))
if (anyNA(metadata$orig.ident) || any(!metadata$orig.ident %in% sample_ids) ||
    any(!is.finite(expected)) || any(expected < 0 | expected != floor(expected)) ||
    !identical(observed, as.integer(expected))) stop("Post-QC per-sample counts differ from the reviewed manifests.")
columns <- paste0("Xenium_snn_res.", resolutions)
if (!all(columns %in% names(metadata))) stop("Missing original residual cluster partitions.")
# Preserve every original metadata column, including pre-QC cluster assignments.
# Compare partitions with these preserved columns after reclustering.
for (column in columns) {
  name <- paste0("residual_preQC_", column)
  if (name %in% names(metadata)) stop("Pre-QC provenance column already exists: ", name)
  object[[name]] <- metadata[[column]]
}
dimensions <- as.integer(config$regional_subsets$integration_dimensions_default)
neighbors <- as.integer(config$regional_subsets$neighbor_k_post_qc)
seed <- as.integer(config$runtime$random_seed)
set.seed(seed)
DefaultAssay(object) <- "Xenium"
if (!no_removals) {
  object[["Xenium"]] <- JoinLayers(object[["Xenium"]])
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  if (length(VariableFeatures(object)) <= dimensions) stop("Too few variable genes for PCA.")
  object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
  object <- RunPCA(object, features = VariableFeatures(object), npcs = dimensions, seed.use = seed, verbose = FALSE)
  object <- RunHarmony(object, group.by.vars = "orig.ident", dims.use = seq_len(dimensions), reduction.save = "harmony", verbose = TRUE)
  object <- RunUMAP(object, reduction = "harmony", dims = seq_len(dimensions), reduction.name = "umap_harmony",
                    n.neighbors = neighbors, seed.use = seed, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "harmony", dims = seq_len(dimensions), k.param = neighbors, verbose = FALSE)
  object <- FindClusters(object, resolution = resolutions, verbose = FALSE)
}
if (!all(c("pca", "harmony", "umap_harmony") %in% Reductions(object)) ||
    !"Xenium_snn" %in% Graphs(object) || !all(columns %in% names(object[[]]))) stop("Incomplete post-QC representations or partitions.")
# Original partition columns are retained unchanged; post-QC partitions have
# separate names. seurat_clusters remains the original source identity.
for (column in columns) {
  object[[paste0("residual_postQC_", column)]] <- object[[column, drop = TRUE]]
  object[[column]] <- metadata[[column]]
}
if ("seurat_clusters" %in% names(metadata)) object$seurat_clusters <- metadata$seurat_clusters
Idents(object) <- paste0("residual_postQC_Xenium_snn_res.0.8")
if (!identical(colnames(object), ids) ||
    !identical(object[[]][, names(metadata), drop = FALSE], metadata)) stop("Post-QC processing changed existing metadata or cell IDs/order.")
validate_residual_qc_object(object, state)
object@misc$residual_post_qc <- list(reintegrated = !no_removals, seed = seed,
  decision_md5 = state$decision_md5, input = normalizePath(input))
for (path in outputs) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(resolutions)) {
  plot <- DimPlot(object, reduction = "umap_harmony", group.by = paste0("residual_postQC_", columns[i]), label = TRUE,
                  repel = TRUE, raster = TRUE) + ggtitle(paste("Residual post-QC resolution", resolutions[i]))
  Cairo::CairoTIFF(plots[i], width = 10, height = 8, units = "in", res = 600)
  print(plot); grDevices::dev.off()
}
if (!residual_qc_file_matches(input, input_info$size, as.numeric(input_info$mtime))) stop("QC input changed during processing.")
saveRDS(object, p$input, compress = FALSE)
info <- file.info(p$input)
write.csv(data.frame(input_rds = normalizePath(input), input_size = input_info$size,
  input_mtime = as.numeric(input_info$mtime), output_rds = normalizePath(p$input), output_size = info$size,
  output_mtime = as.numeric(info$mtime), cells_after_qc = length(ids), resolutions = paste(resolutions, collapse = "|"),
  random_seed = seed, decision_md5 = state$decision_md5, reintegrated = !no_removals,
  cluster_prefix = "residual_postQC_Xenium_snn_res.", dimensions = dimensions, neighbor_k = neighbors,
  merge_manifest_md5 = unname(tools::md5sum(merge_manifest)), stringsAsFactors = FALSE), p$provenance, row.names = FALSE)
message("Saved residual post-QC input. Old cluster labels require fresh marker/label review: ", p$input)
