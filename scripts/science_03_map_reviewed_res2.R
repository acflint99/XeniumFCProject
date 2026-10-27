#!/usr/bin/env Rscript
# Apply the approved 2026-10-22 resolution-2 map, remove reviewed labels, reprocess.
# Expected HPC runtime: 1-8 h excluding queue; very low-confidence static estimate.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/science_03_map_reviewed_res2.R")
source(here::here("scripts", "R", "config.R"))
source(here::here("scripts", "R", "science_reviewed_mapping.R"))
options_science <- commandArgs(trailingOnly = TRUE)
allowed <- c("--map-reviewed", "--resolution=2", "--dry-run", "--overwrite")
if (any(!options_science %in% allowed) || all(c("--dry-run", "--overwrite") %in% options_science)) {
  stop("Use science_03_map_reviewed_res2.R [--dry-run|--overwrite]; mapping is only approved for resolution 2.")
}
config <- load_pipeline_config()
source_root <- here::here(config$project$outputs_dir, "references", "science", "reprocessed", "res2")
input_path <- file.path(source_root, "Science_subset_ccscore_harmony_clusters_res2.rds")
marker_input <- file.path(source_root, "Science_clusters_res2_top20_markers.csv")
review_path <- here::here("scripts", "reviews", "Science_res2_reviewed_labels_20261022.csv")
root <- file.path(source_root, "reviewed_retained")
prefix <- "Science_res2_reviewed"
paths <- setNames(file.path(root, paste0(prefix, c(
  "_retained_reprocessed.rds", "_UMAP_before_removal.pdf", "_UMAP_after_removal_original_coordinates.pdf",
  "_UMAP_after_reprocessing.pdf", "_raw_clusters_after_reprocessing_UMAP.pdf",
  "_broad_marker_DotPlot.pdf", "_cell_mapping.csv", "_cluster_retention.csv",
  "_sample_label_counts.csv", "_postreview_top20_markers.csv", "_postreview_marker_counts.csv",
  "_provenance.csv", "_sessionInfo.txt"
))), c("rds", "before", "after_removal", "after_reprocessing", "raw_umap", "dotplot",
       "cells", "retention", "sample_counts", "markers", "marker_counts", "provenance", "session"))
if (!file.exists(review_path)) stop("Missing approved Science review: ", review_path)
review <- read.csv(review_path, colClasses = c(cluster = "character"), stringsAsFactors = FALSE)
validate_science_review(review)
packages <- c("Seurat", "SeuratObject", "Matrix", "harmony", "dplyr", "ggplot2")
if ("--dry-run" %in% options_science) {
  ok <- compact_dry_run("Reviewed Science map/removal/reprocessing (paths only; no RDS loading)",
    c(input_path, marker_input, review_path), paths,
    vapply(packages, requireNamespace, logical(1), quietly = TRUE))
  cat("Reviewed clusters: ", nrow(review), "; input cells: ", sum(review$expected_cells),
      "; removed: ", sum(review$expected_cells[review$action == "remove"]),
      "; retained: ", sum(review$expected_cells[review$action == "retain"]), "\n", sep = "")
  cat("Removed labels: Brainstem; RBC; Other\nOutput: ", root, "\n", sep = "")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!all(file.exists(c(input_path, marker_input)))) stop("Missing source Science object or marker evidence.")
if (any(file.exists(paths)) && !"--overwrite" %in% options_science) stop("Reviewed Science outputs exist; review before --overwrite.")
suppressPackageStartupMessages({library(Seurat); library(harmony); library(dplyr); library(ggplot2)})
source(here::here("scripts", "color_palette.R"))
source(here::here("scripts", "R", "dotplot_helpers.R"))
check_science_reviewed_markers(read.csv(marker_input), review)
message(Sys.time(), ": loading reviewed Science source: ", input_path)
pbmc <- readRDS(input_path)
metadata <- pbmc[[]]
if (!identical(rownames(metadata), colnames(pbmc)) || anyDuplicated(colnames(pbmc))) stop("Unaligned/duplicate Science cell IDs.")
if (!"originalexp" %in% Assays(pbmc) || !"umap" %in% Reductions(pbmc)) stop("Source must have originalexp assay and original UMAP.")
cell_mapping <- map_science_reviewed_metadata(metadata, review)
label_order <- c("RL", "UBC", "Granule", "eCN", "VZ", "Purkinje", "GABA", "Glia", "OPC",
                 "Meninges", "Vascular", "Immune", "Brainstem", "RBC", "Other")
if (!all(cell_mapping$broad_label %in% label_order)) stop("Unrecognized reviewed Science labels.")
if ("science_broad_label" %in% names(metadata)) {
  if ("science_broad_label_before_review" %in% names(metadata)) stop("Prior-review backup column already exists; wrong source object?")
  pbmc$science_broad_label_before_review <- as.character(pbmc$science_broad_label)
}
pbmc$science_broad_label <- factor(cell_mapping$broad_label, levels = label_order)
Idents(pbmc) <- "science_broad_label"
label_colors <- cluster_colors
label_colors["Vascular"] <- cluster_colors[["Vascular"]]
label_colors[c("Brainstem", "RBC", "Other")] <- c("#A6761D", "#B2182B", "grey70")
source_coordinates <- Embeddings(pbmc, "umap")
if (any(!is.finite(source_coordinates)) || !identical(rownames(source_coordinates), colnames(pbmc))) stop("Invalid source UMAP coordinates.")
fixed_x <- range(source_coordinates[, 1]); fixed_y <- range(source_coordinates[, 2])
umap_review <- function(object, title, original_coordinates = FALSE) {
  plot <- DimPlot(object, reduction = "umap", group.by = "science_broad_label", label = TRUE,
                  repel = TRUE, cols = label_colors) + ggtitle(title)
  if (original_coordinates) plot <- plot + coord_cartesian(xlim = fixed_x, ylim = fixed_y)
  plot
}
dir.create(root, recursive = TRUE, showWarnings = FALSE)
ggsave(paths[["before"]], umap_review(pbmc, "Reviewed Science labels: before removal", TRUE), width = 10, height = 8)
keep <- cell_mapping$cell_id[cell_mapping$action == "retain"]
source_sample_ids <- unique(cell_mapping$sample_id)
pbmc <- subset(pbmc, cells = keep)
if (!identical(colnames(pbmc), keep) || !setequal(unique(as.character(pbmc$sample_id)), source_sample_ids)) {
  stop("Retained Science IDs/order or represented sample set changed unexpectedly.")
}
pbmc$science_broad_label <- droplevels(pbmc$science_broad_label)
if (any(as.character(pbmc$science_broad_label) %in% c("Brainstem", "RBC", "Other"))) stop("Excluded Science labels remain.")
retained_metadata <- pbmc[[]]
ggsave(paths[["after_removal"]], umap_review(pbmc, "After removal: original UMAP coordinates", TRUE), width = 10, height = 8)
# Repeat the established originalexp normalization, cycle scoring/regression,
# PCA-50, Harmony(sample_id; dims 1:40), neighbors, res2 clusters, UMAP dims1:30.
# No additional QC thresholds, stage filtering, or package changes are introduced.
science_resolution <- 2
science_seed <- 42L
set.seed(science_seed)
DefaultAssay(pbmc) <- "originalexp"
pbmc <- NormalizeData(pbmc, verbose = FALSE)

## 1) cell cycle score（using updated.2019）
if (!exists("cc.genes.updated.2019")) {
  suppressWarnings(suppressMessages(
    try(data("cc.genes.updated.2019", package = "Seurat"), silent = TRUE)
  ))
}
if (exists("cc.genes.updated.2019")) {
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
} else {
  suppressWarnings(suppressMessages(
    try(data("cc.genes", package = "Seurat"), silent = TRUE)
  ))
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
}

## 2) using the genes in originalexp
feat <- rownames(pbmc[["originalexp"]])
s.use <- intersect(s.genes, feat)
g.use <- intersect(g2m.genes, feat)

## match using SYMBOL → ENSEMBL 
if ((length(s.use) < 10 || length(g.use) < 10) &&
    any(grepl("^ENSG", feat)) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  s.map <- AnnotationDbi::select(org.Hs.eg.db, keys = s.genes,
                                 keytype = "SYMBOL", columns = "ENSEMBL")
  g.map <- AnnotationDbi::select(org.Hs.eg.db, keys = g2m.genes,
                                 keytype = "SYMBOL", columns = "ENSEMBL")
  s.use <- intersect(unique(na.omit(s.map$ENSEMBL)), feat)
  g.use <- intersect(unique(na.omit(g.map$ENSEMBL)), feat)
}

cat("S genes used:", length(s.use), " | G2M genes used:", length(g.use), "\n")

## 3) cell cycle score（ident）
pbmc <- CellCycleScoring(
  pbmc,
  s.features  = s.use,
  g2m.features = g.use,
  set.ident   = FALSE,
  assay       = "originalexp"
)

pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score","G2M.Score"))
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)

dims_use <- 1:40
if (!"pca" %in% Reductions(pbmc)) {
  pbmc <- RunPCA(pbmc, npcs = max(dims_use), verbose = FALSE)
}

# --- using RunHarmony ---
rh <- getS3method("RunHarmony", "Seurat")  
args <- names(formals(rh))
parms <- list(object = pbmc, group.by.vars = "sample_id")

if ("reduction.use" %in% args) parms$reduction.use <- "pca"
if (!("reduction.use" %in% args) && "reduction" %in% args) parms$reduction <- "pca"
if ("dims.use" %in% args) parms$dims.use <- dims_use
if ("assay.use" %in% args) parms$assay.use <- DefaultAssay(pbmc)
if ("project.dim" %in% args) parms$project.dim <- FALSE
if ("verbose" %in% args) parms$verbose <- TRUE

pbmc <- do.call(harmony::RunHarmony, parms)

# --- use harmony UMAP ---
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = dims_use)

# for SNN ）
avail <- c(names(pbmc@graphs), names(pbmc@neighbors))
snn_name <- grep("snn", avail, value = TRUE)[1]

pbmc <- FindClusters(pbmc, graph.name = snn_name, resolution = science_resolution,
                     random.seed = science_seed)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:30, verbose = FALSE)

# Keep the approved labels attached to cells; new numeric IDs never inherit the old map.
pbmc$science_cluster_res2_postreview <- as.character(pbmc$seurat_clusters)
# FindClusters only intentionally changes seurat_clusters and its graph-resolution column.
# All original IDs, sample IDs, source raw clusters, and other metadata must survive.
protected_columns <- setdiff(names(retained_metadata), c("seurat_clusters", "originalexp_snn_res.2", "S.Score", "G2M.Score", "Phase"))
for (column in protected_columns) {
  if (!identical(pbmc[[]][[column]], retained_metadata[[column]])) stop("Reprocessing changed protected metadata: ", column)
}
if (!identical(colnames(pbmc), keep)) stop("Reprocessing changed retained Science cell IDs/order.")
Idents(pbmc) <- "science_broad_label"
ggsave(paths[["after_reprocessing"]], umap_review(pbmc, "Retained Science: reprocessed UMAP"), width = 10, height = 8)
ggsave(paths[["raw_umap"]], DimPlot(pbmc, group.by = "science_cluster_res2_postreview", label = TRUE) +
  ggtitle("Retained Science: new resolution-2 clusters"), width = 10, height = 8)
marker_sets <- lapply(markers, function(genes) intersect(genes, rownames(pbmc[["originalexp"]])))
marker_sets <- marker_sets[lengths(marker_sets) > 0]
if (!length(marker_sets)) stop("No broad marker symbols available after reprocessing.")
ggsave(paths[["dotplot"]], workflow_dotplot(pbmc, features = marker_sets, assay = "originalexp",
  group.by = "science_broad_label") + RotatedAxis() + ggtitle("Retained Science: reviewed broad labels"),
  width = 16, height = 9, limitsize = FALSE)
message(Sys.time(), ": finding post-review raw-cluster markers")
Idents(pbmc) <- "science_cluster_res2_postreview"
marker_results <- FindAllMarkers(pbmc, assay = "originalexp", only.pos = TRUE,
  min.pct = 0.25, logfc.threshold = 0.25, random.seed = science_seed)
required_marker_columns <- c("cluster", "gene", "avg_log2FC")
if (!nrow(marker_results) || !all(required_marker_columns %in% names(marker_results)) ||
    anyNA(marker_results[, required_marker_columns]) || any(!is.finite(marker_results$avg_log2FC)) ||
    anyDuplicated(marker_results[, c("cluster", "gene")])) stop("Invalid post-review marker results; inspect warnings.")
new_clusters <- sort(unique(as.character(pbmc$science_cluster_res2_postreview)))
if (!setequal(as.character(marker_results$cluster), new_clusters)) stop("Incomplete post-review marker cluster set.")
top20 <- marker_results %>% group_by(cluster) %>% arrange(desc(avg_log2FC), gene, .by_group = TRUE) %>%
  slice_head(n = 20L) %>% mutate(marker_rank = row_number()) %>% ungroup()
marker_counts <- data.frame(cluster = new_clusters,
  cells = as.integer(table(factor(pbmc$science_cluster_res2_postreview, levels = new_clusters))),
  exported_markers = as.integer(table(factor(top20$cluster, levels = new_clusters))))
if (any(marker_counts$exported_markers < 20)) warning("Some post-review clusters have fewer than 20 qualifying markers; see marker counts.")
cell_mapping$postreview_cluster <- NA_character_
cell_mapping$postreview_cluster[match(colnames(pbmc), cell_mapping$cell_id)] <- as.character(pbmc$science_cluster_res2_postreview)
retention <- review[, c("cluster", "broad_label", "expected_cells", "action")]
retention$retained_cells <- ifelse(retention$action == "retain", retention$expected_cells, 0L)
retention$removed_cells <- retention$expected_cells - retention$retained_cells
sample_counts <- aggregate(list(cells = rep(1L, nrow(cell_mapping))),
  cell_mapping[, c("sample_id", "broad_label", "action")], sum)
if (sum(retention$retained_cells) != ncol(pbmc) || sum(sample_counts$cells) != nrow(cell_mapping)) stop("Science retention counts do not reconcile.")
write.csv(cell_mapping, paths[["cells"]], row.names = FALSE)
write.csv(retention, paths[["retention"]], row.names = FALSE)
write.csv(sample_counts, paths[["sample_counts"]], row.names = FALSE)
write.csv(top20, paths[["markers"]], row.names = FALSE)
write.csv(marker_counts, paths[["marker_counts"]], row.names = FALSE)
Idents(pbmc) <- "science_broad_label"
saveRDS(pbmc, paths[["rds"]])
write.csv(data.frame(input_rds = normalizePath(input_path), input_size = file.info(input_path)$size,
  input_mtime = as.character(file.info(input_path)$mtime), review_csv = normalizePath(review_path),
  review_md5 = unname(tools::md5sum(review_path)), marker_evidence_md5 = unname(tools::md5sum(marker_input)),
  source_cells = nrow(cell_mapping), retained_cells = ncol(pbmc), removed_cells = sum(retention$removed_cells),
  removed_labels = "Brainstem;RBC;Other", samples = length(source_sample_ids),
  source_cluster_column = "science_cluster_res2", new_cluster_column = "science_cluster_res2_postreview",
  normalization = "existing NormalizeData defaults", cycle_regression = "S.Score;G2M.Score", pca_npcs = 50,
  harmony_group = "sample_id", harmony_neighbor_dims = "1:40", umap_dims = "1:30", resolution = 2,
  seed = science_seed, marker_min_pct = 0.25, marker_logfc_threshold = 0.25, marker_top_n = 20,
  output_rds = paths[["rds"]], script_md5 = unname(tools::md5sum(here::here("scripts", "science_03_map_reviewed_res2.R"))),
  created = as.character(Sys.time())), paths[["provenance"]], row.names = FALSE)
writeLines(capture.output(sessionInfo()), paths[["session"]])
if (any(is.na(file.info(paths)$size) | file.info(paths)$size <= 0)) stop("Incomplete reviewed Science outputs.")
message("Completed reviewed Science mapping/removal/reprocessing: ", root)
