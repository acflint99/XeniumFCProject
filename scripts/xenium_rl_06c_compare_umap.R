#!/usr/bin/env Rscript
# Visualization-only comparison: fit all cells; use a fixed balanced display subset.
# Never writes a Seurat object or changes approved labels/clusters.
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
render_selected <- "--render-selected" %in% args
upstream_only <- "--upstream-only" %in% args
if (upstream_only && !render_selected) stop("--upstream-only requires --render-selected.")
args <- args[!args %in% c("--render-selected", "--upstream-only")]
if (any(!args %in% c("--list", "--dry-run", "--overwrite"))) {
  stop("Usage: xenium_rl_06c_compare_umap.R [--render-selected [--upstream-only]] [--list|--dry-run|--overwrite]")
}
if (length(args) > 1L) stop("Choose only one option.")
settings <- data.frame(
  comparison = c("nn30_md0.10", "nn30_md0.05", "nn15_md0.10", "nn15_md0.05"),
  n_neighbors = c(30L, 30L, 15L, 15L),
  min_dist = c(0.10, 0.05, 0.10, 0.05)
)
if ("--list" %in% args) {
  print(settings, row.names = FALSE)
  quit(save = "no")
}
input <- here(config$project$outputs_dir, "xenium", "rl", "06_subclusters",
              "rds", "Xenium_RL_subclusters.rds")
out_dir <- here(config$project$outputs_dir, "xenium", "rl", "06c_umap_comparison")
if (render_selected) {
  source(here("scripts", "R", "rl_umap_selected.R"))
  render_selected_rl_umap(config, input, out_dir,
                          dry_run = "--dry-run" %in% args,
                          overwrite = "--overwrite" %in% args,
                          upstream_only = upstream_only)
  quit(save = "no")
}
pdf_path <- file.path(out_dir, "Xenium_RL_UMAP_Parameter_Comparison.pdf")
manifest_path <- file.path(out_dir, "Xenium_RL_UMAP_Comparison_manifest.csv")
display_path <- file.path(out_dir, "Xenium_RL_UMAP_display_cells.csv")
coordinate_paths <- file.path(out_dir, paste0("Xenium_RL_", settings$comparison, "_coordinates.rds"))
session_path <- file.path(out_dir, "sessionInfo.txt")
outputs <- c(pdf_path, manifest_path, display_path, coordinate_paths, session_path)
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("RL UMAP parameter comparison (paths only)", input, outputs)
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!file.exists(input)) stop("Missing approved annotated object: ", input)
existing <- outputs[file.exists(outputs)]
if (length(existing) && !"--overwrite" %in% args) {
  stop("Refusing to overwrite comparison outputs:\n", paste(existing, collapse = "\n"),
       "\nUse --overwrite only after review.")
}
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source(here("scripts", "color_palette.R"))
seed <- as.integer(config$runtime$random_seed)
dims <- seq_len(as.integer(config$regional_subsets$post_qc_dimensions_observed))
if (length(seed) != 1L || is.na(seed) || !length(dims)) stop("Invalid seed or dimensions.")
input_info <- file.info(input)
message("Reading approved object: ", input)
obj <- readRDS(input)
if (!all(c("harmony_clean", "umap_clean") %in% Reductions(obj))) {
  stop("Required reductions: harmony_clean and umap_clean.")
}
if (!"RL_subcluster" %in% colnames(obj[[]])) stop("Missing RL_subcluster metadata.")
cells <- colnames(obj)
if (!length(cells) || anyDuplicated(cells)) stop("Empty or duplicate cell IDs.")
meta <- obj[[]]
if (!identical(rownames(meta), cells)) stop("Metadata/cell order mismatch.")
labels <- as.character(meta$RL_subcluster)
if (anyNA(labels) || any(!nzchar(labels)) || any(!labels %in% names(rl_palette))) {
  stop("Missing, blank, or unregistered RL labels.")
}
latent <- Embeddings(obj, "harmony_clean")
baseline <- Embeddings(obj, "umap_clean")
if (!identical(rownames(latent), cells) || !identical(rownames(baseline), cells)) {
  stop("Reduction/cell order mismatch.")
}
if (ncol(latent) < max(dims) || ncol(baseline) != 2L) stop("Unexpected reduction dimensions.")
latent <- latent[, dims, drop = FALSE]
if (any(!is.finite(latent)) || any(!is.finite(baseline))) stop("Non-finite coordinates.")
rm(obj, meta)
gc()
levels_used <- rl_subcluster_order[rl_subcluster_order %in% labels]
if (!setequal(levels_used, unique(labels))) stop("RL label order is incomplete.")
# Display at most 2,000 cells per identity; all cells participate in every fit.
set.seed(seed)
display_indices <- unlist(lapply(levels_used, function(label) {
  ids <- which(labels == label)
  if (length(ids) > 2000L) ids <- ids[sample.int(length(ids), 2000L)]
  ids
}), use.names = FALSE)
display_indices <- display_indices[sample.int(length(display_indices))]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(cell_id = cells[display_indices],
                     RL_subcluster = labels[display_indices]), display_path, row.names = FALSE)
plot_embedding <- function(coords, title) {
  d <- data.frame(UMAP_1 = coords[display_indices, 1],
                  UMAP_2 = coords[display_indices, 2],
                  RL_subcluster = factor(labels[display_indices], levels = levels_used))
  ggplot(d, aes(UMAP_1, UMAP_2, color = RL_subcluster)) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_manual(values = rl_palette, drop = FALSE) +
    coord_equal() + theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = title, subtitle = paste0("All ", length(cells),
      " cells fitted; same display subset, up to 2,000 cells per identity"), color = NULL)
}
run_comparisons <- function() {
  grDevices::pdf(pdf_path, width = 12, height = 9, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_embedding(baseline, "Existing approved UMAP (coordinates unchanged)"))
  records <- vector("list", nrow(settings))
  for (i in seq_len(nrow(settings))) {
    p <- settings[i, ]
    message(format(Sys.time()), " | Starting ", p$comparison)
    start <- Sys.time()
    fit <- RunUMAP(latent, n.neighbors = p$n_neighbors, min.dist = p$min_dist,
                      spread = 1, metric = "cosine", umap.method = "uwot",
                      seed.use = seed, verbose = TRUE)
    coords <- Embeddings(fit)
    if (!is.matrix(coords) || !identical(dim(coords), c(length(cells), 2L)) ||
        any(!is.finite(coords)) || !identical(rownames(coords), cells)) {
      stop("Invalid UMAP output or changed cell order for ", p$comparison)
    }
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    saveRDS(data.frame(cell_id = cells, RL_subcluster = labels,
                       UMAP_1 = coords[, 1], UMAP_2 = coords[, 2]), coordinate_paths[i])
    print(plot_embedding(coords, paste0("n.neighbors = ", p$n_neighbors,
                                        "; min.dist = ", p$min_dist)))
    records[[i]] <- data.frame(p, input = input, input_bytes = input_info$size,
      input_mtime = as.character(input_info$mtime), n_cells = length(cells),
      reduction = "harmony_clean", dimensions = paste(dims, collapse = "|"),
      seed = seed, spread = 1, metric = "cosine", umap_method = "uwot",
      seurat_version = as.character(packageVersion("Seurat")),
      uwot_version = as.character(packageVersion("uwot")),
      display_n = length(display_indices), fit_minutes = elapsed,
      coordinates = coordinate_paths[i], comparison_pdf = pdf_path)
    rm(coords, fit)
    gc()
  }
  do.call(rbind, records)
}
manifest <- run_comparisons()
write.csv(manifest, manifest_path, row.names = FALSE)
capture.output(sessionInfo(), file = session_path)
message("Completed four comparisons. Review: ", pdf_path)
