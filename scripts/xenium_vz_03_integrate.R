# Clear the environment
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- samples$sample_id

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--dry-run", "--overwrite", "--plots-only", "--include-uncorrected",
  "--resume-checkpoint"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_vz_03_integrate.R ",
    "[--dry-run] [--overwrite] [--plots-only] [--include-uncorrected] ",
    "[--resume-checkpoint]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
plots_only <- "--plots-only" %in% args
include_uncorrected <- "--include-uncorrected" %in% args
resume_checkpoint <- "--resume-checkpoint" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (include_uncorrected && !plots_only) {
  stop("--include-uncorrected requires --plots-only.")
}
if (resume_checkpoint && plots_only) {
  stop("--resume-checkpoint and --plots-only cannot be combined.")
}

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
output_path <- file.path(rds_dir, "Xenium_VZ_integrated.rds")
checkpoint_path <- file.path(rds_dir, "Xenium_VZ_precluster_checkpoint.rds")
batch_plot_path <- file.path(plot_dir, "Xenium_VZ_Batch_Comp_UMAP.tif")
plot_paths <- c(
  file.path(plot_dir, paste0("Xenium_VZ_UMAP_Res_", res_list, ".tif")),
  file.path(plot_dir, "Xenium_VZ_ConsensusLabel_UMAP.tif"),
  file.path(plot_dir, "Xenium_VZ_UnknownClusters_BySample_UMAP.tif")
)
requested_plot_paths <- c(plot_paths, if (include_uncorrected) batch_plot_path)
expected_outputs <- c(output_path, checkpoint_path, plot_paths)

manifest_ready <- FALSE
manifest_rows <- NA_integer_
manifest_ids_match <- FALSE
manifest_columns_match <- FALSE
manifest_values_valid <- FALSE
if (file.exists(merge_manifest_path)) {
  merge_manifest_preview <- read.csv(merge_manifest_path, stringsAsFactors = FALSE)
  manifest_rows <- nrow(merge_manifest_preview)
  manifest_columns_match <- all(c(
    "sample_id", "input_path", "cells", "configured_broad_cells",
    "shared_unknown_cells", "sample_qualified_unknown_labels", "PCW"
  ) %in% names(merge_manifest_preview))
  manifest_ids_match <- "sample_id" %in% names(merge_manifest_preview) &&
    setequal(as.character(merge_manifest_preview$sample_id), sample_ids)
  if (manifest_columns_match) {
    manifest_cells <- suppressWarnings(as.numeric(merge_manifest_preview$cells))
    manifest_pcw <- as.character(merge_manifest_preview$PCW)
    manifest_broad <- suppressWarnings(as.numeric(merge_manifest_preview$configured_broad_cells))
    manifest_unknown <- suppressWarnings(as.numeric(merge_manifest_preview$shared_unknown_cells))
    manifest_qualified <- as.character(merge_manifest_preview$sample_qualified_unknown_labels)
    manifest_values_valid <- all(is.finite(manifest_cells) & manifest_cells > 0) &&
      all(is.finite(manifest_broad) & manifest_broad >= 0) &&
      all(is.finite(manifest_unknown) & manifest_unknown >= 0) &&
      all(manifest_broad + manifest_unknown == manifest_cells) &&
      all(
        (manifest_unknown > 0 & !is.na(manifest_qualified) & nzchar(manifest_qualified)) |
          (manifest_unknown == 0 & (is.na(manifest_qualified) | !nzchar(manifest_qualified)))
      ) &&
      all(!is.na(manifest_pcw) & nzchar(manifest_pcw))
  }
  manifest_ready <- manifest_rows == length(sample_ids) && manifest_columns_match &&
    manifest_ids_match && manifest_values_valid
}

if (dry_run) {
  if (plots_only) {
    compact_dry_run(
      "VZ integration plots only",
      inputs = output_path,
      outputs = requested_plot_paths
    )
    quit(save = "no", status = 0L)
  }
  if (resume_checkpoint) {
    compact_dry_run(
      "VZ integration checkpoint resume",
      inputs = c(checkpoint_path, merged_path, merge_manifest_path),
      outputs = c(output_path, plot_paths),
      checks = c(
        manifest_rows = manifest_rows == length(sample_ids),
        manifest_columns = manifest_columns_match,
        manifest_values = manifest_values_valid,
        manifest_sample_ids = manifest_ids_match
      )
    )
    quit(save = "no", status = 0L)
  }
  compact_dry_run(
    "VZ integration",
    inputs = c(merged_path, merge_manifest_path),
    outputs = expected_outputs,
    checks = c(
      manifest_rows = manifest_rows == length(sample_ids),
      manifest_columns = manifest_columns_match,
      manifest_values = manifest_values_valid,
      manifest_sample_ids = manifest_ids_match
    )
  )
  quit(save = "no", status = 0L)
}

if (plots_only) {
  if (!file.exists(output_path)) stop("Integrated VZ input not found: ", output_path)
  existing_outputs <- requested_plot_paths[file.exists(requested_plot_paths)]
} else {
  if (!file.exists(merged_path)) stop("Merged VZ input not found: ", merged_path)
  if (!manifest_ready) {
    stop("VZ merge manifest is missing, incomplete, or inconsistent: ", merge_manifest_path)
  }
  if (resume_checkpoint && !file.exists(checkpoint_path)) {
    stop("VZ pre-clustering checkpoint not found: ", checkpoint_path)
  }
  candidate_outputs <- if (resume_checkpoint) c(output_path, plot_paths) else expected_outputs
  existing_outputs <- candidate_outputs[file.exists(candidate_outputs)]
}
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ ",
    if (plots_only) "plot" else "processing",
    " outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}
if (!plots_only && file.exists(batch_plot_path)) {
  warning(
    "The existing deferred batch-comparison plot is not regenerated by full stage 3 ",
    "and may be stale until --plots-only --include-uncorrected completes: ",
    batch_plot_path
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
source(here("scripts", "color_palette.R"))

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
workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
cluster_workers <- suppressWarnings(as.integer(
  Sys.getenv("XENIUM_STAGE3_CLUSTER_WORKERS", unset = "1")
))
if (is.na(cluster_workers) || cluster_workers < 1L) {
  stop("XENIUM_STAGE3_CLUSTER_WORKERS must be a positive integer when set.")
}
cluster_workers <- min(cluster_workers, workers, length(res_list))
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(config$runtime$random_seed)

with_fork_workers <- function(expr, max_workers = workers) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  active_workers <- min(workers, max_workers)
  if (active_workers > 1L) {
    if (!future::supportsMulticore()) {
      stop(
        "Fork-based parallelism is unavailable. Resubmit with ",
        "--cpus-per-task=1 to avoid multisession copies of the large object."
      )
    }
    plan(multicore, workers = active_workers)
  }
  force(expr)
}

check_mem("PIPELINE START")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

if (plots_only) {
  message("Plot-only mode: loading the existing integrated VZ object.")
  obj <- readRDS(output_path)
  required_plot_metadata <- c(
    "orig.ident", "consensus_label", "regional_consensus_label",
    "regional_consensus_label_pooled", "regional_unknown_label"
  )
  missing_plot_metadata <- setdiff(required_plot_metadata, colnames(obj[[]]))
  if (length(missing_plot_metadata)) {
    stop("Integrated VZ input lacks plot metadata: ", paste(missing_plot_metadata, collapse = ", "))
  }
  if (!setequal(unique(as.character(obj$orig.ident)), sample_ids)) {
    stop("Integrated VZ input does not contain exactly the configured sample IDs.")
  }
  if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
    stop("Integrated VZ sample-qualified Unknown labels are invalid.")
  }
  check_mem("INTEGRATED OBJECT LOADED FOR PLOTTING")
} else if (resume_checkpoint) {
  message("Resume mode: loading the VZ pre-clustering checkpoint.")
  obj <- readRDS(checkpoint_path)
  provenance <- attr(obj, "xenium_stage3_checkpoint")
  merged_info <- file.info(merged_path)
  checkpoint_valid <- is.list(provenance) &&
    identical(provenance$region, "VZ") &&
    identical(provenance$sample_ids, sample_ids) &&
    identical(provenance$merged_path, normalizePath(merged_path)) &&
    identical(provenance$merged_size, as.numeric(merged_info$size)) &&
    identical(provenance$merged_mtime, as.numeric(merged_info$mtime)) &&
    identical(provenance$random_seed, config$runtime$random_seed) &&
    identical(provenance$resolutions, res_list)
  if (!checkpoint_valid) {
    stop("VZ checkpoint provenance does not match the current merge, manifest, or parameters.")
  }
  if (!all(c("pca", "harmony", "umap_harmony") %in% Reductions(obj))) {
    stop("VZ checkpoint lacks PCA, Harmony, or corrected UMAP reductions.")
  }
  if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
    stop("VZ checkpoint sample-qualified Unknown labels are invalid.")
  }
  check_mem("PRE-CLUSTERING CHECKPOINT LOADED")
} else {
  # 3. LOAD DATA
  obj <- readRDS(merged_path)
  required_metadata <- c(
    "orig.ident", "consensus_label", "regional_consensus_label",
    "regional_consensus_label_pooled", "regional_unknown_label",
    "VZ_subset_reason", "PCW"
  )
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
  if (anyNA(obj$VZ_subset_reason) || any(!obj$VZ_subset_reason %in%
      c("configured_broad_label", "shared_unknown"))) {
    stop("Merged VZ input has missing or invalid VZ_subset_reason values.")
  }
  observed_unknown <- vapply(sample_ids, function(sample_id) {
    sample_cells <- as.character(obj$orig.ident) == sample_id
    sum(obj$VZ_subset_reason[sample_cells] == "shared_unknown")
  }, integer(1))
  expected_unknown <- as.integer(merge_manifest_preview$shared_unknown_cells[manifest_order])
  if (!identical(unname(observed_unknown), unname(expected_unknown)) || !all(
    (obj$VZ_subset_reason == "shared_unknown") ==
      is_unknown_consensus_label(obj$consensus_label)
  )) {
    stop("Merged VZ shared-Unknown provenance does not match its manifest or labels.")
  }
  if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
    stop("Merged VZ sample-qualified Unknown labels are invalid.")
  }
  observed_qualified <- vapply(sample_ids, function(sample_id) {
    sample_cells <- as.character(obj$orig.ident) == sample_id
    paste(sort(unique(as.character(obj$regional_consensus_label[
      sample_cells & is_unknown_consensus_label(obj$consensus_label)
    ]))), collapse = "|")
  }, character(1))
  expected_qualified <- as.character(
    merge_manifest_preview$sample_qualified_unknown_labels[manifest_order]
  )
  expected_qualified[is.na(expected_qualified)] <- ""
  if (!identical(unname(observed_qualified), unname(expected_qualified))) {
    stop("Merged VZ sample-qualified Unknown labels do not match the merge manifest.")
  }
  Idents(obj) <- "regional_consensus_label"
  check_mem("DATA LOADED")

  # 4. STANDARD WORKFLOW
  message("Normalizing and Finding Variable Features...")
  obj <- NormalizeData(obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  message("Scaling data...")
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
  check_mem("POST-SCALE")

  message("Running PCA...")
  obj <- RunPCA(
    obj, npcs = 30, verbose = FALSE,
    seed.use = config$runtime$random_seed
  )

  message(
    "Deferring the uncorrected diagnostic UMAP. Generate it later with ",
    "--plots-only --include-uncorrected."
  )

  # 5. RUN HARMONY
  message("Running Harmony integration...")
  obj <- RunHarmony(obj, group.by.vars = "orig.ident",
                    dims.use = 1:30, reduction.save = "harmony", verbose = TRUE)
  check_mem("POST-HARMONY")

  # 6. POST-HARMONY REDUCTION & NEIGHBORS
  message("Running Integrated UMAP and Neighbors...")
  obj <- with_fork_workers(
    RunUMAP(obj, reduction = "harmony", dims = 1:30,
            reduction.name = "umap_harmony", verbose = TRUE,
            seed.use = config$runtime$random_seed),
    max_workers = workers
  )

  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = TRUE)
  attr(obj, "xenium_stage3_checkpoint") <- list(
    region = "VZ",
    sample_ids = sample_ids,
    merged_path = normalizePath(merged_path),
    merged_size = as.numeric(file.info(merged_path)$size),
    merged_mtime = as.numeric(file.info(merged_path)$mtime),
    random_seed = config$runtime$random_seed,
    resolutions = res_list,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  message("Saving restart checkpoint: ", checkpoint_path)
  saveRDS(obj, checkpoint_path, compress = FALSE)
  check_mem("PRE-CLUSTERING CHECKPOINT SAVED")
}

if (!plots_only) {
  # 7. GENERATE CLUSTERS. This defaults to one process because the 1.25-million-
  # cell VZ graph exceeded 256 GB when forked across three workers (job 40676165).
  message(
    Sys.time(), ": Calculating all clusters with ", cluster_workers,
    if (cluster_workers == 1L) " sequential worker..." else " fork-based workers..."
  )
  obj <- with_fork_workers(
    FindClusters(obj, resolution = res_list, verbose = FALSE),
    max_workers = cluster_workers
  )
  check_mem("POST-BATCH-CLUSTERING")
}

assay_prefix <- DefaultAssay(obj)
required_reductions <- "umap_harmony"
missing_reductions <- setdiff(required_reductions, Reductions(obj))
if (length(missing_reductions)) {
  stop("Integrated VZ object lacks reductions: ", paste(missing_reductions, collapse = ", "))
}
required_resolution_columns <- paste0(assay_prefix, "_snn_res.", res_list)
missing_resolution_columns <- setdiff(required_resolution_columns, colnames(obj[[]]))
if (length(missing_resolution_columns)) {
  stop(
    "Integrated VZ object lacks resolution metadata: ",
    paste(missing_resolution_columns, collapse = ", ")
  )
}

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
    filename = file.path(plot_dir, paste0("Xenium_VZ_UMAP_Res_", res, ".tif")),
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

if (include_uncorrected) {
  message("Generating the deferred uncorrected UMAP and batch comparison...")
  if (!"pca" %in% Reductions(obj)) stop("Integrated VZ object lacks the PCA reduction.")
  if (!"umap_uncorrected" %in% Reductions(obj)) {
    obj <- with_fork_workers(
      RunUMAP(obj, reduction = "pca", dims = 1:30,
              reduction.name = "umap_uncorrected", verbose = TRUE,
              seed.use = config$runtime$random_seed),
      max_workers = workers
    )
  }
  p1 <- DimPlot(
    obj, reduction = "umap_uncorrected", group.by = "orig.ident", raster = TRUE
  ) + NoLegend() + ggtitle("Pre-Harmony")
  p2 <- DimPlot(
    obj, reduction = "umap_harmony", group.by = "orig.ident", raster = TRUE
  ) + ggtitle("Post-Harmony")
  Cairo::CairoTIFF(
    filename = batch_plot_path, width = 16, height = 7, units = "in", res = 600
  )
  print(p1 + p2)
  grDevices::dev.off()
  rm(p1, p2)
  gc()
}

# Display each independently numbered per-sample Unknown identity separately,
# while retaining non-Unknown cells as a neutral background category.
unknown_levels <- unlist(lapply(sample_ids, function(sample_id) {
  values <- unique(as.character(obj$regional_consensus_label[
    as.character(obj$orig.ident) == sample_id &
      is_unknown_consensus_label(obj$consensus_label)
  ]))
  sort(values)
}), use.names = FALSE)
unknown_levels <- unique(unknown_levels[nzchar(unknown_levels)])
if (!length(unknown_levels)) stop("VZ integration contains no sample-qualified Unknown labels.")
obj$regional_unknown_label <- factor(
  as.character(obj$regional_unknown_label),
  levels = c("Non-Unknown", unknown_levels)
)
unknown_colors <- resolve_sample_qualified_unknown_colors(
  c("Non-Unknown", unknown_levels)
)
p_unknown_by_sample <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by = "regional_unknown_label",
  cols = unknown_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  label.box = TRUE,
  raster = TRUE,
  pt.size = 0.5,
  alpha = 0.8
) +
  ggtitle("Xenium VZ Integrated UMAP: Sample-qualified Unknown clusters") +
  theme_classic() +
  guides(
    color = guide_legend(
      ncol = 2,
      byrow = FALSE,
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.spacing.y = grid::unit(0.03, "cm")
  )

Cairo::CairoTIFF(
  filename = file.path(plot_dir, "Xenium_VZ_UnknownClusters_BySample_UMAP.tif"),
  width = 16,
  height = 14,
  units = "in",
  res = 600
)
print(p_unknown_by_sample)
grDevices::dev.off()

# Companion overview: preserve known consensus identities and combine every
# sample-specific Unknown identity into one explicitly grey category.
pooled_levels <- order_consensus_levels(
  obj$regional_consensus_label_pooled, celltype_order
)
obj$regional_consensus_label_pooled <- factor(
  as.character(obj$regional_consensus_label_pooled), levels = pooled_levels
)
pooled_colors <- resolve_cluster_colors(pooled_levels)
p_unknown_pooled <- DimPlot(
  obj,
  reduction = "umap_harmony",
  group.by = "regional_consensus_label_pooled",
  cols = pooled_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  label.box = TRUE,
  raster = TRUE,
  pt.size = 0.5,
  alpha = 0.8
) +
  ggtitle("Xenium VZ Integrated UMAP: Unknown pooled in grey") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

Cairo::CairoTIFF(
  filename = file.path(plot_dir, "Xenium_VZ_ConsensusLabel_UMAP.tif"),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p_unknown_pooled)
grDevices::dev.off()

if (!plots_only) {
  saveRDS(obj, output_path, compress = FALSE)
}

plan("sequential")

if (plots_only) {
  check_mem("PLOT-ONLY MODE COMPLETE")
  message("Plot-only mode replaced the VZ stage-3 plots and did not modify ", output_path, ".")
} else {
  check_mem("PIPELINE COMPLETE - FILE SAVED")
  message(
    "Stage 3 integration does not run marker testing. Run the explicit VZ stage-4 ",
    "QC review to create the QC violin PDF, cluster summary, and top-10 marker table."
  )
  message("Restart checkpoint retained at ", checkpoint_path, ".")
  message("The pre-/post-Harmony comparison remains deferred to plot-only mode.")
}
