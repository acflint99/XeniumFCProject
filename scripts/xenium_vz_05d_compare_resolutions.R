#!/usr/bin/env Rscript
# Protected full-data comparison. Reuse saved 0.5/0.8 labels; calculate 0.6/0.7
# on the identical saved SNN, with production Louvain parameters and seed.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_vz_05d_compare_resolutions.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "regional_cluster_stability.R"))
source(here("scripts", "R", "regional_resolution_comparison.R"))
source(here("scripts", "R", "vz_resolution_efficiency.R"))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  flags <- c("--list", "--dry-run", "--check-packages", "--overwrite", "--resume-checkpoint", "--tiff")
  if (any(!args %in% flags & !grepl("^--run-id=[A-Za-z0-9][A-Za-z0-9_-]*$", args))) {
    stop("Usage: [--list|--dry-run|--check-packages] [--run-id=review1] [--resume-checkpoint] [--tiff] [--overwrite]")
  }
  run_id <- sub("^--run-id=", "", args[startsWith(args, "--run-id=")])
  if (length(run_id) > 1L) stop("Repeated --run-id.")
  if (!length(run_id)) run_id <- "review1"
  if ("--overwrite" %in% args && any(c("--list", "--dry-run", "--check-packages") %in% args)) {
    stop("Inspection and overwrite cannot be combined.")
  }
  resume_checkpoint <- "--resume-checkpoint" %in% args
  render_tiff <- "--tiff" %in% args
  config <- load_pipeline_config()
  samples <- load_sample_manifest(config)$sample_id
  if (!length(samples) || anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) stop("Invalid sample manifest.")
  resolutions <- c(0.5, 0.6, 0.7, 0.8)
  columns <- paste0("VZ_compare_res.", resolutions)
  root <- here(config$project$outputs_dir, "xenium", "vz", "05_post_qc")
  input <- file.path(root, "rds", "Xenium_VZ_postQC.rds")
  input_manifest <- file.path(root, "tables", "Xenium_VZ_postQC_manifest.csv")
  parent <- file.path(root, "diagnostics", "resolution_comparison_05_06_07_08")
  destination <- file.path(parent, run_id)
  checkpoint_path <- file.path(parent, paste0(".", run_id, "_partitions.rds"))
  common <- c("clustree.pdf", "assignments.csv.gz", "cluster_summary.csv", "sample_composition.csv",
    "input_label_composition.csv", "all_pair_transitions.csv", "production0.8_transitions.csv",
    "focus_population_transitions.csv", "provenance.csv", "sessionInfo.txt", "METHODS.txt", "outputs.csv", "partitions_checkpoint.rds", "timings.csv")
  per_resolution <- unlist(lapply(resolutions, function(res) paste0("res", res, "_",
    c("UMAP.pdf", "top50_markers.csv", "marker_status.csv",
      "focused_markers.pdf", if (render_tiff) "focused_markers.tif", "marker_availability.csv"))))
  outputs <- c(common, per_resolution)
  packages <- c("Seurat", "SeuratObject", "Matrix", "future", "ggplot2", "clustree", "ggraph", if (render_tiff) "Cairo", "scattermore")
  if ("--list" %in% args) {
    cat("VZ full-data comparison: 0.5 | 0.6 | 0.7 | 0.8\nManifest samples:", length(samples),
        "\nInput:", input, "\nOutput:", destination, "\nNo array; no production RDS output.\n")
    return(invisible(NULL))
  }
  if ("--dry-run" %in% args) {
    ok <- compact_dry_run("VZ resolution 0.5/0.6/0.7/0.8 comparison (paths only; no RDS load)",
      inputs = c(input, input_manifest), outputs = file.path(destination, outputs),
      checks = c(manifest_valid = TRUE, destination_protected = !file.exists(destination),
        checkpoint_ready = if (resume_checkpoint) file.exists(checkpoint_path) else !file.exists(checkpoint_path)))
    if (!ok) stop("Dry-run failed; inspect missing inputs or choose a fresh run-id.")
    return(invisible(NULL))
  }
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if ("--check-packages" %in% args) print(data.frame(package = packages, available), row.names = FALSE)
  if (!all(available)) stop("Missing dependencies: ", paste(packages[!available], collapse = ", "),
    ". Use an explicit project renv setup before submission; this job never installs packages.")
  regional_stability_check_plot()
  if ("--check-packages" %in% args) return(invisible(NULL))
  if (!all(file.exists(c(input, input_manifest)))) stop("Missing Stage-5 RDS or manifest.")
  if (file.exists(destination) && !"--overwrite" %in% args) stop("Protected output exists: ", destination)
  if (resume_checkpoint && !file.exists(checkpoint_path)) stop("Missing partition checkpoint: ", checkpoint_path)
  if (!resume_checkpoint && file.exists(checkpoint_path) && !"--overwrite" %in% args) {
    stop("Protected checkpoint exists. Use --resume-checkpoint, or choose a fresh --run-id: ", checkpoint_path)
  }
  timing_rows <- list()
  timed <- function(label, expression) {
    started <- proc.time()[[3]]
    value <- force(expression)
    elapsed <- proc.time()[[3]] - started
    timing_rows[[length(timing_rows) + 1L]] <<- data.frame(step = label, elapsed_seconds = elapsed)
    message(label, ": ", round(elapsed, 1), " seconds")
    value
  }
  suppressPackageStartupMessages(library(Seurat))
  source(here("scripts", "R", "subcluster_review_markers.R"))
  source(here("scripts", "color_palette.R"))
  marker_sets <- add_palette_markers_to_review(vz_focused_review_markers, vz_markers, "vz_markers")
  future::plan(future::sequential)
  message("Optional presto backend available: ", requireNamespace("presto", quietly = TRUE),
          ". No packages are installed or changed by this review.")
  seed <- as.integer(config$runtime$random_seed)
  if (length(seed) != 1L || is.na(seed) || seed < 0L) stop("Invalid configured marker/plot seed.")
  set.seed(seed)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(parent, paste0(".lock_", run_id))
  if (!dir.create(lock, showWarnings = FALSE)) stop("Run-id locked; check running/failed jobs: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  provenance <- read.csv(input_manifest, stringsAsFactors = FALSE)
  fields <- c("output_rds", "output_size", "output_mtime", "cells_after_qc", "resolutions", "random_seed")
  if (nrow(provenance) != 1L || !all(fields %in% names(provenance))) stop("Invalid Stage-5 manifest.")
  info <- file.info(input)
  if (!identical(normalizePath(input), normalizePath(provenance$output_rds, mustWork = TRUE)) ||
      !isTRUE(all.equal(as.numeric(info$size), as.numeric(provenance$output_size))) ||
      !isTRUE(all.equal(as.numeric(info$mtime), as.numeric(provenance$output_mtime))) ||
      !identical(seed, as.integer(provenance$random_seed)) ||
      !all(c("0.5", "0.8") %in% strsplit(provenance$resolutions, "|", fixed = TRUE)[[1]])) {
    stop("Stage-5 RDS/provenance mismatch or missing anchor resolutions.")
  }
  message(Sys.time(), ": reading ", input)
  obj <- timed("load_input", readRDS(input))
  ids <- colnames(obj)
  md <- obj[[]]
  anchors <- c("Xenium_snn_res.0.5", "Xenium_snn_res.0.8")
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(rownames(md), ids) || length(ids) != provenance$cells_after_qc ||
      !all(c("orig.ident", "regional_consensus_label_pooled", anchors) %in% names(md))) {
    stop("Missing/invalid metadata, cell IDs/order, count, or anchor partitions.")
  }
  samples_in_object <- as.character(md$orig.ident)
  if (anyNA(samples_in_object) || any(!nzchar(samples_in_object)) ||
      !setequal(samples_in_object, samples)) stop("Object sample coverage differs from authoritative manifest.")
  pooled_labels <- as.character(md$regional_consensus_label_pooled)
  if (anyNA(pooled_labels) || any(!nzchar(pooled_labels))) stop("Missing pooled input labels.")
  if (!"Xenium" %in% SeuratObject::Assays(obj) || !"Xenium_snn" %in% names(obj@graphs) ||
      !"umap_clean" %in% names(obj@reductions)) stop("Missing Xenium assay, saved SNN, or umap_clean.")
  graph <- obj[["Xenium_snn"]]
  if (!identical(rownames(graph), ids) || !identical(colnames(graph), ids) ||
      any(!is.finite(graph@x))) stop("Invalid saved SNN values or cell order.")
  coords <- SeuratObject::Embeddings(obj[["umap_clean"]])
  if (!identical(rownames(coords), ids) || any(!is.finite(coords))) stop("Invalid saved UMAP.")
  commands <- obj@commands[grepl("^FindClusters", names(obj@commands))]
  commands <- Filter(function(x) identical(x@params$graph.name, "Xenium_snn"), commands)
  if (length(commands) != 1L) stop("Cannot identify unique production clustering parameters.")
  cp <- commands[[1]]@params
  required <- c("algorithm", "modularity.fxn", "n.start", "n.iter", "group.singletons", "random.seed")
  if (!all(required %in% names(cp)) || !identical(as.integer(cp$algorithm), 1L) ||
      length(cp$random.seed) != 1L || !is.finite(cp$random.seed)) stop("Unsupported/incomplete saved Louvain parameters.")
  if (!all(as.character(c(8, 12, 17, 19, 20, 21, 24)) %in% as.character(md[[anchors[2]]]))) {
    stop("Expected VZ review targets absent; review the input before clustering.")
  }
  key <- list(input = normalizePath(input), size = as.numeric(info$size), mtime = as.numeric(info$mtime),
    stage5_manifest_md5 = unname(tools::md5sum(input_manifest)),
    config_md5 = unname(tools::md5sum(attr(config, "config_path"))),
    sample_manifest_md5 = unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))),
    implementation_md5 = unname(tools::md5sum(here("scripts", c("xenium_vz_05d_compare_resolutions.R", "R/vz_resolution_efficiency.R")))),
    seurat_version = as.character(utils::packageVersion("Seurat")),
    object_version = as.character(utils::packageVersion("SeuratObject")),
    cell_ids = ids, sample_ids = samples_in_object, clustering_parameters = cp[required])
  if (resume_checkpoint) {
    checkpoint <- readRDS(checkpoint_path)
    vz_comparison_validate_checkpoint(checkpoint, key, columns, md[anchors])
    assignments <- checkpoint$assignments
    for (column in columns) obj[[column]] <- factor(assignments[[column]],
      levels = as.character(sort(unique(as.numeric(assignments[[column]])))))
    message("Reusing validated partitions: 0.6/0.7 clustering skipped.")
    rm(checkpoint)
  } else {
  assignments <- data.frame(cell_id = ids, sample_id = samples_in_object, stringsAsFactors = FALSE)
  for (i in seq_along(resolutions)) {
    if (resolutions[i] %in% c(.5, .8)) {
      values <- as.character(md[[paste0("Xenium_snn_res.", resolutions[i])]])
    } else {
      message(Sys.time(), ": clustering saved SNN at resolution ", resolutions[i])
      partition <- timed(paste0("cluster_res", resolutions[i]), do.call(Seurat::FindClusters, c(list(object = graph, resolution = resolutions[i], verbose = TRUE), cp[required])))
      if (ncol(partition) != 1L || !identical(rownames(partition), ids)) stop("Clustering changed cell IDs/order.")
      values <- as.character(partition[[1]])
    }
    numeric_values <- suppressWarnings(as.numeric(values))
    if (length(values) != length(ids) || anyNA(numeric_values) || any(!is.finite(numeric_values)) ||
        any(numeric_values != floor(numeric_values)) || length(unique(values)) < 2L) stop("Invalid partition.")
    assignments[[columns[i]]] <- values
    obj[[columns[i]]] <- factor(values, levels = as.character(sort(unique(numeric_values))))
  }
    checkpoint <- list(key = key, assignments = assignments)
    vz_comparison_validate_checkpoint(checkpoint, key, columns, md[anchors])
    temporary_checkpoint <- tempfile(pattern = ".partitions_", tmpdir = parent)
    saveRDS(checkpoint, temporary_checkpoint, compress = FALSE)
    if (file.exists(checkpoint_path)) {
      backup <- tempfile(pattern = paste0(run_id, "_previous_partitions_"), tmpdir = parent)
      if (!"--overwrite" %in% args || !file.rename(checkpoint_path, backup)) stop("Cannot archive protected checkpoint.")
    }
    if (!file.rename(temporary_checkpoint, checkpoint_path)) stop("Cannot publish partition checkpoint.")
    rm(checkpoint)
  }
  rm(graph, commands, coords, key)
  for (graph_name in names(obj@graphs)) obj[[graph_name]] <- NULL
  invisible(gc())
  # Joining assay layers is needed for descriptive marker tests, never for clustering.
  obj <- timed("join_layers", Seurat::JoinLayers(obj, assay = "Xenium"))
  obj <- timed("slim_working_object", vz_comparison_slim(obj))
  if (!identical(colnames(obj), ids) || !all(vapply(anchors, function(column)
      identical(as.character(obj[[column, drop = TRUE]]), as.character(md[[column]])), logical(1)))) {
    stop("Layer joining changed IDs/order or saved anchor assignments.")
  }
  staging <- tempfile(pattern = paste0(".staging_", run_id, "_"), tmpdir = parent)
  if (!dir.create(staging)) stop("Cannot create staging directory.")
  published <- FALSE
  on.exit(if (!published) message("Partial diagnostic files retained only at: ", staging), add = TRUE)
  if (!file.copy(checkpoint_path, file.path(staging, "partitions_checkpoint.rds"), overwrite = FALSE)) stop("Cannot copy diagnostic checkpoint.")
  write_table <- function(x, name) write.csv(x, file.path(staging, name), row.names = FALSE)
  connection <- gzfile(file.path(staging, "assignments.csv.gz"), "wt")
  tryCatch(write.csv(assignments, connection, row.names = FALSE), finally = close(connection))
  tree_data <- as.data.frame(lapply(assignments[columns], as.numeric), check.names = FALSE)
  ggplot2::ggsave(file.path(staging, "clustree.pdf"),
    regional_stability_tree(tree_data, prefix = "VZ_compare_res."), width = 12, height = 12)
  transitions <- lapply(combn(seq_along(resolutions), 2, simplify = FALSE), function(pair) {
    x <- regional_resolution_contingency(assignments[[columns[pair[1]]]], assignments[[columns[pair[2]]]], ids)
    x$from_resolution <- resolutions[pair[1]]; x$to_resolution <- resolutions[pair[2]]
    x
  })
  write_table(do.call(rbind, transitions), "all_pair_transitions.csv")
  anchored <- do.call(rbind, lapply(1:3, function(i) {
    x <- regional_resolution_contingency(assignments[[columns[4]]], assignments[[columns[i]]], ids)
    x$from_resolution <- .8; x$to_resolution <- resolutions[i]; x
  }))
  write_table(anchored, "production0.8_transitions.csv")
  targets <- setNames(paste0("Production cluster ", c(8, 12, 17, 19, 20, 21, 24), " (identity/boundary review)"),
    as.character(c(8, 12, 17, 19, 20, 21, 24)))
  if (!all(names(targets) %in% unique(assignments[[columns[4]]]))) stop("Expected VZ review targets absent; review the input.")
  focus <- anchored[anchored$from_cluster %in% names(targets), ]
  focus$review_target <- unname(targets[focus$from_cluster])
  write_table(focus, "focus_population_transitions.csv")
  summaries <- compositions <- label_compositions <- list()
  for (i in seq_along(resolutions)) {
    res <- resolutions[i]; column <- columns[i]; tag <- paste0("res", res, "_")
    values <- assignments[[column]]
    counts <- table(values)
    summaries[[i]] <- data.frame(resolution = res, cluster = names(counts), cells = as.integer(counts),
      samples = vapply(names(counts), function(x) length(unique(samples_in_object[values == x])), integer(1)),
      assignment_source = if (res %in% c(.5, .8)) "saved_production" else "saved_SNN_new_partition")
    composition <- as.data.frame(table(cluster = values, sample_id = factor(samples_in_object, levels = samples)))
    names(composition)[3] <- "cells"; composition$resolution <- res; compositions[[i]] <- composition
    label_composition <- as.data.frame(table(cluster = values, pooled_input_label = pooled_labels))
    names(label_composition)[3] <- "cells"; label_composition$resolution <- res
    label_compositions[[i]] <- label_composition
    p <- Seurat::DimPlot(obj, reduction = "umap_clean", group.by = column, label = TRUE,
      label.size = 4, repel = TRUE, raster = TRUE) + ggplot2::ggtitle(paste("VZ resolution", res))
    ggplot2::ggsave(file.path(staging, paste0(tag, "UMAP.pdf")), p, width = 12, height = 9)
    SeuratObject::Idents(obj) <- column
    message(Sys.time(), ": descriptive markers for resolution ", res)
    markers <- timed(paste0("markers_res", res), Seurat::FindAllMarkers(obj, assay = "Xenium", only.pos = TRUE, min.pct = .25,
      logfc.threshold = .25, test.use = "wilcox", max.cells.per.ident = 1000, random.seed = seed))
    if (!nrow(markers) || !all(c("cluster", "gene", "avg_log2FC") %in% names(markers)) ||
        any(!is.finite(markers$avg_log2FC))) stop("No valid marker table for resolution ", res)
    markers <- markers[order(as.character(markers$cluster), -markers$avg_log2FC, markers$gene), ]
    top <- do.call(rbind, lapply(split(markers, as.character(markers$cluster)), head, n = 50))
    write_table(top, paste0(tag, "top50_markers.csv"))
    marker_counts <- table(factor(as.character(markers$cluster), levels = names(counts)))
    status <- data.frame(cluster = names(counts), cells = as.integer(counts),
      positive_markers = as.integer(marker_counts), status = ifelse(as.integer(counts) < 3L,
        "fewer_than_3_cells_not_testable", ifelse(marker_counts > 0, "markers_available", "no_passing_markers")))
    if (any(status$status == "no_passing_markers")) stop("Unexpected missing markers in a testable cluster; inspect warnings.")
    write_table(status, paste0(tag, "marker_status.csv"))
    timed(paste0("focused_plot_res", res), write_focused_marker_dotplot(obj, marker_sets, assay = "Xenium", group_by = column, branch = "VZ comparison",
      resolution_label = as.character(res), availability_path = file.path(staging, paste0(tag, "marker_availability.csv")),
      tiff_path = if (render_tiff) file.path(staging, paste0(tag, "focused_markers.tif")) else NULL,
      pdf_path = file.path(staging, paste0(tag, "focused_markers.pdf")),
      pdf_device = grDevices::pdf))
    rm(p, markers, top); invisible(gc())
  }
  write_table(do.call(rbind, summaries), "cluster_summary.csv")
  write_table(do.call(rbind, compositions), "sample_composition.csv")
  write_table(do.call(rbind, label_compositions), "input_label_composition.csv")
  parameters <- c(input = normalizePath(input), input_size = info$size, input_mtime = as.numeric(info$mtime),
    cells = length(ids), samples = length(samples), graph = "Xenium_snn (saved; never rebuilt)",
    resolutions = paste(resolutions, collapse = "|"), reused_resolutions = "0.5|0.8", computed_resolutions = "0.6|0.7",
    resume_checkpoint = resume_checkpoint, render_tiff = render_tiff,
    efficiency = "saved SNN and anchors; slim working object; PDF default; validated partition checkpoint",
    optional_presto_available = requireNamespace("presto", quietly = TRUE),
    focused_pdf_device = "grDevices::pdf",
    configured_marker_plot_seed = seed, marker_max_cells_per_ident = 1000, marker_min_pct = .25,
    marker_logfc_threshold = .25, marker_test = "wilcox", marker_only_pos = TRUE,
    production_manifest_md5 = unname(tools::md5sum(input_manifest)),
    samples_md5 = unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))),
    config_md5 = unname(tools::md5sum(attr(config, "config_path"))),
    driver_md5 = unname(tools::md5sum(here("scripts", "xenium_vz_05d_compare_resolutions.R"))),
    run_id = run_id, output_directory = destination, job_id = Sys.getenv("SLURM_JOB_ID", "local"),
    setNames(vapply(cp[required], paste, character(1), collapse = "|"), paste0("clustering_", required)),
    setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), paste0("version_", packages)))
  write_table(data.frame(parameter = names(parameters), value = unname(parameters)), "provenance.csv")
  writeLines(capture.output(sessionInfo()), file.path(staging, "sessionInfo.txt"))
  writeLines(c("Full-data resolution sensitivity only; this is not subsampling stability or automatic annotation.",
    "Saved partitions 0.5 and 0.8 are preserved exactly; 0.6 and 0.7 use the same saved SNN and production clustering seed/settings.",
    "Saved UMAP coordinates are reused for visualization. No normalization, integration, neighbor graph, QC or production labels are changed.",
    "Cell-level Wilcoxon markers (up to 1000 cells per identity) are descriptive, not independent biological-replicate inference.",
    "Marker thresholds match the existing VZ Stage-6 descriptive review: only.pos=TRUE, min.pct=.25, logfc.threshold=.25.",
    "Clusters with fewer than 3 cells are explicitly marked not testable, never removed.",
    "Marker availability tables document panel omissions; gene-wise scaled DotPlot colours cannot establish absolute cross-resolution changes.",
    "Contingency tables include both source recovery and destination purity fractions; cluster numbers have no identity correspondence across resolutions.",
    "All clusters use numeric IDs; current working identities are not reapplied or changed.",
    "The working object is slimmed after clustering: all metadata/cells, joined Xenium counts/data and saved UMAP remain unchanged.",
    "Graphs, unused assays, scaled data and unused reductions are released only in memory; production input is never written.",
    "PDF-only marker plots are the default; --tiff adds the duplicate high-resolution TIFF output.",
    "A diagnostic partition checkpoint permits validated --resume-checkpoint after a downstream failure; the input RDS still loads once.",
    "Inspect VZ, PC/GABA, glial, RL VZ and UBC lineage-marker panels together with exact transitions before choosing a resolution.",
    "No production RDS is written and no outputs become stale from running this diagnostic alone."
  ), file.path(staging, "METHODS.txt"))
  write_table(do.call(rbind, timing_rows), "timings.csv")
  files <- file.path(staging, setdiff(outputs, "outputs.csv"))
  if (!all(file.exists(files)) || any(file.info(files)$size <= 0)) stop("Missing/empty diagnostic outputs.")
  after <- file.info(input)
  if (!identical(info$size, after$size) || !identical(info$mtime, after$mtime)) stop("Input changed during diagnostic.")
  write_table(data.frame(file = basename(files), bytes = file.info(files)$size,
    md5 = unname(tools::md5sum(files))), "outputs.csv")
  if (file.exists(destination)) {
    if (!"--overwrite" %in% args) stop("Destination appeared during run; refusing replacement.")
    backup <- tempfile(pattern = paste0(run_id, "_previous_"), tmpdir = parent)
    if (!file.rename(destination, backup)) stop("Cannot archive previous diagnostic outputs.")
    message("Archived previous diagnostic directory: ", backup)
  }
  if (!file.rename(staging, destination)) stop("Cannot publish complete output set.")
  published <- TRUE
  message("VZ resolution comparison complete: ", destination)
}
main()
