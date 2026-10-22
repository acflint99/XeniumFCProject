#!/usr/bin/env Rscript
# Read-only Stage-5 review: clustree and sample-stratified cell subsampling.
# Fixed production Harmony representation; no expression assays are reclustered.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_regional_05c_resolution_stability.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "regional_cluster_stability.R"))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  flags <- c("--dry-run", "--list", "--overwrite", "--check-packages")
  if (any(!args %in% flags & !grepl("^--(region|repeats|fraction|run-id)=.+$", args))) {
    stop("Usage: --region=VZ|RL [--repeats=5] [--fraction=0.8] [--run-id=pilot5] ",
         "[--list|--dry-run|--check-packages] [--overwrite]")
  }
  option <- function(key, default) {
    values <- sub(paste0("^--", key, "="), "", args[startsWith(args, paste0("--", key, "="))])
    if (length(values) > 1L) stop("Repeated option: ", key)
    if (length(values)) values else default
  }
  region <- option("region", "")
  if (!region %in% c("VZ", "RL")) stop("Specify --region=VZ or --region=RL.")
  repeats <- suppressWarnings(as.numeric(option("repeats", "5")))
  fraction <- suppressWarnings(as.numeric(option("fraction", "0.8")))
  if (length(repeats) != 1L || !is.finite(repeats) || repeats < 2 ||
      repeats > 1000 || repeats != floor(repeats)) stop("repeats must be an integer from 2 to 1000.")
  if (!is.finite(fraction) || fraction <= 0 || fraction >= 1) stop("fraction must be between 0 and 1.")
  run_id <- option("run-id", paste0("pilot", repeats))
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) stop("Invalid run-id.")
  inspect_only <- any(c("--dry-run", "--list", "--check-packages") %in% args)
  if (inspect_only && "--overwrite" %in% args) stop("Inspection and overwrite cannot be combined.")
  config <- load_pipeline_config()
  samples <- load_sample_manifest(config)$sample_id
  if (!length(samples) || anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("Invalid authoritative sample manifest.")
  }
  root <- here(config$project$outputs_dir, "xenium", tolower(region), "05_post_qc")
  input <- file.path(root, "rds", paste0("Xenium_", region, "_postQC.rds"))
  provenance <- file.path(root, "tables", paste0("Xenium_", region, "_postQC_manifest.csv"))
  parent <- file.path(root, "diagnostics", "resolution_stability")
  destination <- file.path(parent, run_id)
  output_names <- c("clustree.pdf", "stability.pdf", "cluster_summary.csv", "recovery_by_repeat.csv",
                    "sample_composition.csv", "sampling_by_repeat.csv", "resolution_summary.csv",
                    "provenance.csv", "samples.csv", "sessionInfo.txt", "METHODS.txt", "outputs.csv")
  packages <- c("Seurat", "SeuratObject", "Matrix", "clustree", "ggraph", "ggplot2", "future")
  if ("--check-packages" %in% args) {
    available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    print(data.frame(package = packages, available = available), row.names = FALSE)
    if (!all(available)) stop("Missing dependencies. Use an explicit project renv installation before submitting.")
    regional_stability_check_plot()
    return(invisible(NULL))
  }
  if ("--list" %in% args) {
    cat("Region:", region, "| manifest samples:", length(samples), "| repeats:", repeats,
        "| fraction:", fraction, "\nInput:", input, "\nOutput:", destination, "\n")
    return(invisible(NULL))
  }
  if ("--dry-run" %in% args) {
    ok <- compact_dry_run(paste(region, "resolution stability (paths only; no RDS load)"),
      inputs = c(input, provenance), outputs = file.path(destination, output_names),
      checks = c(manifest_valid = TRUE, destination_protected =
        !file.exists(destination)))
    if (!ok) stop("Dry-run failed; inspect missing inputs or select a fresh run-id.")
    return(invisible(NULL))
  }
  if (!all(file.exists(c(input, provenance)))) stop("Missing Stage-5 RDS or provenance: ", root)
  if (file.exists(destination) && !"--overwrite" %in% args) {
    stop("Output directory already exists. Use a new --run-id, or --overwrite after review: ", destination)
  }
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing dependencies: ", paste(missing, collapse = ", "),
                            ". Install explicitly through the project renv environment.")
  # Exercise the exact PDF rendering path before loading a large production RDS.
  regional_stability_check_plot()
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(parent, paste0(".lock_", run_id))
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Another run owns this run-id, or a failed job left its lock: ", lock,
         ". Check the corresponding job before removing an abandoned lock.")
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  future::plan(future::sequential)
  base_seed <- as.integer(config$runtime$random_seed)
  if (length(base_seed) != 1L || is.na(base_seed) || base_seed < 0 ||
      base_seed > .Machine$integer.max - repeats) stop("Invalid configured seed.")
  manifest <- read.csv(provenance, stringsAsFactors = FALSE)
  fields <- c("output_rds", "output_size", "output_mtime", "cells_after_qc", "resolutions", "random_seed")
  if (nrow(manifest) != 1L || !all(fields %in% names(manifest))) stop("Invalid Stage-5 provenance.")
  info <- file.info(input)
  if (!identical(normalizePath(input), normalizePath(manifest$output_rds, mustWork = TRUE)) ||
      !isTRUE(all.equal(as.numeric(info$size), as.numeric(manifest$output_size))) ||
      !isTRUE(all.equal(as.numeric(info$mtime), as.numeric(manifest$output_mtime))) ||
      !identical(base_seed, as.integer(manifest$random_seed))) stop("Stage-5 input does not match its manifest.")
  resolutions <- suppressWarnings(as.numeric(strsplit(manifest$resolutions, "|", fixed = TRUE)[[1]]))
  allowed <- if (region == "VZ") c(0.3, 0.5, 0.8, 1) else c(0.3, 0.5, 0.8)
  if (length(resolutions) < 2L || anyNA(resolutions) || anyDuplicated(resolutions) ||
      !all(resolutions %in% allowed) || !0.8 %in% resolutions) stop("Unexpected Stage-5 resolutions.")
  resolutions <- sort(resolutions)
  columns <- paste0("Xenium_snn_res.", resolutions)
  message(Sys.time(), ": reading ", input)
  obj <- readRDS(input)
  metadata <- obj[[]]
  ids <- colnames(obj)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(rownames(metadata), ids) || ncol(obj) != manifest$cells_after_qc ||
      !all(c("orig.ident", columns) %in% names(metadata))) stop("Invalid cell IDs, count, or required metadata.")
  sample_id <- as.character(metadata$orig.ident)
  if (anyNA(sample_id) || any(!nzchar(sample_id)) || !setequal(unique(sample_id), samples)) {
    stop("Object sample coverage differs from config/samples.csv. Missing: ",
         paste(setdiff(samples, sample_id), collapse = ", "), "; unexpected: ",
         paste(setdiff(sample_id, samples), collapse = ", "))
  }
  labels <- lapply(metadata[columns], as.character)
  if (any(vapply(labels, function(x) anyNA(x) || any(!nzchar(x)) || length(unique(x)) < 2L, logical(1)))) {
    stop("Missing, blank, or single-cluster partitions.")
  }
  reduction <- if ("harmony_clean" %in% names(obj@reductions)) "harmony_clean" else "harmony"
  if (!reduction %in% names(obj@reductions)) stop("Missing production Harmony representation.")
  # Read saved command parameters, including k=20 for older zero-removal reuse.
  # Never silently substitute the configured full-reintegration k=30.
  commands <- obj@commands
  neighbor_commands <- commands[grepl("^FindNeighbors", names(commands))]
  neighbor_commands <- Filter(function(x) identical(x@params$reduction, reduction) &&
    "Xenium_snn" %in% x@params$graph.name, neighbor_commands)
  cluster_commands <- commands[grepl("^FindClusters", names(commands))]
  cluster_commands <- Filter(function(x) identical(x@params$graph.name, "Xenium_snn"), cluster_commands)
  if (length(neighbor_commands) != 1L || length(cluster_commands) != 1L) {
    stop("Cannot identify unique production FindNeighbors/FindClusters command records; review provenance before proceeding.")
  }
  np <- neighbor_commands[[1]]@params
  cp <- cluster_commands[[1]]@params
  required_nn <- c("dims", "k.param", "nn.method", "n.trees", "annoy.metric", "nn.eps", "prune.SNN", "l2.norm")
  required_cluster <- c("algorithm", "modularity.fxn", "n.start", "n.iter", "group.singletons")
  if (!all(required_nn %in% names(np)) || !all(required_cluster %in% names(cp)) ||
      !identical(as.integer(cp$algorithm), 1L)) stop("Incomplete or unsupported production graph/clustering parameters.")
  embedding <- SeuratObject::Embeddings(obj[[reduction]])
  dims <- np$dims
  if (!is.numeric(dims) || !length(dims) || anyNA(dims) || anyDuplicated(dims) ||
      any(dims != floor(dims)) || any(dims < 1 | dims > ncol(embedding)) ||
      !identical(rownames(embedding), ids)) stop("Invalid saved graph dimensions or embedding cell order.")
  embedding <- embedding[, dims, drop = FALSE]
  if (any(!is.finite(embedding))) stop("Non-finite Harmony coordinates.")
  if (length(np$k.param) != 1L || !is.finite(np$k.param) || np$k.param < 1) stop("Invalid saved k.param.")
  rm(obj, metadata, commands, neighbor_commands, cluster_commands)
  invisible(gc())
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(pattern = paste0(".staging_", run_id, "_"), tmpdir = parent)
  if (!dir.create(staging)) stop("Cannot create staging directory.")
  published <- FALSE
  on.exit(if (!published) message("Incomplete results retained only at: ", staging), add = TRUE)
  write_table <- function(x, name) write.csv(x, file.path(staging, name), row.names = FALSE)
  write_table(data.frame(sample_id = samples), "samples.csv")
  tree_data <- as.data.frame(lapply(labels, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE)
  if (any(!vapply(tree_data, function(x) all(is.finite(x)) && all(x == floor(x)), logical(1)))) {
    stop("clustree requires numeric production cluster IDs.")
  }
  tree <- regional_stability_tree(tree_data)
  ggplot2::ggsave(file.path(staging, "clustree.pdf"), tree, width = 12, height = 10, limitsize = FALSE)
  rm(tree_data, tree)
  composition <- do.call(rbind, lapply(seq_along(columns), function(i) {
    x <- as.data.frame(table(cluster = labels[[i]], sample_id = factor(sample_id, levels = samples)), stringsAsFactors = FALSE)
    names(x)[3] <- "cells"
    x$resolution <- resolutions[i]
    x
  }))
  write_table(composition, "sample_composition.csv")
  recovery <- list()
  sampling <- list()
  for (repeat_id in seq_len(repeats)) {
    started <- proc.time()[[3]]
    seed <- base_seed + repeat_id
    retained <- regional_stability_sample(sample_id, fraction, seed)
    if (length(retained) <= np$k.param || length(retained) == length(ids)) stop("Subset is too small for k, or removed no cells.")
    sampling[[repeat_id]] <- data.frame(repeat_id = repeat_id, seed = seed, sample_id = samples,
      original_cells = as.integer(table(factor(sample_id, levels = samples))),
      retained_cells = as.integer(table(factor(sample_id[retained], levels = samples))))
    set.seed(seed)
    message(Sys.time(), ": repeat ", repeat_id, "/", repeats, "; cells=", length(retained), "; k=", np$k.param)
    graphs <- do.call(Seurat::FindNeighbors, c(list(object = embedding[retained, , drop = FALSE],
      compute.SNN = TRUE, return.neighbor = FALSE, verbose = TRUE), np[setdiff(required_nn, "dims")]))
    graph <- graphs$snn
    if (is.null(graph) || !identical(rownames(graph), ids[retained]) ||
        !identical(colnames(graph), ids[retained])) stop("Subsample graph cell order changed.")
    rm(graphs)
    for (i in seq_along(resolutions)) {
      # One partition at a time prevents parallel graph/object duplication.
      clustered <- do.call(Seurat::FindClusters, c(list(object = graph, resolution = resolutions[i],
        random.seed = seed, verbose = TRUE), cp[required_cluster]))
      if (ncol(clustered) != 1L || !identical(rownames(clustered), ids[retained])) stop("Reclustered IDs/order changed.")
      result <- regional_stability_jaccard(labels[[i]][retained], as.character(clustered[[1]]), sort(unique(labels[[i]])))
      result$repeat_id <- repeat_id
      result$seed <- seed
      result$resolution <- resolutions[i]
      recovery[[length(recovery) + 1L]] <- result
      rm(clustered)
    }
    rm(graph)
    invisible(gc())
    sampling[[repeat_id]]$repeat_elapsed_seconds <- proc.time()[[3]] - started
    # Incremental evidence remains in staging if a job fails or times out.
    write_table(do.call(rbind, recovery), "recovery_by_repeat.csv")
    write_table(do.call(rbind, sampling), "sampling_by_repeat.csv")
  }
  recovery <- do.call(rbind, recovery)
  summary <- do.call(rbind, lapply(seq_along(resolutions), function(i) {
    do.call(rbind, lapply(sort(unique(labels[[i]])), function(cluster) {
      rows <- recovery[recovery$resolution == resolutions[i] & recovery$cluster == cluster, ]
      scores <- rows$jaccard[!is.na(rows$jaccard)]
      data.frame(resolution = resolutions[i], cluster = cluster,
        full_cells = sum(labels[[i]] == cluster), repeats = repeats,
        recovered_repeats = length(scores), absent_repeats = sum(is.na(rows$jaccard)),
        mean_jaccard = if (length(scores)) mean(scores) else NA_real_,
        median_jaccard = if (length(scores)) median(scores) else NA_real_,
        sd_jaccard = if (length(scores) > 1L) sd(scores) else NA_real_,
        q10_jaccard = if (length(scores)) unname(quantile(scores, .1)) else NA_real_,
        q90_jaccard = if (length(scores)) unname(quantile(scores, .9)) else NA_real_)
    }))
  }))
  write_table(summary, "cluster_summary.csv")
  resolution_summary <- do.call(rbind, lapply(resolutions, function(res) {
    x <- summary[summary$resolution == res, ]
    measured <- is.finite(x$mean_jaccard)
    data.frame(resolution = res, clusters = nrow(x), clusters_measured = sum(measured),
      mean_cluster_jaccard = if (any(measured)) mean(x$mean_jaccard[measured]) else NA_real_,
      minimum_cluster_jaccard = if (any(measured)) min(x$mean_jaccard[measured]) else NA_real_,
      cell_weighted_jaccard = if (any(measured)) weighted.mean(x$mean_jaccard[measured], x$full_cells[measured]) else NA_real_,
      absent_cluster_repeat_pairs = sum(x$absent_repeats))
  }))
  write_table(resolution_summary, "resolution_summary.csv")
  plot <- ggplot2::ggplot(recovery, ggplot2::aes(x = cluster, y = jaccard)) +
    ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .12, height = 0, seed = base_seed),
                        alpha = .4, size = .7, na.rm = TRUE) +
    ggplot2::facet_wrap(~resolution, scales = "free_x", ncol = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) + ggplot2::theme_bw() +
    ggplot2::labs(title = paste(region, "cluster recovery on subsampled cells"),
      subtitle = "Fixed production Harmony; boxplots show repeat variability, not confidence intervals",
      x = "Production cluster", y = "Best-match Jaccard")
  ggplot2::ggsave(file.path(staging, "stability.pdf"), plot, width = 12, height = 3 * length(resolutions))
  parameters <- c(region = region, input = normalizePath(input), input_size = info$size,
    input_mtime = as.numeric(info$mtime), stage5_manifest = normalizePath(provenance),
    stage5_manifest_md5 = unname(tools::md5sum(provenance)),
    sample_manifest_md5 = unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))),
    config_md5 = unname(tools::md5sum(attr(config, "config_path"))),
    driver_md5 = unname(tools::md5sum(here("scripts", "xenium_regional_05c_resolution_stability.R"))),
    helper_md5 = unname(tools::md5sum(here("scripts", "R", "regional_cluster_stability.R"))),
    cells = length(ids), samples = length(samples), reduction = reduction,
    dims = paste(dims, collapse = "|"), resolutions = paste(resolutions, collapse = "|"),
    repeats = repeats, fraction = fraction, base_seed = base_seed,
    run_id = run_id, job_id = Sys.getenv("SLURM_JOB_ID", "local"),
    production_resolution = "0.8 (unchanged)", output_directory = destination)
  parameters <- c(parameters, setNames(vapply(np[setdiff(required_nn, "dims")], paste, character(1), collapse = "|"),
                                      paste0("neighbors_", setdiff(required_nn, "dims"))),
    setNames(vapply(cp[required_cluster], paste, character(1), collapse = "|"), paste0("clustering_", required_cluster)),
    setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), paste0("version_", packages)))
  write_table(data.frame(parameter = names(parameters), value = unname(parameters)), "provenance.csv")
  writeLines(capture.output(sessionInfo()), file.path(staging, "sessionInfo.txt"))
  writeLines(c(
    "Diagnostic review only: no selected resolution, production graph, label, or RDS was changed.",
    "clustree: Zappia & Oshlack (2018), doi:10.1093/gigascience/giy083.",
    "Cluster-wise Jaccard: Hennig (2007), doi:10.1016/j.csda.2006.11.025.",
    "This is not chooseR and does not implement chooseR's optimization rule.",
    "Each repeat samples without replacement within each biological sample at the same fraction (floor, minimum one).",
    "The complete production Harmony embedding is held fixed; neighbors are rebuilt on retained cells.",
    "Production graph/clustering parameters come from saved Seurat commands. Repeat seeds are base_seed + repeat_id.",
    "For each saved resolution, restrict original labels to retained cells and find the best Jaccard match per cluster.",
    "Matching is many-to-one: merges and splits reduce recovery, but this is not a one-to-one label assignment.",
    "Absent clusters have NA recovery and are counted explicitly; quantiles describe repeat variation, not biological uncertainty.",
    "No universal stability cutoff or automatic resolution ranking is applied. Inspect rare clusters separately.",
    "Cell subsampling tests robustness conditional on integration, not reproducibility across biological specimens.",
    "Compare clusters across resolutions using clustree; cluster integer labels do not correspond across resolutions.",
    "Review marker, sample/age, QC and spatial evidence before approving any new production partition.",
    "Five repeats are a runtime pilot; 20 repeats are exploratory and need convergence assessment before final reporting."
  ), file.path(staging, "METHODS.txt"))
  if (!all(file.exists(file.path(staging, setdiff(output_names, "outputs.csv"))))) stop("Partial output set.")
  # Confirm the input was not replaced while the diagnostic ran.
  after <- file.info(input)
  if (!identical(after$size, info$size) || !identical(after$mtime, info$mtime)) stop("Stage-5 input changed during review.")
  files <- file.path(staging, setdiff(output_names, "outputs.csv"))
  write_table(data.frame(file = basename(files), bytes = file.info(files)$size,
                         md5 = unname(tools::md5sum(files))), "outputs.csv")
  # Directory rename publishes a complete set. Explicit overwrite archives the old set.
  if (file.exists(destination)) {
    if (!"--overwrite" %in% args) stop("Output appeared during run; refusing replacement.")
    backup <- tempfile(pattern = paste0(run_id, "_previous_"), tmpdir = parent)
    if (!file.rename(destination, backup)) stop("Cannot archive previous diagnostic outputs.")
    message("Previous diagnostic outputs archived at: ", backup)
  }
  if (!file.rename(staging, destination)) stop("Cannot publish diagnostic directory: ", staging)
  published <- TRUE
  message("Complete ", region, " resolution review: ", destination)
}
main()
