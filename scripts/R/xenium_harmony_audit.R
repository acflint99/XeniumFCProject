# Quantitative, read-only pre-/post-Harmony neighborhood audit.

parse_pcw_numeric <- function(values) {
  values <- trimws(as.character(values))
  # Production metadata can store either numeric ages or labels such as PCW12.
  values <- sub("^PCW[[:space:]]*", "", values, ignore.case = TRUE)
  suppressWarnings(as.numeric(sub("[^0-9.].*$", "", values)))
}

stratified_audit_cells <- function(metadata, max_cells, seed) {
  n <- nrow(metadata)
  if (n <= max_cells) return(seq_len(n))
  strata <- interaction(
    metadata$orig.ident, metadata$audit_biology_label,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(n), strata)
  allocation <- pmax(1L, floor(as.double(max_cells) * lengths(groups) / n))
  while (sum(allocation) > max_cells) {
    eligible <- which(allocation > 1L)
    allocation[eligible[[which.max(allocation[eligible])]]] <-
      allocation[eligible[[which.max(allocation[eligible])]]] - 1L
  }
  while (sum(allocation) < max_cells) {
    capacity <- lengths(groups) - allocation
    eligible <- which(capacity > 0L)
    if (!length(eligible)) break
    allocation[eligible[[which.max(capacity[eligible])]]] <-
      allocation[eligible[[which.max(capacity[eligible])]]] + 1L
  }
  set.seed(seed)
  sort(unlist(Map(function(index, size) {
    sample(index, size = min(size, length(index)), replace = FALSE)
  }, groups, allocation), use.names = FALSE))
}

neighbor_metrics <- function(embedding, metadata, k) {
  if (!is.matrix(embedding) || !is.numeric(embedding) ||
      nrow(embedding) != nrow(metadata) || any(!is.finite(embedding))) {
    stop("Harmony-audit neighbor input must be a finite numeric matrix aligned with metadata.")
  }
  if (length(k) != 1L || !is.finite(k) || k < 1L || k != as.integer(k) ||
      nrow(embedding) <= k) {
    stop("Harmony-audit neighbor count must be a positive integer below the audited cell count.")
  }
  nn <- RANN::nn2(embedding, query = embedding, k = k + 1L)$nn.idx[, -1L, drop = FALSE]
  samples <- as.character(metadata$orig.ident)
  biology <- as.character(metadata$audit_biology_label)
  pcw <- metadata$audit_pcw_numeric
  sample_neighbors <- matrix(samples[nn], nrow = nrow(nn))
  biology_neighbors <- matrix(biology[nn], nrow = nrow(nn))
  pcw_neighbors <- matrix(pcw[nn], nrow = nrow(nn))
  result <- data.frame(
    cell_id = rownames(metadata),
    sample_id = samples,
    biology_label = biology,
    different_sample_neighbor_fraction = rowMeans(sample_neighbors != samples),
    same_biology_neighbor_fraction = rowMeans(biology_neighbors == biology),
    mean_absolute_pcw_neighbor_difference = rowMeans(abs(pcw_neighbors - pcw), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(metrics = result, neighbor_index = nn)
}

summarize_harmony_metric <- function(metrics, stage, group_column = NULL) {
  if (is.null(group_column)) {
    groups <- list(all_cells = seq_len(nrow(metrics)))
  } else {
    groups <- split(seq_len(nrow(metrics)), metrics[[group_column]])
  }
  do.call(rbind, lapply(names(groups), function(group) {
    index <- groups[[group]]
    data.frame(
      stage = stage,
      group_type = if (is.null(group_column)) "all" else group_column,
      group = group,
      n_cells = length(index),
      median_different_sample_neighbor_fraction = stats::median(
        metrics$different_sample_neighbor_fraction[index], na.rm = TRUE
      ),
      median_same_biology_neighbor_fraction = stats::median(
        metrics$same_biology_neighbor_fraction[index], na.rm = TRUE
      ),
      median_absolute_pcw_neighbor_difference = stats::median(
        metrics$mean_absolute_pcw_neighbor_difference[index], na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
}

run_harmony_audit <- function(region,
                              input_path,
                              output_dir,
                              overwrite = FALSE,
                              dimensions = 30L,
                              neighbors = 30L,
                              max_cells = 100000L,
                              seed = 42L) {
  region <- toupper(region)
  output_paths <- c(
    summary = file.path(output_dir, paste0("Xenium_", region, "_Harmony_Audit_Summary.csv")),
    paired = file.path(output_dir, paste0("Xenium_", region, "_Harmony_Audit_Paired_Metrics.csv")),
    provenance = file.path(output_dir, paste0("Xenium_", region, "_Harmony_Audit_Provenance.csv"))
  )
  existing <- output_paths[file.exists(output_paths)]
  if (length(existing) && !overwrite) {
    stop(
      "Refusing to overwrite Harmony-audit output(s):\n- ",
      paste(existing, collapse = "\n- "),
      "\nRerun with --overwrite only after reviewing the existing files."
    )
  }
  if (!file.exists(input_path)) stop("Harmony-audit input not found: ", input_path)
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("The project renv must provide RANN for the Harmony neighborhood audit.")
  }
  log_stage <- function(text) {
    message(format(Sys.time(), tz = "UTC", usetz = TRUE), " | ", region, " | ", text)
  }
  log_stage(paste0("Loading integrated RDS: ", input_path))
  object <- readRDS(input_path)
  total_cells <- ncol(object)
  log_stage(paste0("RDS loaded: ", total_cells, " cells; extracting audit inputs."))
  if (!all(c("pca", "harmony") %in% Seurat::Reductions(object))) {
    stop("Harmony-audit input must contain pca and harmony reductions.")
  }
  required_metadata <- c("orig.ident", "PCW", "regional_consensus_label_pooled")
  missing <- setdiff(required_metadata, names(object[[]]))
  if (length(missing)) stop("Harmony-audit input lacks: ", paste(missing, collapse = ", "))
  # Retain only the three source fields used by the audit, rather than copying
  # every integrated-object metadata column into the sampled data frame.
  metadata <- object[[]][, required_metadata, drop = FALSE]
  metadata$audit_biology_label <- trimws(as.character(
    metadata$regional_consensus_label_pooled
  ))
  metadata$audit_pcw_numeric <- parse_pcw_numeric(metadata$PCW)
  if (anyNA(metadata$audit_biology_label) || any(!nzchar(metadata$audit_biology_label)) ||
      anyNA(metadata$audit_pcw_numeric)) {
    stop("Harmony-audit biology labels or PCW values are invalid.")
  }
  selected <- stratified_audit_cells(metadata, max_cells, seed)
  metadata <- metadata[selected, , drop = FALSE]
  if (ncol(Seurat::Embeddings(object, "pca")) < dimensions ||
      ncol(Seurat::Embeddings(object, "harmony")) < dimensions) {
    stop("Harmony-audit input lacks the requested number of embedding dimensions.")
  }
  pre <- Seurat::Embeddings(object, "pca")[selected, seq_len(dimensions), drop = FALSE]
  post <- Seurat::Embeddings(object, "harmony")[selected, seq_len(dimensions), drop = FALSE]
  if (!identical(rownames(pre), rownames(metadata)) ||
      !identical(rownames(post), rownames(metadata))) {
    stop("Harmony-audit embeddings do not align with selected metadata rows.")
  }
  # Counts, scaled expression and SNN graphs are not used below. Release them
  # before native neighbor searches; the initial read still needs ample memory.
  rm(object)
  invisible(gc())
  log_stage(paste0("Released integrated object; auditing ", nrow(metadata),
                   " cells, ", dimensions, " dimensions and ", neighbors, " neighbors."))
  log_stage("Starting pre-Harmony neighbor search.")
  pre_result <- neighbor_metrics(pre, metadata, neighbors)
  rm(pre)
  invisible(gc())
  log_stage("Starting post-Harmony neighbor search.")
  post_result <- neighbor_metrics(post, metadata, neighbors)
  rm(post)
  invisible(gc())
  log_stage("Neighbor searches complete; compiling paired metrics and summaries.")
  paired <- pre_result$metrics
  metric_columns <- c(
    "different_sample_neighbor_fraction", "same_biology_neighbor_fraction",
    "mean_absolute_pcw_neighbor_difference"
  )
  names(paired)[match(metric_columns, names(paired))] <- paste0("pre_", metric_columns)
  for (column in metric_columns) {
    paired[[paste0("post_", column)]] <- post_result$metrics[[column]]
    paired[[paste0("change_", column)]] <-
      paired[[paste0("post_", column)]] - paired[[paste0("pre_", column)]]
  }
  paired$neighbor_set_overlap_fraction <- vapply(seq_len(nrow(paired)), function(i) {
    length(intersect(
      pre_result$neighbor_index[i, ], post_result$neighbor_index[i, ]
    )) / neighbors
  }, numeric(1))

  summaries <- rbind(
    summarize_harmony_metric(pre_result$metrics, "pre_harmony"),
    summarize_harmony_metric(post_result$metrics, "post_harmony"),
    summarize_harmony_metric(pre_result$metrics, "pre_harmony", "sample_id"),
    summarize_harmony_metric(post_result$metrics, "post_harmony", "sample_id"),
    summarize_harmony_metric(pre_result$metrics, "pre_harmony", "biology_label"),
    summarize_harmony_metric(post_result$metrics, "post_harmony", "biology_label")
  )
  input_info <- file.info(input_path)
  provenance <- data.frame(
    region = region,
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_rds = normalizePath(input_path, mustWork = TRUE),
    input_size = as.numeric(input_info$size),
    input_mtime = as.numeric(input_info$mtime),
    total_cells = total_cells,
    audited_cells = nrow(metadata),
    dimensions = dimensions,
    neighbors = neighbors,
    maximum_audit_cells = max_cells,
    random_seed = seed,
    interpretation = paste(
      "higher sample mixing is desirable only when biology-label purity and",
      "developmental-stage neighborhood structure are retained"
    ),
    integration_or_clustering_modified = FALSE,
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summaries, output_paths[["summary"]], row.names = FALSE)
  write.csv(paired, output_paths[["paired"]], row.names = FALSE)
  write.csv(provenance, output_paths[["provenance"]], row.names = FALSE)
  log_stage("Saved all three Harmony-audit tables; integrated input was not modified.")
  invisible(output_paths)
}
