#!/usr/bin/env Rscript
# Export raw-count-derived mean CPM separately for corrected broad and fine labels.
# Run from the project root or scripts/ using the project R/renv environment:
# Rscript --max-ppsize=500000 "scripts/Gene expression matrix.R" --dry-run
# Rscript --max-ppsize=500000 "scripts/Gene expression matrix.R"
# Count layers are aggregated separately to reduce peak memory; input RDS is unchanged.
# Optional project-relative overrides: --input=path --out=directory --assay=Xenium
# Expected HPC runtime: 5-30 minutes (low confidence; depends on RDS size/I/O).
# --dry-run checks paths only (<1 minute); it does not load the merged object.

mean_cpm_by_label <- function(counts, labels, library_size) {
  # Preserve the original exact exclusions; numbered Unknown-XX fine labels
  # remain separate groups. Missing/blank labels are validated before this step.
  keep_label <- !tolower(labels) %in% c("unknown", "unassigned", "unannotated", "na")
  keep <- keep_label & library_size > 0
  if (!any(keep)) stop("No annotated cells with nonzero UMI counts remain.")
  groups <- droplevels(factor(labels[keep]))
  n_cells <- tabulate(as.integer(groups), nbins = nlevels(groups))

  # Fold per-cell CPM scaling and equal-cell averaging into a sparse design.
  # This avoids allocating a second full gene-by-cell CPM matrix.
  design <- Matrix::sparseMatrix(
    i = which(keep), j = as.integer(groups),
    x = 1e6 / library_size[keep] / n_cells[as.integer(groups)],
    dims = c(ncol(counts), nlevels(groups))
  )
  ref <- counts %*% design
  dimnames(ref) <- list(rownames(counts), levels(groups))
  ref <- ref[Matrix::rowSums(ref) > 0, , drop = FALSE]
  if (!nrow(ref)) stop("No expressed genes remain.")
  list(
    matrix = ref,
    cell_counts = data.frame(celltype = levels(groups), n_cells = n_cells),
    excluded_label = sum(!keep_label),
    excluded_zero_umi = sum(keep_label & library_size == 0),
    retained_cells = sum(keep)
  )
}

# Aggregate one count layer at a time; do not join the full gene-by-cell matrix.
mean_cpm_by_layers <- function(assay_obj, meta, label_columns) {
  layer_names <- if (inherits(assay_obj, "Assay5")) {
    SeuratObject::Layers(assay_obj, search = "^counts($|[.])")
  } else "counts"
  if (!length(layer_names)) stop("No raw count layers found.")
  genes <- rownames(assay_obj)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Assay must contain unique, nonblank gene names.")
  }
  states <- setNames(lapply(label_columns, function(label) {
    labels <- trimws(as.character(meta[[label]]))
    keep_label <- !tolower(labels) %in% c("unknown", "unassigned", "unannotated", "na")
    groups <- levels(factor(labels[keep_label]))
    list(sum_cpm = matrix(0, length(genes), length(groups), dimnames = list(genes, groups)),
         n_cells = setNames(integer(length(groups)), groups),
         excluded_label = sum(!keep_label), excluded_zero_umi = 0L)
  }), label_columns)
  seen <- logical(nrow(meta))
  for (i in seq_along(layer_names)) {
    layer <- layer_names[[i]]
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            " | Count layer ", i, "/", length(layer_names), ": ", layer)
    counts <- SeuratObject::LayerData(assay_obj, layer = layer)
    cells <- colnames(counts)
    layer_genes <- rownames(counts)
    cell_index <- match(cells, rownames(meta))
    gene_index <- match(layer_genes, genes)
    if (!nrow(counts) || !ncol(counts) || is.null(cells) || is.null(layer_genes) ||
        anyNA(cell_index) || anyNA(gene_index) || anyDuplicated(cells) ||
        anyDuplicated(layer_genes) || any(seen[cell_index])) {
      stop("Count layer ", layer, " has missing, unmatched, or duplicate gene/cell IDs; layers must not overlap.")
    }
    seen[cell_index] <- TRUE
    count_values <- if (inherits(counts, "sparseMatrix")) counts@x else counts
    if (any(!is.finite(count_values)) || any(count_values < 0)) {
      stop("Count layer ", layer, " contains non-finite or negative values.")
    }
    rm(count_values)
    library_size <- Matrix::colSums(counts)
    if (any(!is.finite(library_size))) stop("Non-finite library sizes in ", layer)
    for (label in label_columns) {
      labels <- trimws(as.character(meta[[label]][cell_index]))
      keep_label <- !tolower(labels) %in% c("unknown", "unassigned", "unannotated", "na")
      state <- states[[label]]
      state$excluded_zero_umi <- state$excluded_zero_umi + sum(keep_label & library_size == 0)
      if (any(keep_label & library_size > 0)) {
        result <- mean_cpm_by_label(counts, labels, library_size)
        ref <- result$matrix
        group_names <- result$cell_counts$celltype
        n_cells <- result$cell_counts$n_cells
        # Weight each layer's mean by its contributing cell count. Only after
        # summing all layers divide by the global count, preserving equal-cell
        # weighting even when sample sizes differ or a layer lacks a group.
        state$sum_cpm[rownames(ref), group_names] <-
          state$sum_cpm[rownames(ref), group_names, drop = FALSE] +
          sweep(as.matrix(ref), 2L, n_cells, "*")
        state$n_cells[group_names] <- state$n_cells[group_names] + n_cells
        rm(result, ref)
      }
      states[[label]] <- state
    }
    rm(counts, library_size, state)
    invisible(gc(verbose = FALSE))
  }
  if (!all(seen)) stop("Raw count layers are missing ", sum(!seen), " metadata cells.")
  setNames(lapply(label_columns, function(label) {
    state <- states[[label]]
    keep <- state$n_cells > 0
    if (!any(keep)) stop("No annotated cells with nonzero UMI counts remain for ", label)
    n_cells <- state$n_cells[keep]
    ref <- sweep(state$sum_cpm[, keep, drop = FALSE], 2L, n_cells, "/")
    ref <- ref[rowSums(ref) > 0, , drop = FALSE]
    if (!nrow(ref)) stop("No expressed genes remain for ", label)
    if (sum(n_cells) + state$excluded_label + state$excluded_zero_umi != nrow(meta)) {
      stop("Retained and excluded cell counts do not cover the object for ", label)
    }
    list(matrix = Matrix::Matrix(ref, sparse = TRUE),
         cell_counts = data.frame(celltype = names(n_cells), n_cells = unname(n_cells)),
         excluded_label = state$excluded_label, excluded_zero_umi = state$excluded_zero_umi,
         retained_cells = sum(n_cells))
  }), label_columns)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c("--dry-run", "--overwrite", "--help")
  options <- c("input", "out", "assay")
  option_names <- sub("=.*$", "", sub("^--", "", args))
  valid <- args %in% flags | (grepl("^--[^=]+=.+$", args) & option_names %in% options)
  if (any(!valid) || anyDuplicated(option_names)) {
    stop("Use each option once: --input=path --out=directory --assay=name --dry-run --overwrite --help")
  }
  if ("--help" %in% args) {
    cat("Export corrected_broad_label and comb_subcluster mean CPM matrices.\n",
        "Options: --input=path --out=directory --assay=Xenium --dry-run --overwrite\n",
        "Relative paths are resolved from the project root.\n", sep = "")
    return(invisible(NULL))
  }
  dry_run <- "--dry-run" %in% args
  overwrite <- "--overwrite" %in% args
  if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
  value <- function(key, default) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit)) substring(hit, nchar(key) + 4L) else default
  }
  if (!requireNamespace("here", quietly = TRUE)) stop("The project 'here' package is required.")
  here::i_am("scripts/Gene expression matrix.R")
  source(here::here("scripts", "R", "config.R"), local = TRUE)
  config <- load_pipeline_config()
  nonspatial_root <- file.path(config$project$outputs_dir, "xenium", "vz_rl", "nonspatial")
  input_path <- normalizePath(here::here(value("input", file.path(
    nonspatial_root, "03_merged", "rds", "Xenium_VZRL_merged.rds"
  ))), mustWork = FALSE)
  out <- normalizePath(here::here(value("out", file.path(
    nonspatial_root, "gene_expression_matrix"
  ))), mustWork = FALSE)
  assay <- value("assay", "Xenium")
  label_columns <- c("corrected_broad_label", "comb_subcluster")
  outputs <- unlist(lapply(label_columns, function(label) file.path(
    out, paste0(label, c("_mean_CPM.tsv", "_mean_CPM.rds", "_cell_counts.csv"))
  )), use.names = FALSE)
  manifest_path <- file.path(out, "gene_expression_export_manifest.csv")
  outputs <- c(outputs, manifest_path)
  if (input_path %in% normalizePath(outputs, mustWork = FALSE)) {
    stop("An output path would replace the input object.")
  }
  message("Input: ", input_path, "\nAssay: ", assay, "\nOutput directory: ", out)
  if (dry_run) {
    ok <- compact_dry_run(
      "Separate corrected broad/subcluster mean CPM exports (paths only)",
      inputs = input_path, outputs = outputs,
      checks = c(outputs_protected = !any(file.exists(outputs)))
    )
    if (!ok) stop("Dry-run failed; review missing inputs or existing outputs.")
    return(invisible(NULL))
  }
  if (!file.exists(input_path)) stop("Merged object not found: ", input_path)
  existing <- outputs[file.exists(outputs)]
  if (length(existing) && !overwrite) {
    stop("Refusing to overwrite exports:\n", paste(existing, collapse = "\n"),
         "\nUse --overwrite only after review.")
  }
  for (package in c("Seurat", "Matrix")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing project dependency: ", package, ". Restore the approved renv environment before running.")
    }
  }
  expected_samples <- as.character(load_sample_manifest(config)$sample_id)
  if (!length(expected_samples) || anyNA(expected_samples) ||
      any(!nzchar(expected_samples)) || anyDuplicated(expected_samples)) {
    stop("Sample manifest must contain unique, nonblank sample IDs.")
  }
  input_info <- file.info(input_path)
  source(here::here("scripts", "R", "downstream_resource_logging.R"), local = TRUE)
  log_downstream_memory(paste("BEFORE RDS LOAD; file GiB=", round(input_info$size / 1024^3, 2)))
  obj <- readRDS(input_path)
  log_downstream_memory("RDS LOADED", obj)
  if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object.")
  meta <- obj[[]]
  required <- c(label_columns, "orig.ident")
  if (!all(required %in% colnames(meta))) {
    stop("Missing metadata: ", paste(setdiff(required, colnames(meta)), collapse = ", "))
  }
  if (anyDuplicated(colnames(obj)) || anyDuplicated(rownames(meta)) ||
      !identical(rownames(meta), colnames(obj))) {
    stop("Object cells and metadata must have unique, identically ordered cell IDs.")
  }
  sample_ids <- as.character(meta$orig.ident)
  if (anyNA(sample_ids) || !setequal(unique(sample_ids), expected_samples)) {
    stop("Object orig.ident must contain exactly the samples in config/samples.csv.")
  }
  for (label in label_columns) {
    labels <- as.character(meta[[label]])
    if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
      stop("Missing or blank labels in ", label, "; resolve before exporting.")
    }
  }
  if (!assay %in% SeuratObject::Assays(obj)) stop("Assay missing: ", assay)
  assay_obj <- obj[[assay]]
  meta <- meta[, required, drop = FALSE]
  input_cells <- ncol(obj)
  rm(obj)
  log_downstream_memory("BEFORE LAYER AGGREGATION", assay_obj)
  results <- mean_cpm_by_layers(assay_obj, meta, label_columns)
  rm(assay_obj)
  log_downstream_memory("AGGREGATION COMPLETE")
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) stop("Cannot create output directory: ", out)
  manifest <- do.call(rbind, lapply(label_columns, function(label) {
    result <- results[[label]]
    prefix <- file.path(out, label)
    write.table(data.frame(GENE = rownames(result$matrix), as.matrix(result$matrix), check.names = FALSE),
                paste0(prefix, "_mean_CPM.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    saveRDS(result$matrix, paste0(prefix, "_mean_CPM.rds"))
    write.csv(result$cell_counts, paste0(prefix, "_cell_counts.csv"), row.names = FALSE)
    message(label, ": ", nrow(result$matrix), " genes, ", ncol(result$matrix),
            " groups, ", result$retained_cells, " retained cells; excluded labels=",
            result$excluded_label, ", annotated zero-UMI cells=", result$excluded_zero_umi)
    data.frame(
      input_path = input_path, input_size = input_info$size,
      input_mtime = as.character(input_info$mtime), assay = assay, label_column = label,
      method = "arithmetic mean of per-cell raw-count CPM; no log transform",
      excluded_labels = "unknown;unassigned;unannotated;na (case-insensitive exact match)",
      numbered_unknown_policy = "retained as separate fine-label groups",
      input_cells = input_cells, samples = length(expected_samples),
      retained_cells = result$retained_cells, excluded_label = result$excluded_label,
      excluded_zero_umi = result$excluded_zero_umi,
      genes = nrow(result$matrix), groups = ncol(result$matrix),
      tsv_path = paste0(prefix, "_mean_CPM.tsv"), rds_path = paste0(prefix, "_mean_CPM.rds"),
      cell_counts_path = paste0(prefix, "_cell_counts.csv"),
      exported_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  }))
  write.csv(manifest, manifest_path, row.names = FALSE)
  message("Export completed: ", manifest_path)
  invisible(manifest)
}

if (sys.nframe() == 0L) main()
