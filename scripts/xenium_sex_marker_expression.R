#!/usr/bin/env Rscript
# Read-only expression evidence from whole-tissue QC objects; no RDS is written.
# HPC estimate: 1-5 min/sample; 15-70 min for 41 tasks at concurrency 3,
# plus 1-3 min to combine (low confidence; dominated by RDS I/O).
classify_sex_marker_evidence <- function(result, core_y_genes, min_cells) {
  core_y <- result$gene %in% core_y_genes
  y_measured <- sum(result$measured[core_y])
  y_supported <- result$gene[core_y & !is.na(result$supported) & result$supported]
  xist_measured <- isTRUE(result$measured[result$gene == "XIST"])
  xist_supported <- isTRUE(result$supported[result$gene == "XIST"])
  sry_supported <- isTRUE(result$supported[result$gene == "SRY"])
  call <- "unresolved"
  pattern <- "no_supported_markers"
  reason <- "Neither marker pattern meets the exploratory support criteria."
  if (unique(result$n_cells) < min_cells) {
    pattern <- "too_few_cells"
    reason <- "Too few cells for the exploratory sample call."
  } else if ((length(y_supported) > 0L || sry_supported) && xist_supported) {
    pattern <- "conflicting_XIST_and_Y"
    reason <- "Both XIST and Y-marker support: review possible mixture or discordance."
  } else if (length(y_supported) > 0L) {
    call <- "male-like"
    pattern <- if (length(y_supported) >= 2L) "multiple_Y_markers" else "single_Y_marker"
    reason <- paste0("Core Y-marker support: ", paste(y_supported, collapse = ", "), ". ",
      if (!xist_measured) "XIST unavailable; female-associated expression cannot be assessed. " else "XIST not supported. ",
      if (length(y_supported) == 1L) "Single-marker evidence; seek corroboration." else "Multiple Y markers corroborate the pattern.")
  } else if (sry_supported) {
    pattern <- "SRY_only"
    reason <- "Only SRY is supported; seek corroboration from core Y markers."
  } else if (xist_supported && y_measured >= 1L) {
    call <- "female-like"
    pattern <- "XIST_supported_Y_low"
    reason <- "XIST supported without support from measured Y markers."
  } else if (xist_supported) {
    pattern <- "XIST_without_Y_coverage"
    reason <- "XIST supported, but no core Y markers are measured to assess discordance."
  } else if (y_measured == 0L) {
    pattern <- "no_core_Y_coverage"
    reason <- "No core Y markers are measured; insufficient marker coverage."
  } else if (!xist_measured) {
    pattern <- "low_Y_XIST_unavailable"
    reason <- "Measured Y markers lack support and XIST is unavailable; compatible with female origin but insufficient for a female-like call."
  }
  list(provisional_sex = call, evidence_pattern = pattern, reason = reason,
    n_y_markers_measured = y_measured, n_y_markers_supported = length(y_supported),
    supported_y_markers = if (length(y_supported)) paste(y_supported, collapse = ";") else "none",
    missing_markers = if (any(!result$measured)) paste(result$gene[!result$measured], collapse = ";") else "none")
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat(paste0(
      "Usage: Rscript --max-ppsize=500000 scripts/xenium_sex_marker_expression.R\n",
      "  --list | --dry-run | [--dry-run] --task N | [--dry-run] --combine\n",
      "Options: --min-positive-cells 10 --min-positive-pct 1 --min-cells 100\n",
      "         --output-dir PATH (default: outputs/xenium/sex_marker_expression)\n",
      "Calls are exploratory expression patterns, not confirmed sex. No overwrite.\n"
    ))
    return(invisible(NULL))
  }
  flags <- c("--list", "--dry-run", "--combine")
  values <- c("--task", "--min-positive-cells", "--min-positive-pct",
              "--min-cells", "--output-dir")
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c(flags, values) || !is.null(opts[[key]])) {
      stop("Unknown or repeated option: ", key)
    }
    if (key %in% flags) opts[[key]] <- TRUE else {
      i <- i + 1L
      if (i > length(args) || startsWith(args[[i]], "--")) stop("Missing value: ", key)
      opts[[key]] <- args[[i]]
    }
    i <- i + 1L
  }
  number <- function(key, default, integer = FALSE, max_value = Inf) {
    x <- if (is.null(opts[[key]])) default else suppressWarnings(as.numeric(opts[[key]]))
    if (length(x) != 1L || !is.finite(x) || x <= 0 || x > max_value ||
        (integer && x != floor(x))) stop("Invalid positive value for ", key)
    x
  }
  min_positive_cells <- number("--min-positive-cells", 10, TRUE)
  min_positive_pct <- number("--min-positive-pct", 1, max_value = 100)
  min_cells <- number("--min-cells", 100, TRUE)
  for (p in c("here", "yaml")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing dependency: ", p,
      ". Use the existing project renv environment; no packages are installed here.")
  }
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
  source(file.path(script_dir, "R", "config.R"), local = TRUE)
  config <- load_pipeline_config()
  samples <- load_sample_manifest(config)
  if (!"sample_id" %in% names(samples) || !nrow(samples) || anyNA(samples$sample_id) ||
      any(!nzchar(samples$sample_id)) || anyDuplicated(samples$sample_id) ||
      any(basename(samples$sample_id) != samples$sample_id)) stop("Invalid sample manifest.")
  task_map <- data.frame(task = seq_len(nrow(samples)), sample_id = samples$sample_id)
  if (isTRUE(opts[["--list"]])) {
    if (length(opts) != 1L) stop("Use --list alone.")
    write.table(task_map, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  combine <- isTRUE(opts[["--combine"]])
  preflight <- isTRUE(opts[["--dry-run"]]) && !combine && is.null(opts[["--task"]])
  if (!preflight && combine == !is.null(opts[["--task"]])) {
    stop("Choose exactly one of --task N or --combine (or use --dry-run alone).")
  }
  core_y_genes <- c("DDX3Y", "RPS4Y1", "KDM5D", "UTY", "EIF1AY", "USP9Y")
  genes <- c("SRY", core_y_genes, "XIST")
  marker_set_version <- "expanded_Y_v2"
  out_dir <- if (is.null(opts[["--output-dir"]])) {
    here::here(config$project$outputs_dir, "xenium", "sex_marker_expression")
  } else here::here(opts[["--output-dir"]])
  inputs <- file.path(here::here(config$project$outputs_dir), "xenium", "preprocess",
                      "02_qc", "rds", paste0(samples$sample_id, "_CB_QC.rds"))
  per_sample <- file.path(out_dir, "samples", paste0(samples$sample_id, "_sex_markers.csv"))
  combined_outputs <- file.path(out_dir, c("sample_sex_evidence.csv", "gene_expression.csv",
    "sex_marker_expression.pdf", "provenance.txt", "marker_coverage.csv"))
  if (preflight) {
    required <- inputs
    destinations <- c(per_sample, combined_outputs)
  } else if (combine) {
    required <- per_sample
    destinations <- combined_outputs
  } else {
    task <- number("--task", NA_real_, TRUE, nrow(samples))
    required <- inputs[[task]]
    destinations <- per_sample[[task]]
  }
  if (isTRUE(opts[["--dry-run"]])) {
    message("Markers (", marker_set_version, "): ", paste(genes, collapse = ", "))
    ok <- compact_dry_run("sex-marker expression (paths only; no RDS loaded)",
      inputs = required, outputs = destinations,
      checks = c(outputs_protected = !any(file.exists(destinations))))
    if (!ok) stop("Preflight failed; inspect missing inputs or choose a fresh output directory.")
    return(invisible(NULL))
  }
  if (any(!file.exists(required))) stop("Missing required inputs: ",
    paste(required[!file.exists(required)], collapse = ", "))
  if (any(file.exists(destinations))) stop("Outputs already exist; choose a fresh --output-dir: ",
    paste(destinations[file.exists(destinations)], collapse = ", "))
  # Exclusive creation also protects against concurrent duplicate submissions.
  write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    con <- file(path, open = "wx")
    on.exit(close(con))
    write.csv(x, con, row.names = FALSE, na = "NA")
  }
  if (!combine) {
    for (p in c("SeuratObject", "Matrix")) {
      if (!requireNamespace(p, quietly = TRUE)) stop("Missing project dependency: ", p)
    }
    sample_id <- samples$sample_id[[task]]
    message("Task ", task, ": ", sample_id, " <- ", required)
    obj <- readRDS(required)
    if (!inherits(obj, "Seurat") || !"Xenium" %in% SeuratObject::Assays(obj)) {
      stop("Expected a Seurat object with a Xenium assay: ", required)
    }
    cells <- colnames(obj)
    if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
      stop("Empty object or invalid/duplicate cell IDs: ", sample_id)
    }
    metadata <- obj[[]]
    if (!identical(rownames(metadata), cells)) stop("Metadata/cell order mismatch: ", sample_id)
    for (column in intersect(c("sample_id", "sample_name"), names(metadata))) {
      if (anyNA(metadata[[column]]) || any(as.character(metadata[[column]]) != sample_id)) {
        stop("Unexpected sample labels in metadata column ", column, ": ", sample_id)
      }
    }
    assay <- obj[["Xenium"]]
    layers <- SeuratObject::Layers(assay, search = "^counts($|\\.)")
    if (length(layers) != 1L || layers[[1]] != "counts") {
      stop("Expected one unsplit raw counts layer; found: ", paste(layers, collapse = ", "))
    }
    counts <- SeuratObject::LayerData(assay, layer = "counts")
    if (!identical(colnames(counts), cells) || anyDuplicated(rownames(counts)) ||
        anyNA(rownames(counts)) || any(!nzchar(rownames(counts)))) {
      stop("Invalid feature IDs or counts/cell order mismatch: ", sample_id)
    }
    # Check raw counts without converting a sparse matrix into a dense matrix.
    stored <- if (inherits(counts, "sparseMatrix") && "x" %in% slotNames(counts)) {
      counts@x
    } else as.vector(counts)
    if (any(!is.finite(stored)) || any(stored < 0) || any(stored != floor(stored))) {
      stop("Counts contain non-finite, negative, or non-integer values: ", sample_id)
    }
    totals <- Matrix::colSums(counts)
    if (any(!is.finite(totals)) || sum(totals) <= 0) stop("Invalid/empty transcript totals.")
    measured <- genes %in% rownames(counts)
    result <- data.frame(task = task, sample_id = sample_id, gene = genes,
      measured = measured, n_cells = length(cells), total_transcripts = sum(totals),
      gene_transcripts = NA_real_, positive_cells = NA_real_, positive_pct = NA_real_,
      mean_counts_per_cell = NA_real_, transcripts_per_10000 = NA_real_)
    for (j in which(measured)) {
      x <- as.numeric(counts[genes[[j]], , drop = TRUE])
      result$gene_transcripts[[j]] <- sum(x)
      result$positive_cells[[j]] <- sum(x > 0)
      result$positive_pct[[j]] <- 100 * mean(x > 0)
      result$mean_counts_per_cell[[j]] <- mean(x)
      result$transcripts_per_10000[[j]] <- 1e4 * sum(x) / sum(totals)
    }
    result$supported <- ifelse(measured,
      result$positive_cells >= min_positive_cells & result$positive_pct >= min_positive_pct, NA)
    evidence <- classify_sex_marker_evidence(result, core_y_genes, min_cells)
    for (column in names(evidence)) result[[column]] <- evidence[[column]]
    result$marker_set_version <- marker_set_version
    result$min_positive_cells <- min_positive_cells
    result$min_positive_pct <- min_positive_pct
    result$min_cells <- min_cells
    result$input_path <- normalizePath(required)
    info <- file.info(required)
    result$input_bytes <- info$size
    result$input_mtime_utc <- format(info$mtime, tz = "UTC", usetz = TRUE)
    result$analysis_time_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    result$r_version <- R.version.string
    result$seuratobject_version <- as.character(utils::packageVersion("SeuratObject"))
    write_csv(result, destinations)
    message(sample_id, ": ", evidence$provisional_sex, "; ", evidence$reason)
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing project dependency: ggplot2")
  tables <- lapply(seq_along(per_sample), function(k) {
    x <- read.csv(per_sample[[k]], stringsAsFactors = FALSE)
    required_columns <- c("task", "sample_id", "gene", "measured", "n_cells",
      "total_transcripts", "gene_transcripts", "positive_cells", "positive_pct",
      "mean_counts_per_cell", "transcripts_per_10000", "supported", "provisional_sex",
      "reason", "min_positive_cells", "min_positive_pct", "min_cells", "input_path",
      "input_bytes", "input_mtime_utc", "analysis_time_utc", "r_version", "seuratobject_version",
      "marker_set_version", "evidence_pattern", "n_y_markers_measured", "n_y_markers_supported",
      "supported_y_markers", "missing_markers")
    if (length(setdiff(required_columns, names(x))) ||
        anyNA(x[setdiff(required_columns, c("gene_transcripts", "positive_cells", "positive_pct",
          "mean_counts_per_cell", "transcripts_per_10000", "supported"))])) {
      stop("Incomplete sample table: ", per_sample[[k]])
    }
    if (nrow(x) != length(genes) || !identical(x$gene, genes) ||
        !all(x$marker_set_version == marker_set_version) ||
        !all(x$sample_id == samples$sample_id[[k]]) || !all(x$task == k) ||
        !all(x$min_positive_cells == min_positive_cells) ||
        !all(x$min_positive_pct == min_positive_pct) || !all(x$min_cells == min_cells) ||
        anyNA(x$measured) || any(!is.finite(x$n_cells)) || any(x$n_cells <= 0) ||
        any(!is.finite(x$positive_pct[x$measured])) ||
        any(!x$provisional_sex %in% c("male-like", "female-like", "unresolved")) ||
        length(unique(x$provisional_sex)) != 1L) {
      stop("Invalid or inconsistent sample table: ", per_sample[[k]])
    }
    expected_support <- ifelse(x$measured,
      x$positive_cells >= min_positive_cells & x$positive_pct >= min_positive_pct, NA)
    if (!identical(x$supported, expected_support)) stop("Inconsistent marker support: ", per_sample[[k]])
    evidence <- classify_sex_marker_evidence(x, core_y_genes, min_cells)
    for (column in names(evidence)) {
      if (!all(x[[column]] == evidence[[column]])) stop("Inconsistent ", column, ": ", per_sample[[k]])
    }
    x
  })
  all_genes <- do.call(rbind, tables)
  summary_columns <- c("task", "sample_id", "n_cells", "total_transcripts", "provisional_sex",
    "reason", "min_positive_cells", "min_positive_pct", "min_cells", "input_path",
    "input_bytes", "input_mtime_utc", "analysis_time_utc", "r_version", "seuratobject_version",
    "marker_set_version", "evidence_pattern", "n_y_markers_measured", "n_y_markers_supported",
    "supported_y_markers", "missing_markers")
  summary <- all_genes[all_genes$gene == "SRY", summary_columns]
  for (g in genes) {
    x <- all_genes[all_genes$gene == g, ]
    if (!identical(x$sample_id, summary$sample_id)) stop("Sample alignment error.")
    for (col in c("measured", "supported", "positive_cells", "positive_pct", "mean_counts_per_cell")) {
      summary[[paste(g, col, sep = "_")]] <- x[[col]]
    }
  }
  plot_data <- all_genes
  plot_data$gene <- factor(plot_data$gene, levels = genes)
  plot_data$sample_id <- factor(plot_data$sample_id, levels = rev(samples$sample_id))
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(gene, sample_id)) +
    ggplot2::geom_point(data = plot_data[plot_data$measured, ],
      ggplot2::aes(size = positive_pct, colour = log1p(mean_counts_per_cell))) +
    ggplot2::geom_text(data = plot_data[!plot_data$measured, ], label = "not measured", size = 2) +
    ggplot2::scale_size_area(max_size = 8, limits = c(0, 100), name = "% cells >0 counts") +
    ggplot2::scale_colour_viridis_c(name = "log1p(mean counts)") +
    ggplot2::labs(x = NULL, y = NULL, title = "Sample sex-marker expression",
      subtitle = "Raw Xenium counts in whole-tissue QC cells; zero expression has no dot") +
    ggplot2::theme_minimal(base_size = 10)
  write_csv(summary, destinations[[1]])
  write_csv(all_genes, destinations[[2]])
  ggplot2::ggsave(destinations[[3]], plot, width = 12,
    height = max(5, 0.23 * nrow(samples) + 2), limitsize = FALSE)
  writeLines(c(
    "Exploratory expression evidence; provisional calls require review, not confirmed sex.",
    paste("Core Y markers:", paste(core_y_genes, collapse = ", ")),
    paste("Marker set version:", marker_set_version),
    "SRY is supporting evidence; XIST is associated with X inactivation.",
    "Missing SRY expression in cerebellum does not exclude a male sample.",
    "Low Y counts may reflect background. Mixed signals need review; these criteria are not validated.",
    "Core Y support can yield male-like evidence without XIST; single and multiple Y support are distinguished.",
    "Low Y expression without XIST remains unresolved. SRY-only support remains unresolved.",
    "No filtering, normalization, integration, annotation, metadata updates, or RDS writes performed.",
    paste("Manifest:", resolve_config_path(config$manifests$samples, config)),
    paste("Samples:", nrow(samples)),
    paste("Outputs:", paste(destinations, collapse = "; ")),
    paste("Command:", paste(commandArgs(), collapse = " ")),
    capture.output(sessionInfo())
  ), destinations[[4]])
  coverage <- do.call(rbind, lapply(genes, function(g) {
    x <- all_genes[all_genes$gene == g, ]
    data.frame(gene = g, n_samples_measured = sum(x$measured),
      n_samples_missing = sum(!x$measured),
      missing_samples = if (any(!x$measured)) paste(x$sample_id[!x$measured], collapse = ";") else "none")
  }))
  write_csv(coverage, destinations[[5]])
  message("Combined ", nrow(samples), " complete sample tables in ", out_dir)
}

if (sys.nframe() == 0L) main()
