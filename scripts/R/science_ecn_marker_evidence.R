# Descriptive evidence only; no identity calls, filtering, or object writes.
# Expected runtime: 1-10 min after loading (low-confidence static estimate).
science_ecn_paths <- function(resolution = "1.2") {
  if (length(resolution) != 1L || is.na(resolution) ||
      !grepl("^[0-9]+([.][0-9]+)?$", as.character(resolution))) stop("Invalid Science resolution tag.")
  root <- here::here("outputs", "references", "science", "diagnostics", "raw_cluster_ecn")
  stats::setNames(file.path(root, paste0("Science_res", resolution, "_eCN_", c(
    "DotPlot.pdf", "availability.csv", "cluster_expression.csv",
    "sample_cluster_expression.csv", "coexpression.csv", "provenance.csv",
    "interpretation.txt"
  ))), c("plot", "availability", "pooled", "sample", "pairs", "provenance", "notes"))
}

check_science_ecn_outputs <- function(paths, overwrite = FALSE) {
  if (any(file.exists(paths)) && !overwrite) {
    stop("eCN diagnostic outputs exist; review before using --overwrite: ", dirname(paths[[1]]))
  }
}

write_science_ecn_evidence <- function(object, input_path, overwrite = FALSE, resolution = "1.2") {
  paths <- science_ecn_paths(resolution)
  check_science_ecn_outputs(paths, overwrite)
  metadata <- object[[]]
  cluster_column <- paste0("science_cluster_res", resolution)
  required <- c(cluster_column, "sample_id")
  if (!nrow(metadata) || !all(required %in% names(metadata)) ||
      anyDuplicated(colnames(object)) || !identical(rownames(metadata), colnames(object))) {
    stop("eCN evidence requires unique aligned cells and metadata: ", paste(required, collapse = ", "))
  }
  for (column in required) {
    if (anyNA(metadata[[column]]) || any(!nzchar(trimws(as.character(metadata[[column]]))))) {
      stop("Missing/blank eCN grouping metadata: ", column)
    }
  }
  metadata$cluster <- as.character(metadata[[cluster_column]])
  metadata$sample_id <- as.character(metadata$sample_id)
  panels <- list(
    "eCN core" = c("SLC17A6", "MEIS2"),
    "eCN subtype / developmental" = c("TBR1", "LHX9", "LMX1A", "BARHL1", "NEUROD6"),
    "eCN supporting (shared)" = c("NTNG1", "ZNF804A", "NR4A2", "NRP1", "C1QL3", "NTNG2", "CREG2", "NEUROD2"),
    "GABA context" = c("GAD1", "GAD2", "SLC32A1"),
    "Other excitatory context" = c("SLC17A7", "PAX6", "ATOH1", "EOMES")
  )
  availability <- data.frame(panel = rep(names(panels), lengths(panels)),
                             gene = unlist(panels, use.names = FALSE))
  if (!"originalexp" %in% SeuratObject::Assays(object)) stop("Missing originalexp assay.")
  availability$present <- availability$gene %in% rownames(object[["originalexp"]])
  features <- availability$gene[availability$present]
  if (!length(features)) stop("No eCN/context marker symbols found in originalexp; check feature naming.")
  if (any(!availability$present)) message("Unavailable eCN/context genes: ",
    paste(availability$gene[!availability$present], collapse = ", "))
  # Work on an assay copy; never normalize or modify the caller's object.
  assay <- subset(object[["originalexp"]], features = features)
  # Legacy Assay::Layers matches search literally; Assay5 supports regex.
  # Enumerate first and filter ourselves so both formats use the same lookup.
  available_layers <- SeuratObject::Layers(assay, search = NA)
  matrices <- lapply(c("counts", "data"), function(layer) {
    pattern <- paste0("^", layer, "($|[.])")
    layers <- available_layers[grepl(pattern, available_layers)]
    if (!length(layers)) stop("Missing originalexp ", layer, " layer. Available layers: ",
                              paste(available_layers, collapse = ", "))
    local_assay <- assay
    if (length(layers) > 1L) local_assay <- SeuratObject::JoinLayers(
      local_assay, layers = pattern, new = layer)
    x <- SeuratObject::LayerData(local_assay, layer = if (length(layers) > 1L) layer else layers[[1]])
    if (anyDuplicated(colnames(x)) || !setequal(colnames(x), rownames(metadata)) ||
        !all(features %in% rownames(x))) stop("Incomplete eCN layer cells/features: ", layer)
    x <- x[features, rownames(metadata), drop = FALSE]
    values <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
    if (any(!is.finite(values)) || any(values < 0)) stop("Invalid eCN expression values: ", layer)
    x
  })
  detected <- matrices[[1]] > 0
  linear <- expm1(matrices[[2]])
  if (any(!is.finite(if (inherits(linear, "sparseMatrix")) linear@x else as.vector(linear)))) {
    stop("Non-finite linearized originalexp expression.")
  }
  summarize <- function(columns) {
    keys <- unique(metadata[, columns, drop = FALSE])
    key <- function(frame) do.call(paste, c(lapply(columns, function(column)
      match(frame[[column]], unique(metadata[[column]]))), sep = ":"))
    index <- match(key(metadata), key(keys))
    if (anyNA(index)) stop("Incomplete eCN grouping map.")
    design <- Matrix::sparseMatrix(i = seq_len(nrow(metadata)), j = index, x = 1,
                                   dims = c(nrow(metadata), nrow(keys)))
    n <- as.numeric(Matrix::colSums(design))
    out <- keys[rep(seq_len(nrow(keys)), each = length(features)), , drop = FALSE]
    out$gene <- rep(features, nrow(keys))
    out$cells <- rep(n, each = length(features))
    out$mean_linear <- as.vector(as.matrix(linear %*% design)) / out$cells
    out$log1p_mean <- log1p(out$mean_linear)
    out$pct_detected <- 100 * as.vector(as.matrix(detected %*% design)) / out$cells
    out
  }
  pooled <- summarize("cluster")
  by_sample <- summarize(c("sample_id", "cluster"))
  pair_available <- all(c("SLC17A6", "MEIS2") %in% features)
  pair_frame <- metadata[, c("sample_id", "cluster")]
  pair_frame$cells <- 1L
  pair_frame$both_cells <- if (pair_available) as.integer(
    as.numeric(detected["SLC17A6", ]) > 0 & as.numeric(detected["MEIS2", ]) > 0) else NA_integer_
  pairs <- aggregate(pair_frame[, c("cells", "both_cells")],
                     pair_frame[, c("sample_id", "cluster")], sum)
  pairs$available <- pair_available
  pairs$gene1 <- "SLC17A6"
  pairs$gene2 <- "MEIS2"
  pairs$pct_both <- 100 * pairs$both_cells / pairs$cells
  if (sum(pairs$cells) != nrow(metadata)) stop("eCN cell counts do not reconcile.")
  cluster_levels <- unique(metadata$cluster)
  numbers <- suppressWarnings(as.numeric(cluster_levels))
  cluster_levels <- cluster_levels[if (all(is.finite(numbers))) order(numbers) else order(cluster_levels)]
  pooled$cluster <- factor(pooled$cluster, levels = rev(cluster_levels))
  pooled$gene <- factor(pooled$gene, levels = features)
  pooled$panel <- factor(availability$panel[match(pooled$gene, availability$gene)], levels = names(panels))
  plot <- ggplot2::ggplot(pooled, ggplot2::aes(gene, cluster)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_detected, color = log1p_mean)) +
    ggplot2::facet_grid(cols = ggplot2::vars(panel), scales = "free_x", space = "free_x") +
    ggplot2::scale_size_area(max_size = 6, limits = c(0, 100), name = "% raw count > 0") +
    ggplot2::scale_color_gradient(low = "grey90", high = "#B2182B", name = "log1p(mean\nnormalized expression)") +
    ggplot2::labs(title = "Reprocessed Science: eCN evidence in every raw cluster",
                  x = NULL, y = paste0("Raw resolution-", resolution, " cluster"),
                  caption = "Descriptive evidence; shared/developmental markers do not establish eCN identity.") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  dir.create(dirname(paths[[1]]), recursive = TRUE, showWarnings = FALSE)
  for (name in c("availability", "pooled", "sample", "pairs")) {
    write.csv(switch(name, availability = availability, pooled = pooled,
                     sample = by_sample, pairs = pairs), paths[[name]], row.names = FALSE)
  }
  ggplot2::ggsave(paths[["plot"]], plot, width = 14,
                  height = max(6, length(cluster_levels) * 0.25 + 3), limitsize = FALSE)
  write.csv(data.frame(input_path = input_path, cells = nrow(metadata),
    samples = length(unique(metadata$sample_id)), clusters = length(cluster_levels),
    cluster_column = cluster_column, resolution = resolution, assay = "originalexp",
    layers = "counts;data", detection = "raw count > 0", created = as.character(Sys.time())),
    paths[["provenance"]], row.names = FALSE)
  writeLines(c(
    paste0("All raw Science resolution-", resolution, " clusters are evaluated before broad-label mapping; no labels are changed."),
    "Panels reuse the project's canonical eCN diagnostic markers; SLC17A6/MEIS2 coexpression is supporting evidence, not a classifier.",
    "Subtype/developmental markers need not be universal; supporting genes and excitatory context genes are shared with other lineages.",
    "Absent genes are unavailable evidence, not zero expression. No thresholds or significance tests are applied.",
    "Means are cell-pooled existing log-normalized originalexp expression, linearized before averaging; no normalization is rerun by this check.",
    "Sample-cluster tables and within-cell SLC17A6/MEIS2 detection support replication review. Absent sample-cluster groups are not imputed.",
    "Coexpression cannot exclude segmentation artifacts/doublets. Inspect expression, detection, and sample counts together."
  ), paths[["notes"]])
  if (any(is.na(file.info(paths)$size) | file.info(paths)$size <= 0)) stop("Incomplete eCN diagnostic outputs.")
  message("Completed Science raw-cluster eCN evidence: ", dirname(paths[[1]]))
  invisible(paths)
}
