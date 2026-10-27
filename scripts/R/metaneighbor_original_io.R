# Original-source access for MetaNeighborUS. All source objects are read-only.
mn_text <- function(x) length(x) > 0L && !anyNA(x) && all(nzchar(trimws(as.character(x))))
mn_sce <- function(x) !is.data.frame(attr(x, "meta.data", exact = TRUE)) && inherits(x, "SummarizedExperiment")
mn_metadata <- function(x) {
  m <- attr(x, "meta.data", exact = TRUE)
  if (!is.data.frame(m)) {
    if (!mn_sce(x)) stop("Expected an original Seurat or SummarizedExperiment object.")
    m <- as.data.frame(SummarizedExperiment::colData(x))
  }
  if (!nrow(m) || !mn_text(rownames(m)) || anyDuplicated(rownames(m))) stop("Invalid metadata cell IDs.")
  m
}
mn_panel <- function(file, expected = NULL) {
  d <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  targets <- Filter(function(x) identical(x$type$descriptor, "gene"), d$payload$targets)
  get_field <- function(k) vapply(targets, function(x) {
    a <- x$type$data[[k]]
    if (!is.character(a) || length(a) != 1L || is.na(a)) stop("Invalid panel annotation: ", k)
    a
  }, "")
  p <- data.frame(gene = get_field("name"), ensembl_id = get_field("id"))
  if (!mn_text(p$gene) || anyDuplicated(p$gene) || anyDuplicated(p$ensembl_id) ||
      any(!grepl("^ENSG[0-9]{11}$", p$ensembl_id))) stop("Panel names and unversioned Ensembl IDs must be nonblank and unique.")
  if (!is.null(expected) && !setequal(expected, p$gene)) stop("Original panel differs from configured panel names.")
  p
}
mn_cohort <- function(md, ref) {
  if (ref == "aldinger") {
    if (!"age" %in% names(md) || !mn_text(md$age) || any(!grepl("^[0-9]+(\\.[0-9]+)? PCW$", md$age))) stop("Original Aldinger age must contain explicit PCW strings.")
    age <- as.numeric(sub(" PCW$", "", md$age)); keep <- age >= 9 & age <= 21
  } else if (ref == "sepp") {
    stages <- c("7 wpc", "8 wpc", "9 wpc", "11 wpc", "17 wpc", "20 wpc")
    tissues <- c("cerebellar cx", "cerebellar dn", "cerebellar flocculonodular lobe", "cerebellar vermis", "cerebellum")
    if (!all(c("Stage", "Tissue") %in% names(md)) || !mn_text(md$Stage) || !mn_text(md$Tissue) ||
        any(!md$Stage %in% c(stages, "newborn", "infant", "toddler", "adult")) || any(!md$Tissue %in% tissues)) stop("Original Sepp Stage/Tissue differs from the approved inventory.")
    keep <- md$Stage %in% stages
  } else stop("Unsupported original RDS reference.")
  if (!any(keep)) stop("Approved cohort is empty.")
  keep
}
mn_raw_blocks <- function(obj, assay) {
  md <- mn_metadata(obj)
  if (mn_sce(obj)) {
    if (!assay %in% SummarizedExperiment::assayNames(obj)) stop("Missing original assay: ", assay)
    m <- SummarizedExperiment::assay(obj, assay)
    if (!identical(colnames(m), rownames(md))) stop("Original SCE cells do not align to metadata.")
    return(list(counts = m))
  }
  assays <- attr(obj, "assays", exact = TRUE)
  if (!assay %in% names(assays)) stop("Missing original assay: ", assay)
  a <- assays[[assay]]; m <- attr(a, "counts", exact = TRUE)
  if (!is.null(m)) return(list(counts = m))
  layers <- SeuratObject::Layers(a, search = "^counts($|\\.)")
  if (!length(layers)) stop("No raw count layers.")
  setNames(lapply(layers, function(layer) SeuratObject::LayerData(a, layer = layer)), layers)
}
mn_extract_counts <- function(obj, assay, cells, panel, ensembl = FALSE) {
  md <- mn_metadata(obj)
  if (!mn_text(cells) || anyDuplicated(cells) || anyNA(match(cells, rownames(md)))) stop("Invalid selected cells.")
  blocks <- mn_raw_blocks(obj, assay)
  original_genes <- rownames(blocks[[1L]])
  if (!mn_text(original_genes) || anyDuplicated(original_genes)) stop("Original feature IDs are invalid.")
  alignment <- if (ensembl) panel_reference_alignment(original_genes, panel) else
    data.frame(original_feature = original_genes, stable_id = original_genes,
               comparison_feature = ifelse(original_genes %in% panel$gene, original_genes, NA_character_),
               usable = original_genes %in% panel$gene)
  gi <- which(alignment$usable)
  if (length(gi) < 5L || anyDuplicated(alignment$comparison_feature[gi])) stop("Fewer than five one-to-one panel genes.")
  seen <- integer(nrow(md)); pieces <- list(); piece_ids <- character()
  for (layer in names(blocks)) {
    m <- blocks[[layer]]
    if (!identical(rownames(m), original_genes) || anyDuplicated(colnames(m))) stop("Original count-layer features/cells differ.")
    ix <- match(colnames(m), rownames(md))
    if (anyNA(ix)) stop("Count-layer cells absent from metadata.")
    seen[ix] <- seen[ix] + 1L
    k <- colnames(m) %in% cells
    if (!any(k)) next
    m <- methods::as(Matrix::Matrix(m[gi, k, drop = FALSE], sparse = TRUE), "dgCMatrix")
    if (any(!is.finite(m@x)) || any(m@x < 0) || any(m@x != round(m@x))) stop("Original counts must be finite, nonnegative integers.")
    pieces[[length(pieces) + 1L]] <- Matrix::drop0(m)
    piece_ids <- c(piece_ids, colnames(m))
  }
  if (any(seen != 1L)) stop("Every original cell must occur in exactly one raw count layer.")
  result <- do.call(cbind, pieces)
  if (anyDuplicated(piece_ids) || !setequal(piece_ids, cells)) stop("Selected count cells are incomplete/duplicated.")
  result <- result[, match(cells, colnames(result)), drop = FALSE]
  rownames(result) <- alignment$comparison_feature[gi]
  list(counts = result, alignment = alignment)
}
mn_science_counts <- function(plan, panel, output_dir) {
  a <- resolve_science_genes(plan, panel$gene, panel = panel)
  pieces <- vector("list", nrow(plan$blocks)); cells <- plan$cells
  tmp <- tempfile("mn_science_", fileext = ".bin"); log <- tempfile("mn_science_log_")
  on.exit(unlink(c(tmp, log)), add = TRUE)
  seen <- integer(nrow(cells)); ng <- length(plan$genes)
  for (i in seq_len(nrow(plan$blocks))) {
    start <- plan$blocks$source_start_zero_based[i]; nr <- plan$blocks$cells[i]
    ci <- match(seq.int(start + 1L, length.out = nr), cells$source_row_one_based)
    if (anyNA(ci)) stop("Science count block contains excluded cells.")
    status <- system2(Sys.which("h5dump"), c("-d", "/shoji/Expression", "-s", paste0(start, ",0"),
      "-c", paste0(nr, ",", ng), "-b", "LE", "-o", shQuote(tmp), shQuote(plan$input)), stdout = log, stderr = log)
    if (status != 0L || !file.exists(tmp) || file.info(tmp)$size != as.double(nr) * ng * 2) stop("Science extraction failed at block ", i)
    v <- readBin(tmp, integer(), n = nr * ng, size = 2L, signed = FALSE, endian = "little")
    m <- matrix(v, nr, ng, byrow = TRUE)
    if (any(rowSums(m) != cells$TotalUMIs[ci])) stop("Original Science count totals differ.")
    pieces[[i]] <- Matrix::Matrix(t(m[, a$original_gene_index_one_based, drop = FALSE]), sparse = TRUE)
    colnames(pieces[[i]]) <- cells$CellID[ci]
    seen[ci] <- seen[ci] + 1L
    if (i == 1L || i %% 20L == 0L || i == nrow(plan$blocks)) message("Science cell count blocks: ", i, "/", nrow(plan$blocks))
  }
  if (any(seen != 1L)) stop("Science selected cell coverage is not exactly once.")
  if (file.info(plan$input)$size != plan$source_size || as.numeric(file.info(plan$input)$mtime) != plan$source_mtime) stop("Original Science source changed during reading.")
  result <- methods::as(do.call(cbind, pieces), "dgCMatrix")
  result <- result[, match(cells$CellID, colnames(result)), drop = FALSE]
  rownames(result) <- a$gene
  write.csv(a, file.path(output_dir, "science_gene_alignment.csv"), row.names = FALSE)
  write.csv(plan$lookup, file.path(output_dir, "science_original_cluster_selection.csv"), row.names = FALSE)
  write.csv(plan$blocks, file.path(output_dir, "science_read_plan.csv"), row.names = FALSE)
  result
}
mn_rankable_cells <- function(m, dataset) {
  m <- Matrix::drop0(methods::as(m, "dgCMatrix"))
  nnz <- diff(m@p)
  reason <- rep(NA_character_, ncol(m))
  reason[nnz == 0L] <- "all_zero_across_selected_HVGs"
  full <- which(nnz == nrow(m))
  for (j in full) {
    values <- m@x[seq.int(m@p[j] + 1L, m@p[j + 1L])]
    if (length(unique(values)) == 1L) reason[j] <- "constant_nonzero_across_selected_HVGs"
  }
  data.frame(study = dataset, source_cell_id = colnames(m), rankable = is.na(reason),
             exclusion_reason = reason, stringsAsFactors = FALSE, check.names = FALSE)
}
