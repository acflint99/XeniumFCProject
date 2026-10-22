# Remove only a terminal Ensembl version; preserve separate _PAR_Y records.
reference_stable_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))

panel_reference_alignment <- function(ids, panel) {
  stable <- reference_stable_id(ids)
  relevant <- stable %in% panel$ensembl_id
  if (anyDuplicated(stable[relevant])) stop("Multiple original reference rows match one panel Ensembl ID; explicit alignment required.")
  comparison <- panel$gene[match(stable, panel$ensembl_id)]
  data.frame(original_feature = ids, stable_id = stable, comparison_feature = comparison,
             usable = !is.na(comparison), stringsAsFactors = FALSE)
}

# Original published Science HDF5 reader. No package installs or source writes.
# Metadata preflight: ~1–3 min from a 67-second observed audit. Expression pass:
# 1–6 h or longer, low confidence until HDF5 chunk behavior is measured.
science_h5_vector <- function(input, name) {
  allowed <- c("CellID", "Age", "Region", "Clusters", "CellClass", "Donor", "SampleID",
               "Gene", "Accession", "ClusterID", "Class", "NCells", "TotalUMIs")
  if (!name %in% allowed) stop("Not an allowed Science metadata vector: ", name)
  if (!nzchar(Sys.which("h5dump")) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Science reader needs existing h5dump and jsonlite; no dependencies are installed.")
  tmp <- tempfile("science_vector_"); log <- tempfile("science_header_")
  on.exit(unlink(c(tmp, log)), add = TRUE)
  status <- system2(Sys.which("h5dump"), c("-e", "-y", "-w", "0", "-d", shQuote(paste0("/shoji/", name)),
                                            "-o", shQuote(tmp), shQuote(input)), stdout = log, stderr = log)
  if (status != 0L || !file.exists(tmp)) stop("Science metadata extraction failed: ", name)
  v <- tryCatch(jsonlite::fromJSON(paste0("[", paste(readLines(tmp, warn = FALSE), collapse = "\n"), "]")),
                error = function(e) stop("Cannot safely parse Science ", name, ": ", conditionMessage(e)))
  if (!is.atomic(v) || !is.null(dim(v))) stop("Expected a metadata vector: ", name)
  v
}

plan_science_original <- function(input, block_rows = 128L) {
  if (!is.finite(block_rows) || block_rows < 1 || block_rows != floor(block_rows) || block_rows > 1024L) stop("Science block rows must be an integer from 1 to 1024.")
  text_ok <- function(x) !anyNA(x) && all(nzchar(trimws(as.character(x))))
  ids <- science_h5_vector(input, "CellID"); n <- length(ids)
  if (!n || !text_ok(ids) || anyDuplicated(ids)) stop("Invalid original Science CellIDs.")
  md <- setNames(lapply(c("Age", "Region", "Clusters", "CellClass", "Donor", "SampleID", "TotalUMIs"),
                        function(k) science_h5_vector(input, k)),
                 c("Age", "Region", "Clusters", "CellClass", "Donor", "SampleID", "TotalUMIs"))
  if (any(lengths(md) != n)) stop("Science cell metadata lengths differ.")
  if (!is.numeric(md$Age) || any(!is.finite(md$Age))) stop("Science source Age must be finite numeric PCW values.")
  if (!all(vapply(md[c("Region", "CellClass", "Donor", "SampleID")], text_ok, FALSE))) stop("Missing Science cell metadata.")
  if (!is.numeric(md$Clusters) || any(!is.finite(md$Clusters)) || any(md$Clusters != floor(md$Clusters))) stop("Invalid original Science cluster IDs.")
  if (!is.numeric(md$TotalUMIs) || any(!is.finite(md$TotalUMIs)) || any(md$TotalUMIs < 0) || any(md$TotalUMIs != floor(md$TotalUMIs))) stop("Invalid original Science TotalUMIs.")
  lookup <- data.frame(ClusterID = science_h5_vector(input, "ClusterID"),
                       Class = science_h5_vector(input, "Class"),
                       original_cells = science_h5_vector(input, "NCells"))
  if (anyNA(lookup) || anyDuplicated(lookup$ClusterID) || !text_ok(lookup$Class)) stop("Invalid Science cluster lookup.")
  li <- match(md$Clusters, lookup$ClusterID)
  if (anyNA(li) || any(as.character(md$CellClass) != as.character(lookup$Class[li]))) stop("Original Science cluster/class mapping mismatch.")
  observed <- tabulate(li, nbins = nrow(lookup))
  if (any(observed != lookup$original_cells)) stop("Science original NCells mismatch.")
  # User approved 2026-10-21: explicit cerebellum, continuous 9–15 PCW inclusive.
  keep <- md$Region == "Cerebellum" & md$Age >= 9 & md$Age <= 15
  selected <- which(keep)
  if (!length(selected)) stop("No original Science cells match approved cerebellum/9–15 cohort.")
  cells <- data.frame(source_row_one_based = selected, CellID = ids[selected],
                       as.data.frame(lapply(md, `[`, selected), stringsAsFactors = FALSE), check.names = FALSE)
  rm(md, ids)
  genes <- science_h5_vector(input, "Gene"); accession <- science_h5_vector(input, "Accession")
  if (!length(genes) || length(genes) != length(accession) || !text_ok(genes) || !text_ok(accession)) stop("Invalid Science gene metadata alignment.")
  # Repeated symbols are allowed in source metadata. They are ambiguous only
  # when used to align comparison features; never silently choose the first row.
  duplicate <- duplicated(genes) | duplicated(genes, fromLast = TRUE)
  gene_duplicates <- data.frame(source_gene_index_one_based = which(duplicate),
                                 Gene = genes[duplicate], Accession = accession[duplicate])
  if (any(duplicate)) {
    cat("\nORIGINAL SCIENCE DUPLICATE-SYMBOL AUDIT:\n")
    cat("Repeated symbols:", length(unique(genes[duplicate])), "; source rows:", sum(duplicate), "\n")
    print(head(gene_duplicates, 80L), row.names = FALSE)
    if (nrow(gene_duplicates) > 80L) cat("First 80 duplicate rows shown. The full table is saved with a completed comparison.\n")
    cat("Metadata planning can proceed. Shared ambiguous symbols require panel-ID or explicit accession alignment before expression reading (unless omission is explicitly selected).\n")
  }
  header <- system2(Sys.which("h5dump"), c("-H", "-d", "/shoji/Expression", shQuote(input)), stdout = TRUE, stderr = TRUE)
  status <- attr(header, "status")
  if (!is.null(status) && status != 0L) stop("Cannot inspect original Science Expression header.")
  h <- paste(header, collapse = " ")
  dims <- regmatches(h, regexec("SIMPLE[[:space:]]*\\{[[:space:]]*\\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*\\)", h))[[1L]]
  if (length(dims) != 3L || !identical(as.numeric(dims[2:3]), c(as.numeric(n), as.numeric(length(genes)))) || !grepl("H5T_STD_U16LE", h, fixed = TRUE)) stop("Science Expression must be the audited cells x genes uint16 LE matrix.")
  # Split only consecutive retained cells; never read excluded cells to bridge gaps.
  runs <- split(selected, cumsum(c(TRUE, diff(selected) != 1L)))
  blocks <- unlist(lapply(runs, function(x) split(x, ceiling(seq_along(x) / block_rows))), recursive = FALSE)
  starts <- vapply(blocks, function(x) x[1L], 0L); sizes <- lengths(blocks)
  block_plan <- data.frame(source_start_zero_based = starts - 1L, cells = sizes)
  if (sum(sizes) != nrow(cells)) stop("Science chunk plan lost selected cells.")
  counts <- tabulate(match(cells$Clusters, lookup$ClusterID), nbins = nrow(lookup))
  donors <- vapply(lookup$ClusterID, function(k) length(unique(cells$Donor[cells$Clusters == k])), 0L)
  lookup$selected_cells <- counts; lookup$selected_donors <- donors
  summary <- data.frame(region = "Cerebellum", age_min = 9, age_max = 15,
                         cells = nrow(cells), donors = length(unique(cells$Donor)),
                         libraries = length(unique(cells$SampleID)), clusters = sum(counts > 0),
                         source_genes = length(genes), expression_blocks = nrow(block_plan),
                         max_rows_per_block = max(sizes),
                         uncompressed_selected_GiB = nrow(cells) * length(genes) * 2 / 1024^3)
  cat("\nAPPROVED SCIENCE COHORT / EXPRESSION READ PLAN (metadata only):\n")
  print(summary, row.names = FALSE)
  print(stats::aggregate(list(cells = rep.int(1L, nrow(cells))), cells[c("Age", "Donor", "SampleID")], length), row.names = FALSE)
  if (nrow(block_plan) > 10000L) stop("More than 10,000 disjoint expression chunks; review I/O plan before proceeding.")
  list(cells = cells, genes = genes, accession = accession, gene_duplicates = gene_duplicates, lookup = lookup, blocks = block_plan,
       summary = summary, input = input, source_size = file.info(input)$size,
       source_mtime = as.numeric(file.info(input)$mtime))
}

resolve_science_genes <- function(plan, shared_genes, feature_map = "", duplicate_policy = "error", panel = NULL) {
  if (!duplicate_policy %in% c("error", "omit")) stop("Science duplicate policy must be error or omit.")
  if (nzchar(feature_map) && duplicate_policy != "error") stop("Use an explicit Science feature map OR duplicate omission, not both.")
  if (nzchar(feature_map)) {
    mapping <- read.csv(feature_map, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("reference_gene", "query_gene") %in% names(mapping))) stop("Science map CSV requires reference_gene (original Accession) and query_gene columns.")
    good_text <- function(x) !anyNA(x) && all(nzchar(trimws(as.character(x))))
    if (!nrow(mapping) || !good_text(mapping$reference_gene) || !good_text(mapping$query_gene) ||
        anyDuplicated(mapping$reference_gene) || anyDuplicated(mapping$query_gene)) stop("Science feature map must be nonempty, nonblank, and one-to-one.")
    if (any(!mapping$reference_gene %in% plan$accession)) stop("Science feature map contains accessions absent from the source.")
    ambiguous_accession <- unique(plan$accession[duplicated(plan$accession)])
    if (any(mapping$reference_gene %in% ambiguous_accession)) stop("Mapped Science accessions themselves are duplicated in the original source.")
    mapping <- mapping[mapping$query_gene %in% shared_genes, , drop = FALSE]
    gi <- match(mapping$reference_gene, plan$accession)
    keep_genes <- mapping$query_gene
  } else if (!is.null(panel)) {
    if (duplicate_policy != "error") stop("Panel ID alignment cannot be combined with duplicate-symbol omission.")
    a <- panel_reference_alignment(plan$accession, panel)
    gi <- which(a$usable & a$comparison_feature %in% shared_genes)
    keep_genes <- a$comparison_feature[gi]
  } else {
    ambiguous <- unique(plan$genes[duplicated(plan$genes)])
    shared_ambiguous <- intersect(shared_genes, ambiguous)
    if (length(shared_ambiguous) && duplicate_policy == "error") stop("Ambiguous original Science symbols shared with Xenium: ",
       paste(head(shared_ambiguous, 30L), collapse = ", "),
       ". Supply --science-feature-map=CSV with original Accession in reference_gene and Xenium symbol in query_gene, or explicitly select --science-duplicate-policy=omit after reviewing the scientific effect. No duplicate rows are combined or selected automatically.")
    if (length(shared_ambiguous) && duplicate_policy == "omit") message("Explicitly omitting ambiguous Science comparison symbols: ", paste(shared_ambiguous, collapse = ", "))
    keep_genes <- intersect(setdiff(shared_genes, shared_ambiguous), plan$genes)
    gi <- match(keep_genes, plan$genes)
  }
  if (length(keep_genes) < 3L) stop("Fewer than three unambiguous shared Science/Xenium genes.")
  data.frame(gene = keep_genes, original_gene_index_one_based = gi,
             original_symbol = plan$genes[gi], accession = plan$accession[gi], stringsAsFactors = FALSE)
}

aggregate_science_original <- function(plan, shared_genes, output_dir, feature_map = "", duplicate_policy = "error", panel = NULL) {
  if (file.info(plan$input)$size != plan$source_size || as.numeric(file.info(plan$input)$mtime) != plan$source_mtime) stop("Original Science input changed after metadata preflight.")
  alignment <- resolve_science_genes(plan, shared_genes, feature_map, duplicate_policy, panel)
  if (!is.null(panel) && !nzchar(feature_map)) {
    audit <- panel_reference_alignment(plan$accession, panel)
    audit$original_symbol <- plan$genes
    audit$selected_for_comparison <- seq_len(nrow(audit)) %in% alignment$original_gene_index_one_based
    write.csv(audit, file.path(output_dir, "science_accession_alignment_audit.csv"), row.names = FALSE)
    coverage <- panel
    coverage$matched_source <- coverage$gene %in% alignment$gene
    write.csv(coverage, file.path(output_dir, "science_panel_id_coverage.csv"), row.names = FALSE)
  }
  keep_genes <- alignment$gene; gi <- alignment$original_gene_index_one_based
  cluster_ids <- sort(plan$lookup$ClusterID[plan$lookup$selected_cells > 0])
  sums <- matrix(0, length(keep_genes), length(cluster_ids), dimnames = list(keep_genes, as.character(cluster_ids)))
  cells <- plan$cells; n_genes <- length(plan$genes); seen <- integer(nrow(cells))
  write.csv(cells, file.path(output_dir, "science_selected_cells.csv"), row.names = FALSE)
  write.csv(plan$lookup, file.path(output_dir, "science_original_cluster_selection.csv"), row.names = FALSE)
  write.csv(plan$blocks, file.path(output_dir, "science_expression_read_plan.csv"), row.names = FALSE)
  write.csv(plan$summary, file.path(output_dir, "science_cohort_summary.csv"), row.names = FALSE)
  write.csv(alignment, file.path(output_dir, "science_gene_alignment.csv"), row.names = FALSE)
  write.csv(plan$gene_duplicates, file.path(output_dir, "science_original_duplicate_symbols.csv"), row.names = FALSE)
  omitted <- if (duplicate_policy == "omit") plan$gene_duplicates[plan$gene_duplicates$Gene %in% shared_genes, , drop = FALSE] else plan$gene_duplicates[FALSE, , drop = FALSE]
  write.csv(omitted, file.path(output_dir, "science_omitted_ambiguous_genes.csv"), row.names = FALSE)
  tmp <- tempfile("science_count_block_", fileext = ".bin"); log <- tempfile("science_block_log_")
  on.exit(unlink(c(tmp, log)), add = TRUE)
  for (i in seq_len(nrow(plan$blocks))) {
    start <- plan$blocks$source_start_zero_based[i]; nr <- plan$blocks$cells[i]
    ci <- match(seq.int(start + 1L, length.out = nr), cells$source_row_one_based)
    if (anyNA(ci)) stop("Science expression chunk would include an excluded cell.")
    # h5dump exports a bounded hyperslab in C row-major order, uint16 little endian.
    status <- system2(Sys.which("h5dump"), c("-d", "/shoji/Expression", "-s", paste0(start, ",0"),
      "-c", paste0(nr, ",", n_genes), "-b", "LE", "-o", shQuote(tmp), shQuote(plan$input)), stdout = log, stderr = log)
    expected_bytes <- as.double(nr) * n_genes * 2
    if (status != 0L || !file.exists(tmp) || file.info(tmp)$size != expected_bytes) stop("Science binary extraction failed or unexpected byte count at chunk ", i)
    v <- readBin(tmp, what = integer(), n = nr * n_genes, size = 2L, signed = FALSE, endian = "little")
    m <- matrix(v, nrow = nr, ncol = n_genes, byrow = TRUE)
    # A full-gene sum check detects transpose, byte-order, and raw-count inconsistencies.
    if (any(rowSums(m) != cells$TotalUMIs[ci])) stop("Science Expression row totals differ from original TotalUMIs in chunk ", i, "; stop for source/count semantics review.")
    by_cluster <- rowsum(m[, gi, drop = FALSE], group = cells$Clusters[ci], reorder = FALSE)
    jj <- match(rownames(by_cluster), colnames(sums))
    sums[, jj] <- sums[, jj, drop = FALSE] + t(by_cluster)
    seen[ci] <- seen[ci] + 1L
    if (i == 1L || i %% 10L == 0L || i == nrow(plan$blocks)) message("Science expression chunks: ", i, "/", nrow(plan$blocks))
  }
  if (any(seen != 1L)) stop("Not every retained Science cell was read exactly once.")
  if (file.info(plan$input)$size != plan$source_size || as.numeric(file.info(plan$input)$mtime) != plan$source_mtime) stop("Science source changed during expression reading.")
  list(counts = sums, cell_counts = data.frame(reference_label = as.character(cluster_ids),
         cells = plan$lookup$selected_cells[match(cluster_ids, plan$lookup$ClusterID)]))
}
