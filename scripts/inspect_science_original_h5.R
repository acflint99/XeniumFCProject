#!/usr/bin/env Rscript
# Read-only HDF5 inventory and bounded previews. Never reads a complete matrix.
# Expected HPC runtime: 1–5 min, low confidence (schema/storage dependent).
# --dry-run: path-only check, under 10 seconds; no HDF5 contents loaded.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript --max-ppsize=500000 scripts/inspect_science_original_h5.R\n",
      "       [--file=PATH] [--reader=auto|hdf5r|rhdf5|h5dump] [--dry-run]\n",
      "Prints HDF5 structure, identifier/annotation previews, and matrix storage metadata.\n",
      "Source opened read-only; no output files written. Full matrices never loaded.\n", sep = "")
  quit(save = "no")
}
file_args <- args[startsWith(args, "--file=")]
reader_args <- args[startsWith(args, "--reader=")]
if (length(file_args) > 1L || length(reader_args) > 1L || any(!args %in% c(file_args, reader_args, "--dry-run"))) stop("Unknown/repeated option; use --help.")
requested_reader <- if (length(reader_args)) sub("^--reader=", "", reader_args) else "auto"
if (!requested_reader %in% c("auto", "hdf5r", "rhdf5", "h5dump")) stop("Invalid --reader; use --help.")
if (length(file_args)) {
  input <- sub("^--file=", "", file_args)
} else {
  for (p in c("here", "yaml")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing project dependency: ", p)
  here::i_am("scripts/inspect_science_original_h5.R")
  source(here::here("scripts", "R", "config.R"))
  config <- load_pipeline_config()
  input <- file.path(config$inputs$published_data_root, "ScienceBraunFC", "Science_HumanFetalBrainPool.h5")
}
if (!nzchar(input)) stop("Empty input path.")
cat("Original Science candidate:", input, "\n")
cat("Exists:", file.exists(input), "\n")
available <- c(hdf5r = requireNamespace("hdf5r", quietly = TRUE),
               rhdf5 = requireNamespace("rhdf5", quietly = TRUE),
               h5dump = nzchar(Sys.which("h5dump")))
cat("Available readers:", paste(names(available), available, sep = "=", collapse = "; "), "\n")
reader <- if (requested_reader == "auto") names(available)[available][1L] else requested_reader
reader_ok <- length(reader) == 1L && !is.na(reader) && isTRUE(available[[reader]])
if ("--dry-run" %in% args) {
  cat("Path/dependency check; HDF5 contents not loaded. Selected reader:", if (reader_ok) reader else "NONE", "\n")
  quit(save = "no", status = if (file.exists(input) && reader_ok) 0L else 1L)
}
if (!file.exists(input)) stop("Input not found; run on the HPC where this source is available.")
if (!reader_ok) stop("No selected HDF5 reader is available. Readers tried: hdf5r, rhdf5, h5dump. Check existing HPC HDF5 tools with 'module spider HDF5' and share the output. No packages are installed and renv is unchanged.")
cat("Selected reader:", reader, "\n")
cat("File size (GiB):", round(file.info(input)$size / 1024^3, 3), "\n")

inspect <- function(input) {
  f <- hdf5r::H5File$new(input, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  # Enumerate metadata breadth-first, with explicit limits for unusual files.
  queue <- list(list(group = f, path = "", depth = 0L))
  rows <- list(); groups <- list("/" = f); nodes <- 0L; truncated <- FALSE
  while (length(queue)) {
    entry <- queue[[1L]]; queue <- queue[-1L]
    listing <- entry$group$ls(recursive = FALSE)
    if (!nrow(listing)) next
    for (i in seq_len(nrow(listing))) {
      nodes <- nodes + 1L
      if (nodes > 1000L) { truncated <- TRUE; break }
      p <- paste0(entry$path, "/", listing$name[i])
      is_group <- as.character(listing$obj_type[i]) == "H5I_GROUP"
      dims <- if ("dataset.dims" %in% names(listing)) as.character(listing$dataset.dims[i]) else ""
      rows[[length(rows) + 1L]] <- data.frame(path = p, type = if (is_group) "group" else as.character(listing$obj_type[i]), dimensions = dims)
      if (is_group) {
        g <- entry$group[[listing$name[i]]]; groups[[p]] <- g
        if (entry$depth < 5L) queue[[length(queue) + 1L]] <- list(group = g, path = p, depth = entry$depth + 1L) else truncated <- TRUE
      }
    }
    if (nodes > 1000L) break
  }
  inventory <- if (length(rows)) do.call(rbind, rows) else data.frame(path = character(), type = character(), dimensions = character())
  cat("\n--- FILE STRUCTURE (metadata only; at most 1000 entries, depth 6) ---\n")
  print(inventory, row.names = FALSE, max = 4000)
  if (truncated) cat("Inventory truncated by entry/depth limit; unlisted fields may exist.\n")

  preview <- function(p, n = 8L) {
    d <- f[[sub("^/", "", p)]]
    if (!inherits(d, "H5D")) return(invisible(NULL))
    dims <- d$dims
    cat("\n", p, " dimensions=", paste(dims, collapse = " x "), "\n", sep = "")
    if (!length(dims) || any(dims == 0L)) { cat("Empty/scalar dataset; no preview.\n"); return(invisible(NULL)) }
    if (length(dims) != 1L) { cat("Multidimensional dataset: contents not read.\n"); return(invisible(NULL)) }
    ix <- seq_len(min(n, dims[1L]))
    v <- tryCatch(d$read(args = list(ix)), error = function(e) paste("Preview unavailable:", conditionMessage(e)))
    print(v)
    if (dims[1L] > n) cat("Preview only; remaining entries not read.\n")
    invisible(NULL)
  }
  cat("\n--- FORMAT / MATRIX ATTRIBUTES ---\n")
  # AnnData shape and encoding attributes distinguish CSR/CSC and table layouts.
  for (p in names(groups)) {
    if (p != "/" && !grepl("^/(X|raw|layers|matrix|obs|var)(/|$)", p)) next
    g <- groups[[p]]
    attrs <- hdf5r::h5attr_names(g)
    for (a in intersect(attrs, c("encoding-type", "encoding-version", "shape", "_index"))) {
      cat(p, "attribute", a, ": ")
      print(hdf5r::h5attr(g, a))
    }
  }
  cat("\n--- CELL IDS, GENE IDS, AND ANNOTATION PREVIEWS ---\n")
  data_paths <- inventory$path[inventory$type == "H5I_DATASET"]
  metadata_paths <- data_paths[
    grepl("^/(obs|var|col_attrs|row_attrs)(/|$)", data_paths) |
    grepl("(barcode|gene|feature|cell.?type|annotation|cluster|donor|sample|label|_index)", data_paths, ignore.case = TRUE)]
  # Exclude numerical matrix storage even when nested under an annotated group.
  metadata_paths <- metadata_paths[!grepl("/(data|indices|indptr)$", metadata_paths)]
  for (p in head(metadata_paths, 150L)) {
    n <- if (grepl("/(categories|levels)$", p)) 100L else 8L
    preview(p, n)
  }
  if (length(metadata_paths) > 150L) cat("Metadata previews limited to 150 datasets.\n")
  cat("\n--- MATRIX STORAGE PREVIEWS (at most 8 stored values per vector) ---\n")
  storage <- data_paths[grepl("^/(X|raw/X|layers/[^/]+|matrix)(/(data|indices|indptr|shape))?$", data_paths)]
  for (p in head(storage, 30L)) preview(p)
  if (length(storage) > 30L) cat("Storage previews limited to 30 datasets.\n")
  cat("\n--- INTERPRETATION ---\n",
      "obs/col_attrs may contain cell annotations; var/row_attrs or matrix/features may contain gene IDs.\n",
      "Categorical annotations may store integer codes plus a separate categories/levels vector.\n",
      "matrix/barcodes and matrix/features with data/indices/indptr often indicate a 10x count matrix.\n",
      "X, raw/X, and layers are candidate matrices, not proof that values are original raw counts.\n",
      "Integer-looking previews do not establish raw-count provenance.\n",
      "If annotations are absent, the original barcode-to-label metadata file is also required.\n",
      "Previews do not validate complete ID uniqueness, metadata joins, or all annotation values.\n",
      "No original labels/IDs are changed. No expression matrix is reconstructed or fully read.\n", sep = "")
}
inspect_rhdf5 <- function(input) {
  # h5ls reads object metadata only, bounded to depth 6; vectors are hypersliced.
  listing <- rhdf5::h5ls(input, recursive = 6)
  inventory <- head(listing, 1000L)
  inventory$path <- paste0(ifelse(inventory$group == "/", "", inventory$group), "/", inventory$name)
  cat("\n--- FILE STRUCTURE (metadata only; up to 1000 entries, depth 6) ---\n")
  print(inventory[, c("path", "otype", "dclass", "dim")], row.names = FALSE, max = 5000)
  if (nrow(listing) > 1000L) cat("Inventory truncated to 1000 entries.\n")
  cat("Depth limit: deeper fields may exist and are not enumerated.\n")
  cat("\n--- FORMAT / MATRIX ATTRIBUTES ---\n")
  groups <- c("/", inventory$path[inventory$otype == "H5I_GROUP"])
  for (p in groups) {
    if (p != "/" && !grepl("^/(X|raw|layers|matrix|obs|var)(/|$)", p)) next
    attrs <- rhdf5::h5readAttributes(input, p)
    for (a in intersect(names(attrs), c("encoding-type", "encoding-version", "shape", "_index"))) {
      cat(p, "attribute", a, ": "); print(attrs[[a]])
    }
  }
  datasets <- inventory[inventory$otype == "H5I_DATASET", , drop = FALSE]
  preview <- function(i, n = 8L) {
    row <- datasets[i, , drop = FALSE]
    dims <- suppressWarnings(as.numeric(strsplit(as.character(row$dim), "x", fixed = TRUE)[[1L]]))
    cat("\n", row$path, " dimensions=", as.character(row$dim), "\n", sep = "")
    if (length(dims) != 1L || !is.finite(dims) || dims < 1L) {
      cat("Empty, scalar, or multidimensional dataset: contents not read.\n"); return(invisible(NULL))
    }
    result <- tryCatch(rhdf5::h5read(input, row$path, index = list(seq_len(min(n, dims)))),
                       error = function(e) paste("Preview unavailable:", conditionMessage(e)))
    print(result)
    if (dims > n) cat("Preview only; remaining entries not read.\n")
  }
  cat("\n--- CELL IDS, GENE IDS, AND ANNOTATION PREVIEWS ---\n")
  candidates <- which((grepl("^/(obs|var|col_attrs|row_attrs)(/|$)", datasets$path) |
    grepl("(barcode|gene|feature|cell.?type|annotation|cluster|donor|sample|label|_index)", datasets$path, ignore.case = TRUE)) &
    !grepl("/(data|indices|indptr)$", datasets$path))
  for (i in head(candidates, 150L)) preview(i, if (grepl("/(categories|levels)$", datasets$path[i])) 100L else 8L)
  cat("\n--- MATRIX STORAGE PREVIEWS (at most 8 stored values per vector) ---\n")
  storage <- which(grepl("^/(X|raw/X|layers/[^/]+|matrix)(/(data|indices|indptr|shape))?$", datasets$path))
  for (i in head(storage, 30L)) preview(i)
  cat("\nPreviews do not prove raw-count provenance or validate complete cell IDs/joins.\n",
      "Categorical labels may use codes with categories; annotations may require a separate original metadata file.\n",
      "Source opened for reading only; full expression matrices are never read or reconstructed.\n", sep = "")
}

if (reader == "hdf5r") {
  inspect(input)
} else if (reader == "rhdf5") {
  inspect_rhdf5(input)
} else {
  # System HDF5 header dump is a minimal fallback, with no dataset values read.
  cat("\n--- HDF5 HEADERS ONLY (dataset values not read) ---\n")
  report <- system2(Sys.which("h5dump"), c("-H", shQuote(input)), stdout = TRUE, stderr = TRUE)
  status <- attr(report, "status")
  if (!is.null(status) && status != 0L) stop("h5dump failed: ", paste(head(report, 20L), collapse = "\n"))
  cat(paste(head(report, 1500L), collapse = "\n"), "\n")
  if (length(report) > 1500L) cat("Header listing truncated at 1500 lines.\n")
  cat("This fallback shows field names, shapes, and types only. Share this output to identify targeted fields for value previews.\n")
}
