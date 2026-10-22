#!/usr/bin/env Rscript
# Metadata-only audit: 2–10 min on HPC, low confidence for 1.67 million cells.
# --dry-run checks paths/tools/packages only (<10 sec). Expression is never read.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript --max-ppsize=500000 scripts/audit_science_shoji_metadata.R\n",
      "       [--file=ORIGINAL_H5] [--out=NEW_DIRECTORY] [--dry-run]\n",
      "Reads original metadata using existing h5dump and jsonlite. No expression data read.\n",
      "Writes cohort count tables only, with no region/age/cell filtering.\n", sep = "")
  quit(save = "no")
}
opts <- list()
for (a in args) {
  if (a == "--dry-run") next
  if (!grepl("^--(file|out)=.+$", a)) stop("Unknown option: ", a)
  k <- sub("=.*$", "", substring(a, 3L))
  if (!is.null(opts[[k]])) stop("Repeated option: ", k)
  opts[[k]] <- sub("^[^=]+=", "", a)
}
if (is.null(opts$file) || is.null(opts$out)) {
  for (p in c("here", "yaml")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing project dependency: ", p)
  here::i_am("scripts/audit_science_shoji_metadata.R")
  source(here::here("scripts", "R", "config.R"))
  config <- load_pipeline_config()
  if (is.null(opts$file)) opts$file <- file.path(config$inputs$published_data_root, "ScienceBraunFC", "Science_HumanFetalBrainPool.h5")
  if (is.null(opts$out)) opts$out <- here::here(config$project$outputs_dir, "references", "science", "original_metadata_audit")
}
checks <- c(input_exists = file.exists(opts$file), new_output_directory = !file.exists(opts$out),
            h5dump_available = nzchar(Sys.which("h5dump")), jsonlite_available = requireNamespace("jsonlite", quietly = TRUE))
print(data.frame(check = names(checks), pass = unname(checks)), row.names = FALSE)
cat("Input:", opts$file, "\nOutput:", opts$out, "\n")
if ("--dry-run" %in% args) {
  cat("No HDF5 data read; output directory not created.\n")
  quit(save = "no", status = if (all(checks)) 0L else 1L)
}
if (!all(checks)) stop("Preflight failed. No packages are installed; choose a fresh output directory and resolve listed checks.")
started <- Sys.time()
read_vector <- function(name) {
  # Allowlist prevents accidental extraction of the enormous expression matrices.
  allowed <- c("CellID", "Age", "Region", "Subregion", "Subdivision", "Tissue", "Donor", "SampleID", "Clusters", "CellClass", "ClusterID", "Class", "NCells")
  if (!name %in% allowed) stop("Dataset is not an allowed metadata vector: ", name)
  data_path <- paste0("/shoji/", name)
  message("Reading metadata: ", data_path)
  tmp <- tempfile("shoji_metadata_", fileext = ".txt")
  ddl <- tempfile("shoji_ddl_", fileext = ".txt")
  on.exit(unlink(c(tmp, ddl)), add = TRUE)
  status <- system2(Sys.which("h5dump"), c("-e", "-y", "-w", "0", "-d", shQuote(data_path),
                                           "-o", shQuote(tmp), shQuote(opts$file)), stdout = ddl, stderr = ddl)
  if (status != 0L || !file.exists(tmp)) stop("h5dump failed for ", data_path, ": ", paste(head(readLines(ddl, warn = FALSE), 10L), collapse = "\n"))
  text <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  result <- tryCatch(jsonlite::fromJSON(paste0("[", text, "]")),
                     error = function(e) stop("Cannot safely parse metadata export for ", data_path, ": ", conditionMessage(e)))
  if (!is.atomic(result) || !is.null(dim(result))) stop("Expected one-dimensional metadata: ", data_path)
  result
}
ids <- read_vector("CellID")
if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Missing, blank, or duplicate original CellID values.")
ncells <- length(ids); rm(ids)
fields <- c("Age", "Region", "Subregion", "Subdivision", "Tissue", "Donor", "SampleID", "Clusters", "CellClass")
md <- setNames(lapply(fields, read_vector), fields)
if (any(lengths(md) != ncells)) stop("Metadata vector lengths differ from CellID count.")
md <- as.data.frame(md, stringsAsFactors = FALSE)
if (!is.numeric(md$Age) || any(!is.finite(md$Age))) stop("Age must contain finite numeric values. Units still require independent verification.")
cluster_ids <- read_vector("ClusterID"); classes <- read_vector("Class"); n_per_cluster <- read_vector("NCells")
if (!is.numeric(cluster_ids) || !is.numeric(md$Clusters) || anyNA(cluster_ids) || anyDuplicated(cluster_ids) ||
    length(classes) != length(cluster_ids) || length(n_per_cluster) != length(cluster_ids)) stop("Invalid cluster metadata cardinality/types.")
idx <- match(md$Clusters, cluster_ids)
if (anyNA(idx)) stop("Cell cluster IDs absent from original ClusterID table.")
if (anyNA(md$CellClass) || anyNA(classes) || any(as.character(md$CellClass) != as.character(classes[idx]))) stop("CellClass does not match original Class through ClusterID; resolve mapping before analysis.")
observed <- tabulate(idx, nbins = length(cluster_ids))
if (any(!is.finite(n_per_cluster)) || any(observed != n_per_cluster)) stop("Observed per-cluster cells differ from original NCells.")
# This is an audit, not an analysis filter. Explicitly report blank/missing fields.
missing <- data.frame(field = names(md), missing = vapply(md, function(x) sum(is.na(x) | !nzchar(trimws(as.character(x)))), 0L))
key_missing <- "<MISSING_AUDIT_ONLY>"
if (any(vapply(md, function(x) any(as.character(x) == key_missing, na.rm = TRUE), FALSE))) stop("Audit missing-value marker collides with source data.")
for (k in names(md)) {
  if (is.character(md[[k]])) {
    bad <- is.na(md[[k]]) | !nzchar(trimws(md[[k]]))
    md[[k]][bad] <- key_missing
  }
}
counts <- function(columns) {
  # Aggregate only observed combinations; no Cartesian expansion of all factors.
  key <- md[, columns, drop = FALSE]
  result <- stats::aggregate(list(cells = rep.int(1L, ncells)), by = key, FUN = length)
  result[do.call(order, result[columns]), , drop = FALSE]
}
if (!dir.create(opts$out, recursive = TRUE)) stop("Cannot create new audit directory.")
write_table <- function(x, name) write.csv(x, file.path(opts$out, name), row.names = FALSE)
for (k in c("Region", "Subregion", "Subdivision", "Tissue", "Age")) {
  t <- counts(c(k, if (k != "Age") "Age"))
  write_table(t, paste0(k, "_age_counts.csv"))
  cat("\n", k, " values (up to 100 shown):\n", sep = "")
  print(head(t, 100L), row.names = FALSE)
}
cohorts <- counts(c("Region", "Subregion", "Subdivision", "Tissue", "Age", "Donor", "SampleID"))
write_table(cohorts, "cohort_counts.csv")
# Name-based candidates are for review only, never an automatic anatomical filter.
candidate <- Reduce(`|`, lapply(cohorts[c("Region", "Subregion", "Subdivision", "Tissue")], function(x) grepl("cerebell|(^|[^a-z])cb([^a-z]|$)", x, ignore.case = TRUE)))
write_table(cohorts[candidate, , drop = FALSE], "cerebellum_name_candidates.csv")
cat("\nCerebellum-name candidates (heuristic for review, not filtering):\n")
print(cohorts[candidate, , drop = FALSE], row.names = FALSE, max = 2000)
write_table(data.frame(ClusterID = cluster_ids, Class = classes, original_NCells = n_per_cluster, observed_cells = observed), "original_cluster_labels.csv")
write_table(counts(c("Clusters", "CellClass", "Region", "Subregion", "Age")), "cluster_region_age_counts.csv")
write_table(missing, "missing_metadata_counts.csv")
writeLines(c(paste("Input:", normalizePath(opts$file)), paste("Bytes:", file.info(opts$file)$size),
             paste("Modified:", file.info(opts$file)$mtime), paste("Cells:", ncells), paste("Clusters:", length(cluster_ids)),
             paste("Started:", started), paste("Finished:", Sys.time()),
             "Metadata only; no Expression/MeanExpression read; no cells filtered; no labels inferred.",
             "Age values are source values: no PCW conversion or age filter has been applied.",
             "Missing annotation fields are counted and shown with an audit-only marker; resolve before filtering.",
             "Cerebellum-name candidates require review; abbreviations or unexpected region names may be missed.",
             capture.output(sessionInfo())), file.path(opts$out, "provenance.txt"))
writeLines("Metadata audit complete; review cohort_counts and cerebellum_name_candidates before selecting cells.", file.path(opts$out, "COMPLETED.txt"))
message("Metadata audit complete: ", opts$out)
