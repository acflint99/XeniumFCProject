#!/usr/bin/env Rscript
# Lightweight synthetic regression check: under 1 minute; no real RDS inputs.
# Run from scripts/ using the existing R environment; installs no packages.
original <- normalizePath("xenium_sex_marker_expression.R", mustWork = TRUE)
config_original <- normalizePath("config/config.yml", mustWork = TRUE)
manifest_original <- normalizePath("config/samples.csv", mustWork = TRUE)
source(original)
fixture <- tempfile("expanded_sex_marker_test_")
dir.create(fixture)
dir.create(file.path(fixture, "config"))
invisible(file.copy(config_original, file.path(fixture, "config/config.yml")))
invisible(file.copy(manifest_original, file.path(fixture, "config/samples.csv")))
invisible(file.create(file.path(fixture, "fixture.R")))
setwd(fixture)
here::i_am("fixture.R")
samples <- read.csv("config/samples.csv", stringsAsFactors = FALSE)
stopifnot(nrow(samples) >= 16L)
input_dir <- file.path(fixture, "outputs/xenium/preprocess/02_qc/rds")
dir.create(input_dir, recursive = TRUE)
genes <- c("SRY", "DDX3Y", "RPS4Y1", "KDM5D", "UTY", "EIF1AY", "USP9Y", "XIST")
for (k in seq_len(nrow(samples))) {
  n <- if (k == 16L) 10000L else 100L
  counts <- matrix(0, nrow = 10, ncol = n,
    dimnames = list(c(genes, "ACTB", "TUBB"), paste0("cell", seq_len(n))))
  counts["ACTB", ] <- 20
  counts["TUBB", ] <- 5
  if (k == 1L) counts["XIST", 1:20] <- 1
  if (k == 2L) counts["DDX3Y", 1:20] <- 1
  if (k == 3L) counts[c("DDX3Y", "XIST"), 1:20] <- 1
  if (k == 4L) counts <- counts[rownames(counts) != "XIST", ]
  if (k == 6L) counts["DDX3Y", 1:9] <- 1
  if (k == 7L) counts["SRY", 1:20] <- 1
  if (k == 8L) counts <- counts[, 1:50]
  if (k == 9L) {
    counts[c("DDX3Y", "RPS4Y1"), 1:20] <- 1
    counts <- counts[rownames(counts) != "XIST", ]
  }
  if (k == 10L) {
    counts["DDX3Y", 1:20] <- 1
    counts <- counts[rownames(counts) != "XIST", ]
  }
  if (k == 11L) {
    counts["RPS4Y1", 1:20] <- 1
    counts <- counts[!rownames(counts) %in% c("DDX3Y", "XIST"), ]
  }
  if (k == 12L) {
    counts["XIST", 1:20] <- 1
    counts <- counts[rownames(counts) %in% c("XIST", "ACTB", "TUBB"), ]
  }
  if (k == 13L) counts[c("SRY", "XIST"), 1:20] <- 1
  if (k == 14L) counts <- counts[c("ACTB", "TUBB"), , drop = FALSE]
  if (k == 15L) counts["KDM5D", 1:10] <- 1
  if (k == 16L) counts["DDX3Y", 1:10] <- 1
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "Xenium")
  saveRDS(object, file.path(input_dir, paste0(samples$sample_id[k], "_CB_QC.rds")))
}
# Override arguments inside this isolated test process to exercise the driver.
args <- character()
commandArgs <- function(trailingOnly = FALSE) {
  if (trailingOnly) args else c(paste0("--file=", original))
}
args <- "--dry-run"
main()
for (k in seq_len(nrow(samples))) {
  args <- c("--task", as.character(k))
  main()
}
args <- "--combine"
main()
output_dir <- "outputs/xenium/sex_marker_expression"
x <- read.csv(file.path(output_dir, "sample_sex_evidence.csv"), stringsAsFactors = FALSE)
expected <- c("female-like", "male-like", rep("unresolved", 6),
              rep("male-like", 3), rep("unresolved", 3), "male-like", "unresolved")
stopifnot(nrow(x) == nrow(samples), identical(x$sample_id, samples$sample_id),
  identical(x$provisional_sex[1:16], expected),
  x$evidence_pattern[9] == "multiple_Y_markers",
  x$n_y_markers_supported[9] == 2L,
  x$evidence_pattern[10] == "single_Y_marker",
  x$evidence_pattern[4] == "low_Y_XIST_unavailable",
  x$evidence_pattern[7] == "SRY_only",
  x$evidence_pattern[12] == "XIST_without_Y_coverage",
  x$evidence_pattern[14] == "no_core_Y_coverage")
coverage <- read.csv(file.path(output_dir, "marker_coverage.csv"))
stopifnot(identical(coverage$gene, genes),
  all(coverage$n_samples_measured + coverage$n_samples_missing == nrow(samples)),
  coverage$n_samples_missing[coverage$gene == "XIST"] == 5L)
rows <- read.csv(file.path(output_dir, "gene_expression.csv"))
stopifnot(nrow(rows) == 8L * nrow(samples),
  all(is.na(rows$positive_pct[!rows$measured])),
  all(is.na(rows$supported[!rows$measured])))
args <- c("--task", "1")
stopifnot(inherits(try(main(), silent = TRUE), "try-error"))
args <- "--combine"
stopifnot(inherits(try(main(), silent = TRUE), "try-error"))
# An incompatible table must fail combination before any report is created.
bad_dir <- file.path(fixture, "incompatible")
dir.create(file.path(bad_dir, "samples"), recursive = TRUE)
paths <- list.files(file.path(output_dir, "samples"), full.names = TRUE)
invisible(file.copy(paths, file.path(bad_dir, "samples")))
bad_path <- file.path(bad_dir, "samples", basename(paths[[1]]))
bad <- read.csv(bad_path)
bad$marker_set_version <- "original_three_markers"
write.csv(bad, bad_path, row.names = FALSE)
args <- c("--combine", "--output-dir", bad_dir)
stopifnot(inherits(try(main(), silent = TRUE), "try-error"),
  !file.exists(file.path(bad_dir, "sample_sex_evidence.csv")))
cat("PASS: expanded marker calls, coverage, manifest order, incompatible-table rejection, and overwrite protection.\n")
