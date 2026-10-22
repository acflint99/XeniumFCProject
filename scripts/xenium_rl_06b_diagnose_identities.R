#!/usr/bin/env Rscript

# Read-only Stage-5/Stage-6 identity review: within-cell marker detection,
# cluster QC, and original cell-level predictions from all three references.
rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_rl_identity_diagnostics.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite"))) {
  stop("Usage: xenium_rl_06b_diagnose_identities.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
root <- here(config$project$outputs_dir, "xenium", "rl")
input <- file.path(root, "05_post_qc", "rds", "Xenium_RL_postQC.rds")
review_path <- file.path(root, "06_subclusters", "tables",
                         "Xenium_RL_subcluster_review_manifest_Res0.8.csv")
cluster_column <- "Xenium_snn_res.0.8"
focus <- c("5", "6", "10", "12", "18", "19")
references <- c("Aldinger", "Sepp", "Science")
transfer_root <- get_production_output_paths(config)$label_transfer_dir
task_map <- expand.grid(sample_id = samples, reference = references,
                        stringsAsFactors = FALSE)
task_map$input_path <- file.path(
  transfer_root, tolower(task_map$reference), "rds",
  paste0(task_map$sample_id, "_", task_map$reference, "_annotated.rds")
)
diagnostic_root <- file.path(root, "06_subclusters", "diagnostics", "identity_review")
table_dir <- file.path(diagnostic_root, "tables")
plot_dir <- file.path(diagnostic_root, "plots")
names_tables <- c(
  "Cluster_QC.csv", "Cluster_BySample_QC.csv", "Marker_Availability.csv",
  "Marker_Pairs.csv", "Coexpression_ByCluster.csv", "Coexpression_BySample.csv",
  "Focus_Cell_Marker_Counts.csv.gz", "Reference_Label_Composition.csv",
  "Reference_BySample_Composition.csv", "Focus_Cell_Reference_Predictions.csv.gz",
  "Reference_Input_Join_Audit.csv", "Provenance.csv"
)
table_paths <- setNames(file.path(table_dir, paste0("Xenium_RL_", names_tables)), names_tables)
plot_paths <- file.path(plot_dir, paste0("Xenium_RL_", c(
  "Cluster_QC.pdf", "Marker_Coexpression.pdf", "Reference_Composition.pdf"
)))
outputs <- c(unname(table_paths), plot_paths)
inputs <- c(input, review_path, task_map$input_path)
if (dry_run) {
  compact_dry_run("RL identity diagnostics [resolution 0.8, read-only]",
                  inputs = inputs, outputs = outputs)
  cat("  Path-only: no RDS objects loaded. Reference inputs: ", nrow(task_map),
      "; samples: ", length(samples), ".\n", sep = "")
  quit(save = "no", status = if (all(file.exists(inputs))) 0L else 1L)
}
if (any(!file.exists(inputs))) stop("Missing input(s):\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
existing <- outputs[file.exists(outputs)]
if (length(existing) && !overwrite) {
  stop("Refusing to overwrite diagnostic outputs:\n", paste(existing, collapse = "\n"),
       "\nUse --overwrite only after review.")
}
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2)
})
review <- read.csv(review_path, stringsAsFactors = FALSE)
required_review <- c("input_rds", "input_size", "input_mtime", "cells", "cluster_column", "clusters")
info <- file.info(input)
if (nrow(review) != 1L || !all(required_review %in% names(review)) ||
    !identical(review$cluster_column, cluster_column) ||
    !identical(review$input_rds, normalizePath(input)) ||
    !isTRUE(all.equal(review$input_size, as.numeric(info$size))) ||
    !isTRUE(all.equal(review$input_mtime, as.numeric(info$mtime)))) {
  stop("Stage-6 review provenance does not match the current post-QC object.")
}
message(Sys.time(), ": loading RL post-QC object...")
obj <- readRDS(input)
meta <- obj[[]]
required_meta <- c("orig.ident", cluster_column, "PCW", "nCount_Xenium", "nFeature_Xenium")
if (!all(required_meta %in% names(meta))) stop("Required QC/sample/cluster metadata is missing.")
if (anyDuplicated(colnames(obj)) || !identical(rownames(meta), colnames(obj))) {
  stop("Post-QC cell IDs are duplicated or metadata order is inconsistent.")
}
meta <- data.frame(
  cell_id = rownames(meta), sample_id = as.character(meta$orig.ident),
  cluster = as.character(meta[[cluster_column]]), PCW = as.character(meta$PCW),
  nCount_Xenium = meta$nCount_Xenium, nFeature_Xenium = meta$nFeature_Xenium,
  stringsAsFactors = FALSE
)
if (anyNA(meta) || any(!nzchar(meta$cluster)) || any(!nzchar(meta$PCW)) ||
    any(!is.finite(meta$nCount_Xenium)) || any(!is.finite(meta$nFeature_Xenium)) ||
    any(meta$nCount_Xenium <= 0) || any(meta$nFeature_Xenium <= 0) ||
    !setequal(meta$sample_id, samples) || nrow(meta) != review$cells ||
    !setequal(meta$cluster, strsplit(review$clusters, "|", fixed = TRUE)[[1]])) {
  stop("Post-QC metadata, sample coverage, or reviewed cell/cluster set is invalid.")
}
pairs <- data.frame(
  gene1 = c("EOMES", "EOMES", "EOMES", "EOMES", "SOX2", "EOMES", "EOMES",
            "CALB1", "FOXP2", "FOXP4", "CALB1", "CALB1", "GAD1", "EOMES", "OTX2", "EOMES"),
  gene2 = c("SLC1A3", "SOX2", "SOX9", "BCAN", "SLC1A3", "RSPO3", "TRPC3",
            "GAD1", "GAD1", "GAD1", "FOXP2", "FOXP4", "GAD2", "OTX2", "NEUROD1", "NEUROD1"),
  stringsAsFactors = FALSE
)
genes <- unique(c(unlist(pairs), "NES", "FABP7", "PTPRZ1", "HOPX", "GDF10",
                  "COL2A1", "CALCB", "PAX6", "CNTN2", "PRPH", "GRIA2",
                  "SFRP1", "CXCR4", "NEFM", "GABRA6", "SLC17A7",
                  "ITPR1", "DAB1", "SKOR2", "RORA", "PCP4", "CA8", "PAX2"))
if (!"Xenium" %in% Assays(obj)) stop("Post-QC object lacks Xenium assay.")
# Extract only the small targeted count matrix; do not join/copy the full assay.
layer_names <- SeuratObject::Layers(obj[["Xenium"]])
if (!"counts" %in% layer_names) stop("Expected an exact joined counts layer in Stage-5 object.")
available_genes <- genes[genes %in% rownames(obj[["Xenium"]])]
counts <- SeuratObject::LayerData(obj[["Xenium"]], layer = "counts", features = available_genes)
if (!identical(colnames(counts), meta$cell_id) ||
    !setequal(rownames(counts), available_genes)) stop("Targeted counts do not align with cells/genes.")
if (any(!is.finite(counts@x)) || any(counts@x < 0)) stop("Targeted counts contain invalid values.")
rm(obj); invisible(gc(full = TRUE))
message(Sys.time(), ": summarizing QC and within-cell detection...")
qc <- rl_diagnostic_qc(meta, "cluster")
qc_sample <- rl_diagnostic_qc(meta, c("sample_id", "cluster"))
sample_baseline <- rl_diagnostic_qc(meta, "sample_id") |>
  select(sample_id, sample_median_counts = median_counts,
         sample_median_features = median_features)
qc_sample <- left_join(qc_sample, sample_baseline, by = "sample_id")
if (anyNA(qc_sample$sample_median_counts)) stop("Within-sample QC baseline join failed.")
qc_sample$median_counts_ratio_to_sample <- qc_sample$median_counts / qc_sample$sample_median_counts
qc_sample$median_features_ratio_to_sample <- qc_sample$median_features / qc_sample$sample_median_features
availability <- data.frame(gene = genes, available = genes %in% available_genes)
coexpression <- rl_diagnostic_pairs(counts, meta, pairs, "cluster")
coexpression_sample <- rl_diagnostic_pairs(counts, meta, pairs, c("sample_id", "cluster"))
focus_index <- which(meta$cluster %in% focus)
if (!length(focus_index)) stop("None of the requested focus clusters are present.")
cell_counts <- cbind(meta[focus_index, , drop = FALSE],
                     as.data.frame(as.matrix(t(counts[, focus_index, drop = FALSE]))))
rm(counts); invisible(gc(full = TRUE))

reference_summary <- reference_sample <- reference_focus <- vector("list", length(references))
join_audit <- vector("list", nrow(task_map))
for (r in seq_along(references)) {
  reference <- references[r]
  labels <- rep(NA_character_, nrow(meta)); scores <- rep(NA_real_, nrow(meta))
  for (sample_id in samples) {
    task <- which(task_map$reference == reference & task_map$sample_id == sample_id)
    if (length(task) != 1L) stop("Reference task map is not one-to-one.")
    path <- task_map$input_path[task]
    message(Sys.time(), ": reading ", reference, " / ", sample_id)
    original <- readRDS(path)
    original_meta <- original[[]]
    if (!identical(rownames(original_meta), colnames(original)) ||
        !all(c("predicted.id", "prediction.score.max") %in% names(original_meta))) {
      stop("Original prediction columns or cell alignment are missing: ", path)
    }
    index <- which(meta$sample_id == sample_id)
    matched <- rl_diagnostic_match_cells(meta$cell_id[index], rownames(original_meta), sample_id)
    selected_labels <- as.character(original_meta$predicted.id[matched])
    selected_scores <- original_meta$prediction.score.max[matched]
    if (anyNA(selected_labels) || any(!nzchar(trimws(selected_labels))) ||
        !is.numeric(selected_scores) || any(!is.finite(selected_scores)) ||
        any(selected_scores < 0 | selected_scores > 1)) {
      stop("Original predictions have blank labels or invalid scores: ", path)
    }
    labels[index] <- selected_labels; scores[index] <- selected_scores
    source_info <- file.info(path)
    join_audit[[task]] <- data.frame(
      sample_id = sample_id, reference = reference, retained_rl_cells = length(index),
      matched_cells = length(matched), unmatched_cells = 0L,
      original_cells = nrow(original_meta), input_path = normalizePath(path),
      input_size = as.numeric(source_info$size), input_mtime = as.numeric(source_info$mtime)
    )
    rm(original, original_meta); invisible(gc(full = TRUE))
  }
  if (anyNA(labels) || anyNA(scores)) stop("Incomplete reference coverage: ", reference)
  frame <- meta[, c("cell_id", "sample_id", "cluster"), drop = FALSE]
  frame$reference <- reference; frame$predicted_label <- labels; frame$prediction_score <- scores
  reference_summary[[r]] <- rl_diagnostic_reference(frame, c("reference", "cluster"))
  reference_sample[[r]] <- rl_diagnostic_reference(frame, c("reference", "sample_id", "cluster"))
  reference_focus[[r]] <- frame[focus_index, , drop = FALSE]
  rm(frame, labels, scores); invisible(gc(full = TRUE))
}
reference_summary <- bind_rows(reference_summary)
reference_sample <- bind_rows(reference_sample)
reference_focus <- bind_rows(reference_focus)
join_audit <- bind_rows(join_audit)
if (nrow(join_audit) != nrow(task_map) ||
    any(join_audit$matched_cells != join_audit$retained_rl_cells)) stop("Reference audit is incomplete.")
for (reference in references) {
  if (sum(reference_summary$cells[reference_summary$reference == reference]) != nrow(meta)) {
    stop("Reference summary cell totals are inconsistent: ", reference)
  }
}
cluster_levels <- sort(unique(as.integer(meta$cluster)))
if (anyNA(cluster_levels)) stop("Expected numeric resolution-0.8 cluster IDs.")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
tables <- list(qc, qc_sample, availability, pairs, coexpression, coexpression_sample,
               cell_counts, reference_summary, reference_sample, reference_focus, join_audit)
for (i in seq_along(tables)) {
  path <- table_paths[i]
  if (endsWith(path, ".gz")) {
    connection <- gzfile(path, "wt")
    tryCatch(write.csv(tables[[i]], connection, row.names = FALSE), finally = close(connection))
  } else write.csv(tables[[i]], path, row.names = FALSE)
}
# All plotting is descriptive and uses every retained cell; no random sampling.
qc$cluster <- factor(qc$cluster, levels = cluster_levels)
p_qc <- bind_rows(
  transmute(qc, cluster, metric = "Median transcripts", value = median_counts),
  transmute(qc, cluster, metric = "Median detected genes", value = median_features),
  transmute(qc, cluster, metric = "Cells", value = cells)
) |>
  ggplot(aes(cluster, value)) + geom_col() + facet_wrap(~metric, scales = "free_y", ncol = 1) +
  theme_bw() + labs(x = "Resolution-0.8 cluster", y = NULL)
coexpression$pair <- paste(coexpression$gene1, coexpression$gene2, sep = " + ")
coexpression$cluster <- factor(coexpression$cluster, levels = cluster_levels)
p_pairs <- ggplot(coexpression, aes(pair, cluster, fill = 100 * both_fraction)) +
  geom_tile() + scale_fill_viridis_c(limits = c(0, 100), name = "% both detected") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = NULL, y = "Resolution-0.8 cluster", subtitle = "Detection = at least one Xenium transcript for each gene in the same cell")
reference_summary$cluster <- factor(reference_summary$cluster, levels = cluster_levels)
p_reference <- ggplot(reference_summary, aes(cluster, fraction_of_cluster, fill = predicted_label)) +
  geom_col() + facet_wrap(~reference, ncol = 1) + theme_bw() +
  labs(x = "Resolution-0.8 cluster", y = "Fraction of cells", fill = "Original prediction")
for (i in seq_along(plot_paths)) {
  ggsave(plot_paths[i], plot = list(p_qc, p_pairs, p_reference)[[i]],
         device = grDevices::cairo_pdf, width = 12, height = 9)
}
provenance <- data.frame(
  input_rds = normalizePath(input), input_size = as.numeric(info$size),
  input_mtime = as.numeric(info$mtime), review_manifest = normalizePath(review_path),
  cluster_column = cluster_column, cells = nrow(meta), samples = length(samples),
  clusters = paste(cluster_levels, collapse = "|"),
  focus_clusters_present = paste(intersect(focus, as.character(cluster_levels)), collapse = "|"),
  reference_inputs = nrow(task_map), detection_rule = "raw Xenium counts > 0; descriptive only",
  prediction_labels = "original predicted.id; no confidence filtering or consensus replacement",
  prediction_scores = "original prediction.score.max; interpret within each reference",
  source_cell_matching = "exact Stage-2 sample_id prefix; one-to-one match required",
  random_sampling = FALSE, analysis_objects_modified = FALSE,
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write.csv(provenance, table_paths[["Provenance.csv"]], row.names = FALSE)
if (any(!file.exists(outputs)) || any(file.info(outputs)$size <= 0)) stop("Partial diagnostic output set.")
message("RL identity diagnostics complete: ", diagnostic_root,
        ". No analysis RDS files, cluster identities, or QC decisions were changed.")
