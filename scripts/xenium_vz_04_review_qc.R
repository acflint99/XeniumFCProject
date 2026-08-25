#!/usr/bin/env Rscript

# Review merged VZ clusters, then explicitly apply the reviewed cell removal
# decision to every manifest-defined whole-tissue object.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(future)
  library(future.apply)
  library(ggplot2)
})
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
cluster_column <- "Xenium_snn_res.0.8"
args <- commandArgs(trailingOnly = TRUE)

remove_args <- grep("^--remove-clusters=", args, value = TRUE)
if (length(remove_args) > 1L) stop("Supply --remove-clusters only once.")
remove_supplied <- length(remove_args) == 1L
remove_value <- if (remove_supplied) sub("^--remove-clusters=", "", remove_args) else ""
exact_options <- c("--qc-only", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% exact_options &
                          !grepl("^--remove-clusters=", args)]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!startsWith(args, "--"))) {
  stop("All arguments must be named options. See the usage message in this script.")
}
qc_only <- "--qc-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (qc_only && remove_supplied) stop("Choose --qc-only or --remove-clusters, not both.")

clusters_to_remove <- character()
if (remove_supplied && !identical(tolower(trimws(remove_value)), "none")) {
  clusters_to_remove <- unique(trimws(strsplit(remove_value, ",", fixed = TRUE)[[1]]))
  if (any(!nzchar(clusters_to_remove))) stop("--remove-clusters contains a blank cluster ID.")
}

output_root <- here(config$project$outputs_dir)
merged_path <- file.path(output_root, "XenAld_VZ_Res1.5_RDS", "Xenium_VZ_Res1.5.rds")
input_dir <- file.path(output_root, "Xenium_ConsensusABT_Res1.5_RDS")
output_dir <- file.path(output_root, "Xenium_AldingerABT_VZ_QC_Res1.5_RDS")
plot_dir <- file.path(output_root, "XenAld_VZ_QC_Res1.5_Plots")
table_dir <- file.path(output_root, "XenAld_VZ_QC_Res1.5_Tables")

input_names <- paste0(sample_ids, "_Consensus_annotated.rds")
input_paths <- file.path(input_dir, input_names)
output_names <- paste0(sample_ids, "_Ald_VZ_QC.rds")
output_paths <- file.path(output_dir, output_names)
summary_path <- file.path(table_dir, "XenAld_VZ_RawSubcluster_QC_Summary.csv")
marker_path <- file.path(table_dir, "VZ_RawSubcluster_top10_Markers_Res0.8.csv")
review_path <- file.path(table_dir, "XenAld_VZ_QC_review_manifest.csv")
pdf_path <- file.path(plot_dir, "XenAld_VZ_RawSubcluster_QC_Violins.pdf")
removal_path <- file.path(table_dir, "XenAld_VZ_QC_removal_manifest.csv")
qc_outputs <- c(summary_path, marker_path, review_path, pdf_path)
removal_outputs <- c(output_paths, removal_path)

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  cat("Merged VZ object:", merged_path, "\n")
  cat("Merged VZ object exists:", file.exists(merged_path), "\n")
  cat("Expected whole-tissue inputs:", length(input_paths), "\n")
  cat("Existing whole-tissue inputs:", sum(file.exists(input_paths)), "\n")
  cat("Unexpected top-level RDS files:", length(unexpected_names), "\n")
  cat("QC review outputs existing:", sum(file.exists(qc_outputs)), "of", length(qc_outputs), "\n")
  cat("Filtered sample outputs existing:", sum(file.exists(output_paths)), "of", length(output_paths), "\n")
  cat("Requested mode:", if (qc_only) "QC review" else if (remove_supplied) "apply removal" else "inspection only", "\n")
  if (remove_supplied) {
    cat("Clusters to remove:", if (length(clusters_to_remove)) paste(clusters_to_remove, collapse = ", ") else "none", "\n")
  }
  write.table(
    data.frame(sample_id = sample_ids, input_exists = file.exists(input_paths),
               filtered_output_exists = file.exists(output_paths)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  if (length(unexpected_names)) cat("Unexpected files:\n- ", paste(unexpected_names, collapse = "\n- "), "\n")
  quit(save = "no", status = 0L)
}

if (!qc_only && !remove_supplied) {
  stop(
    "Choose a QC stage explicitly:\n",
    "- --qc-only to create evidence for review\n",
    "- --remove-clusters=7,8 after review\n",
    "- --remove-clusters=none after review when no cells should be removed"
  )
}
if (!file.exists(merged_path)) stop("Merged VZ object not found: ", merged_path)

if (qc_only) {
  existing_qc <- qc_outputs[file.exists(qc_outputs)]
  if (length(existing_qc) && !overwrite) {
    stop("Refusing to overwrite existing VZ QC review outputs:\n- ",
         paste(existing_qc, collapse = "\n- "), "\nUse --overwrite only after review.")
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(merged_path)
  if (!cluster_column %in% colnames(obj[[]])) stop("Merged VZ object lacks ", cluster_column, ".")
  if (ncol(obj) == 0L) stop("Merged VZ object contains zero cells.")
  obj <- JoinLayers(obj)
  Idents(obj) <- cluster_column
  set.seed(config$runtime$random_seed)

  qc_stats <- obj@meta.data %>%
    group_by(.data[[cluster_column]]) %>%
    summarise(cell_count = n(), median_counts = median(nCount_Xenium),
              mean_counts = mean(nCount_Xenium),
              median_features = median(nFeature_Xenium), .groups = "drop")
  write.csv(qc_stats, summary_path, row.names = FALSE)

  grDevices::cairo_pdf(pdf_path, width = 12, height = 8)
  for (feature in c("nCount_Xenium", "nFeature_Xenium")) {
    print(
      VlnPlot(obj, features = feature, group.by = cluster_column, pt.size = 0) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        ylab(feature) + ggtitle(paste("Xenium VZ Subcluster QC -", feature))
    )
  }
  grDevices::dev.off()

  all_markers <- FindAllMarkers(
    obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "wilcox", max.cells.per.ident = 1000, random.seed = config$runtime$random_seed
  )
  top_markers <- all_markers %>% group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE)
  write.csv(top_markers, marker_path, row.names = FALSE)

  merged_info <- file.info(merged_path)
  review_manifest <- data.frame(
    merged_path = merged_path, merged_size = as.numeric(merged_info$size),
    merged_mtime = as.numeric(merged_info$mtime), cells = ncol(obj),
    cluster_column = cluster_column,
    clusters = paste(levels(Idents(obj)), collapse = "|"),
    random_seed = config$runtime$random_seed, stringsAsFactors = FALSE
  )
  write.csv(review_manifest, review_path, row.names = FALSE)
  message("VZ QC review complete. Inspect the PDF, summary, and marker table before removal.")
  quit(save = "no", status = 0L)
}

if (length(missing_names)) stop("Cannot apply VZ QC; missing files:\n- ", paste(missing_names, collapse = "\n- "))
if (length(unexpected_names)) {
  stop("Cannot apply VZ QC; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
missing_review <- qc_outputs[!file.exists(qc_outputs)]
if (length(missing_review)) {
  stop("Run --qc-only and review its outputs first. Missing:\n- ", paste(missing_review, collapse = "\n- "))
}
existing_removal <- removal_outputs[file.exists(removal_outputs)]
if (length(existing_removal) && !overwrite) {
  stop("Refusing to overwrite existing VZ QC-filtered outputs:\n- ",
       paste(existing_removal, collapse = "\n- "), "\nUse --overwrite only after review.")
}

obj <- readRDS(merged_path)
if (!cluster_column %in% colnames(obj[[]])) stop("Merged VZ object lacks ", cluster_column, ".")
review <- read.csv(review_path, stringsAsFactors = FALSE)
merged_info <- file.info(merged_path)
review_is_current <- nrow(review) == 1L &&
  isTRUE(all.equal(as.numeric(review$merged_size), as.numeric(merged_info$size))) &&
  isTRUE(all.equal(as.numeric(review$merged_mtime), as.numeric(merged_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj)))
if (!review_is_current) stop("VZ merged object changed after QC review. Rerun --qc-only before removal.")
available_clusters <- unique(as.character(obj[[cluster_column, drop = TRUE]]))
unknown_clusters <- setdiff(clusters_to_remove, available_clusters)
if (length(unknown_clusters)) stop("Requested VZ cluster(s) do not exist: ", paste(unknown_clusters, collapse = ", "))
cells_to_remove <- if (length(clusters_to_remove)) {
  colnames(obj)[as.character(obj[[cluster_column, drop = TRUE]]) %in% clusters_to_remove]
} else character()

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
options(future.globals.maxSize = 200 * 1024^3)
plan(multisession, workers = workers)
on.exit(plan(sequential), add = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rows <- future_lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  temp_obj <- readRDS(input_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(temp_obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  integrated_barcodes <- paste0(sample_id, "_", colnames(temp_obj))
  remove_cells <- integrated_barcodes %in% cells_to_remove
  removed <- sum(remove_cells)
  temp_obj <- subset(temp_obj, cells = colnames(temp_obj)[!remove_cells])
  saveRDS(temp_obj, output_paths[[i]], compress = FALSE)
  data.frame(sample_id = sample_id, input_cells = length(remove_cells),
             removed_cells = removed, remaining_cells = ncol(temp_obj),
             clusters_removed = if (length(clusters_to_remove)) paste(clusters_to_remove, collapse = "|") else "none",
             stringsAsFactors = FALSE)
}, future.seed = TRUE)
plan(sequential)
removal_manifest <- do.call(rbind, rows)
if (sum(removal_manifest$removed_cells) != length(cells_to_remove)) {
  stop("VZ removal count across individual objects does not match the merged object.")
}
write.csv(removal_manifest, removal_path, row.names = FALSE)
message("Applied reviewed VZ removal decision to all ", length(sample_ids), " samples.")
