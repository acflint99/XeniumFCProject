#!/usr/bin/env Rscript

# Merge Xenium cluster-label comparison tables from the Aldinger, Sepp, and
# Science reference runs.  Run from the project root, for example:
#   Rscript outputs/merge_xenium_comparison_tables.R
# Or supply the directory that contains the three reference table directories:
#   Rscript outputs/merge_xenium_comparison_tables.R /path/to/outputs

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1) args[[1]] else "outputs"

references <- c("Aldinger", "Sepp", "Science")
table_dirs <- setNames(
  file.path(output_root, paste0("Xenium_", references, "ABT_Res1.5_Tables")),
  references
)
merged_dir <- file.path(output_root, "Xenium_Comp_ABT_Res1.5_Tables")

missing_dirs <- table_dirs[!dir.exists(table_dirs)]
if (length(missing_dirs) > 0) {
  stop("These reference table directories do not exist:\n", paste(missing_dirs, collapse = "\n"))
}
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)

file_suffix <- function(reference) {
  paste0("_", reference, "_majority_vs_weighted.csv")
}

comparison_files <- lapply(references, function(reference) {
  list.files(
    table_dirs[[reference]],
    pattern = paste0("_", reference, "_majority_vs_weighted\\.csv$"),
    full.names = TRUE
  )
})
names(comparison_files) <- references

sample_names <- sort(unique(unlist(Map(
  function(files, reference) sub(file_suffix(reference), "", basename(files), fixed = TRUE),
  comparison_files,
  references
))))

if (length(sample_names) == 0) {
  stop("No '*_majority_vs_weighted.csv' files were found in the reference table directories.")
}

read_reference_table <- function(sample_name, reference) {
  path <- file.path(table_dirs[[reference]], paste0(sample_name, file_suffix(reference)))
  if (!file.exists(path)) return(NULL)
  
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("seurat_clusters", "cluster_majority", "cluster_weighted")
  if (!all(required %in% names(tab))) {
    stop("Unexpected columns in ", path, ". Expected: ", paste(required, collapse = ", "))
  }
  
  tab <- tab[required]
  names(tab) <- c(
    "seurat_clusters",
    paste0(tolower(reference), "_cluster_majority"),
    paste0(tolower(reference), "_cluster_weighted")
  )
  tab
}

for (sample_name in sample_names) {
  tables <- lapply(references, function(reference) read_reference_table(sample_name, reference))
  names(tables) <- references
  
  absent <- names(tables)[vapply(tables, is.null, logical(1))]
  if (length(absent) > 0) {
    warning(sample_name, ": no comparison table found for ", paste(absent, collapse = ", "), 
            "; corresponding columns will be NA.")
    for (reference in absent) {
      empty_table <- data.frame(seurat_clusters = character(), check.names = FALSE)
      empty_table[[paste0(tolower(reference), "_cluster_majority")]] <- character()
      empty_table[[paste0(tolower(reference), "_cluster_weighted")]] <- character()
      tables[[reference]] <- empty_table
    }
  }
  
  # Full joins preserve a cluster even if a reference produced no label for it.
  merged <- Reduce(function(x, y) merge(x, y, by = "seurat_clusters", all = TRUE, sort = FALSE), tables)
  
  # Choose the most common annotation across all reference/voting methods,
  # ignoring Unknown and missing labels. In a tie, the first label in the
  # column order (Aldinger, Sepp, Science; majority before weighted) is used.
  label_columns <- setdiff(names(merged), "seurat_clusters")
  merged$consensus_label <- vapply(seq_len(nrow(merged)), function(i) {
    labels <- trimws(as.character(unlist(merged[i, label_columns, drop = FALSE])))
    labels <- labels[!is.na(labels) & nzchar(labels) & tolower(labels) != "unknown"]
    if (length(labels) == 0) return("Unknown")
    
    counts <- table(labels)
    winners <- names(counts)[counts == max(counts)]
    labels[match(winners, labels)][1]
  }, character(1))
  
  cluster_number <- suppressWarnings(as.numeric(as.character(merged$seurat_clusters)))
  merged <- merged[order(is.na(cluster_number), cluster_number, merged$seurat_clusters), , drop = FALSE]
  
  out_path <- file.path(merged_dir, paste0(sample_name, "_comparison_merged.csv"))
  write.csv(merged, out_path, row.names = FALSE, na = "")
  message("Wrote: ", out_path)
}

###Check if any samples had consensus label = Unknown
library(dplyr)

merged_dir <- "outputs/Xenium_Comp_ABT_Res1.5_Tables"

unknown_consensus <- list.files(
  merged_dir,
  pattern = "_comparison_merged\\.csv$",
  full.names = TRUE
) %>%
  lapply(function(file) {
    tab <- read.csv(file, stringsAsFactors = FALSE)
    
    tab %>%
      filter(consensus_label == "Unknown") %>%
      transmute(
        sample = sub("_comparison_merged\\.csv$", "", basename(file)),
        seurat_clusters
      )
  }) %>%
  bind_rows()

if (nrow(unknown_consensus) == 0) {
  message("No samples had clusters with consensus_label = Unknown.")
} else {
  print(unknown_consensus)
}