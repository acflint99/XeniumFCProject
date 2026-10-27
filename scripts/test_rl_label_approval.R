#!/usr/bin/env Rscript
# Synthetic approval fixtures only; expected <10 s, no production RDS loads.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_rl_label_approval.R")
scripts_dir <- here::here("scripts")
source(file.path(scripts_dir, "R", "rl_label_approval.R"))
fixture <- tempfile("rl_label_approval_test_")
dir.create(fixture)
input <- file.path(fixture, "source.rds")
review_path <- file.path(fixture, "review.csv")
approval_path <- file.path(fixture, "approved.csv")
saveRDS(list(synthetic = TRUE), input)
cell_ids <- paste0("S1_cell", 1:6)
cluster_ids <- c("0", "0", "0", "1", "1", "1")
column <- "Xenium_snn_res.0.8"
write.csv(data.frame(cells = 6, clusters = "0|1", cluster_column = column),
          review_path, row.names = FALSE)
template <- rl_label_approval_template(input, review_path, cell_ids, cluster_ids, column)
stopifnot(all(!template$approved), all(template$subcluster_label == ""),
          all(template$review_notes == ""), nrow(template) == 2L)
labels <- c("0" = "Differentiating GCs", "1" = "Mature UBCs")
check <- function(full = TRUE) validate_rl_label_approval(
  approval_path, input, review_path, column, labels = if (full) labels else NULL,
  cell_ids = if (full) cell_ids else NULL, cluster_ids = if (full) cluster_ids else NULL
)
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
expect_error(check()) # Missing decision.
write.csv(template, approval_path, row.names = FALSE)
expect_error(check()) # Evidence generation never approves labels.
approved <- template
approved$subcluster_label <- unname(labels[approved$cluster_id])
approved$approved <- TRUE
approved$review_notes <- "Synthetic fixture decision only"
write_approval <- function(x = approved) write.csv(x, approval_path, row.names = FALSE)
write_approval()
check(); check(full = FALSE) # Full validation and the path/small-table preflight.
write_approval(approved[2:1, ])
check() # Reviewer row order must not change cluster-to-label mapping.
for (field in c("input_size", "input_mtime", "cells")) {
  bad <- approved; bad[[field]] <- bad[[field]] + 1
  write_approval(bad); expect_error(check())
}
for (field in c("input_rds", "cluster_column", "source_partition_md5", "review_manifest_md5")) {
  bad <- approved; bad[[field]] <- "stale"
  write_approval(bad); expect_error(check())
}
bad <- approved; bad$approved[1] <- FALSE
write_approval(bad); expect_error(check())
bad <- approved; bad$approved[1] <- NA
write_approval(bad); expect_error(check())
bad <- approved; bad$review_notes[1] <- " "
write_approval(bad); expect_error(check())
bad <- approved; bad$subcluster_label[1] <- "Different label"
write_approval(bad); expect_error(check())
bad <- approved; bad$cluster_id[2] <- bad$cluster_id[1]
write_approval(bad); expect_error(check())
write_approval(approved[-1, , drop = FALSE]); expect_error(check())
write_approval()
old_clusters <- cluster_ids
cluster_ids <- c("1", "0", "0", "0", "1", "1")
expect_error(check()) # Same cell total and integer IDs; a different partition.
cluster_ids <- old_clusters
cell_ids[1] <- "changed_cell_id"
expect_error(check())
cell_ids[1] <- "S1_cell1"
expect_error(rl_partition_md5(rep("duplicate", 6), cluster_ids))
expect_error(rl_partition_md5(cell_ids, c(cluster_ids[-1], NA)))
old_review <- readLines(review_path)
writeLines(c(old_review, ""), review_path)
expect_error(check()) # Rewritten review evidence invalidates prior approval.
writeLines(old_review, review_path)
check()
# Exercise the actual driver's dry-run before its Seurat/RDS-loading branch.
dir.create(file.path(fixture, "scripts", "R"), recursive = TRUE)
dir.create(file.path(fixture, "config"))
file.create(file.path(fixture, ".here"))
for (name in c("config.R", "consensus_labels.R", "subcluster_review_markers.R",
               "dotplot_helpers.R", "xenium_replicate_marker_evidence.R", "rl_label_approval.R")) {
  stopifnot(file.copy(file.path(scripts_dir, "R", name), file.path(fixture, "scripts", "R", name)))
}
driver <- "xenium_rl_06_annotate_subclusters.R"
stopifnot(file.copy(file.path(scripts_dir, driver), file.path(fixture, "scripts", driver)),
          file.copy(file.path(scripts_dir, "config", "config.yml"),
                    file.path(fixture, "config", "config.yml")))
table_dir <- file.path(fixture, "outputs", "xenium", "rl", "06_subclusters", "tables")
plot_dir <- file.path(dirname(table_dir), "plots")
source_path <- file.path(fixture, "outputs", "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC.rds")
dir.create(dirname(source_path), recursive = TRUE)
# Invalid RDS content deliberately proves that inspection does not readRDS.
writeLines("Path-only fixture; not an RDS", source_path)
dir.create(table_dir, recursive = TRUE); dir.create(plot_dir)
replicate_dir <- file.path(table_dir, "replicate_marker_evidence")
dir.create(replicate_dir)
for (suffix in c("Evidence", "Positive_Marker_Candidates", "AntiMarker_Candidates",
                 "Marker_Sample_Coverage")) {
  name <- if (suffix == "Evidence") "ReplicateAware_Marker_Evidence" else suffix
  writeLines("fixture", file.path(replicate_dir, paste0("Xenium_RL_", name, ".csv")))
}
info <- file.info(source_path)
write.csv(data.frame(input_rds = normalizePath(source_path), input_size = info$size,
  input_mtime = as.numeric(info$mtime), cluster_column = column,
  inference_unit = "biological sample", cell_level_p_values_used_for_inference = FALSE),
  file.path(replicate_dir, "Xenium_RL_ReplicateAware_Marker_Provenance.csv"), row.names = FALSE)
for (name in c("RL_RawSubcluster_postQC_top10_Markers_Res0.8.csv",
               "RL_RawSubcluster_postQC_top50_Markers_Res0.8.csv",
               "RL_RawSubcluster_postQC_top50_AntiMarkers_Res0.8.csv",
               "RL_RawSubcluster_FocusedMarker_Availability_Res0.8.csv")) {
  writeLines("fixture", file.path(table_dir, name))
}
for (name in c("Xenium_RL_RawSubcluster_postQC_UMAP_Res0.8.tif",
               "Xenium_RL_RawSubcluster_FocusedMarker_DotPlot_Res0.8.tif",
               "Xenium_RL_RawSubcluster_FocusedMarker_DotPlot_Res0.8.pdf")) {
  writeLines("fixture", file.path(plot_dir, name))
}
review_fixture <- file.path(table_dir, "Xenium_RL_subcluster_review_manifest_Res0.8.csv")
write.csv(data.frame(cells = 6, clusters = "0|1", cluster_column = column,
  input_size = info$size, input_mtime = as.numeric(info$mtime)), review_fixture, row.names = FALSE)
approval_fixture <- rl_label_approval_template(source_path, review_fixture, cell_ids, cluster_ids, column)
write.csv(approval_fixture, file.path(table_dir, "Xenium_RL_Subcluster_Label_Approval_Template_Res0.8.csv"), row.names = FALSE)
approval_fixture$subcluster_label <- unname(labels[approval_fixture$cluster_id])
approval_fixture$approved <- TRUE; approval_fixture$review_notes <- "Synthetic only"
dir.create(file.path(fixture, "scripts", "reviews"))
approval_destination <- file.path(fixture, "scripts", "reviews", "RL_subcluster_labels_Res0.8.csv")
write.csv(approval_fixture, approval_destination, row.names = FALSE)
run_inspection <- function(mode, success) {
  previous <- getwd(); setwd(fixture); on.exit(setwd(previous))
  log <- tempfile()
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--max-ppsize=500000", file.path("scripts", driver), mode, "--dry-run"), stdout = log, stderr = log)
  output <- paste(readLines(log), collapse = "\n")
  unlink(log)
  if ((status == 0L) != success || !grepl(if (success) "DRY-RUN PASS" else "DRY-RUN FAIL", output)) stop(output)
}
run_inspection("--review-only", TRUE)
run_inspection("--apply-reviewed-labels", TRUE)
approval_fixture$approved[1] <- FALSE
write.csv(approval_fixture, approval_destination, row.names = FALSE)
run_inspection("--apply-reviewed-labels", FALSE)
annotated <- file.path(dirname(table_dir), "rds", "Xenium_RL_subclusters.rds")
dir.create(dirname(annotated)); writeLines("Not an RDS", annotated)
run_inspection("--plots-only", TRUE) # Existing plot refresh bypasses label application.
unlink(fixture, recursive = TRUE)
cat("PASS: blank/unapproved decisions, source/review binding, exact partitions, ",
    "cell IDs, cluster coverage, configured labels and actual driver preflights; ",
    "no production data loaded.\n", sep = "")
