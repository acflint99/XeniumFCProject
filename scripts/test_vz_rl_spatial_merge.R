#!/usr/bin/env Rscript
# Tiny synthetic FOVs only; no project inputs or outputs are read or written.
suppressPackageStartupMessages({
  library(here)
  library(SeuratObject)
})
here::i_am("scripts/test_vz_rl_spatial_merge.R")
source(here("scripts", "R", "vz_rl_batch_merge.R"))

expect_error <- function(expression) {
  result <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(result, "error"))
}
fixture <- function(i) {
  sample_id <- paste0("sample", i)
  cells <- paste0(sample_id, "_cell", 1:3)
  counts <- matrix(seq_len(12L) + i, nrow = 4L,
                   dimnames = list(paste0("gene", 1:4), cells))
  object <- CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE),
                               assay = "Xenium", project = "Xenium")
  object$sample_id <- sample_id
  object$approved_label <- rep(c("VZ", "RL", "Residual"), length.out = 3L)
  object[["fov"]] <- CreateFOV(
    coords = data.frame(x = 1:3 + i, y = 4:6 + i, cell = cells),
    type = "centroids", assay = "Xenium"
  )
  object
}

raw <- lapply(seq_len(41L), fixture)
groups <- split(seq_along(raw), ceiling(seq_along(raw) / 6L))
old_merge <- function(objects) {
  if (length(objects) == 1L) return(objects[[1L]])
  merge(objects[[1L]], y = objects[-1L], collapse = FALSE,
        merge.data = TRUE, merge.dr = FALSE, project = "Xenium_VZRL_Spatial")
}
old_batches <- lapply(groups, function(index) old_merge(raw[index]))
old_result <- suppressWarnings(old_merge(old_batches))
cat("Unqualified-name reproduction: expected 41 FOVs; observed",
    length(Images(old_result)), "\n")
stopifnot(length(Images(old_result)) < length(raw))
expect_error(merge_vz_rl_objects(raw[1:2]))

qualified <- lapply(seq_along(raw), function(i) {
  before <- raw[[i]]
  after <- namespace_vz_rl_images(before, paste0("sample", i))
  stopifnot(identical(before[["fov"]], after[[Images(after)[[1L]]]]),
            identical(before[[]], after[[]]),
            identical(colnames(before), colnames(after)),
            identical(Images(before), "fov"))
  validObject(after)
  after
})
batched <- merge_vz_rl_objects(lapply(groups, function(index) {
  merge_vz_rl_objects(qualified[index])
}))
pairwise <- Reduce(function(x, y) merge_vz_rl_objects(list(x, y)), qualified)
expected_cells <- unlist(lapply(qualified, colnames), use.names = FALSE)
expected_images <- unlist(lapply(qualified, Images), use.names = FALSE)
stopifnot(length(Images(batched)) == 41L,
          identical(Images(batched), expected_images),
          identical(Images(pairwise), expected_images),
          identical(colnames(batched), expected_cells),
          identical(colnames(pairwise), expected_cells),
          identical(batched$sample_id, pairwise$sample_id),
          identical(batched$approved_label, pairwise$approved_label))
for (i in seq_along(qualified)) {
  image <- expected_images[[i]]
  stopifnot(identical(Cells(batched[[image]]), Cells(qualified[[i]][[image]])),
            identical(GetTissueCoordinates(batched[[image]]),
                      GetTissueCoordinates(qualified[[i]][[image]])))
}
batched_counts <- LayerData(JoinLayers(batched[["Xenium"]]), layer = "counts")
pairwise_counts <- LayerData(JoinLayers(pairwise[["Xenium"]]), layer = "counts")
expected_counts <- do.call(cbind, lapply(qualified, function(object) {
  as.matrix(LayerData(object, assay = "Xenium", layer = "counts"))
}))
stopifnot(identical(as.matrix(batched_counts), expected_counts),
          identical(as.matrix(pairwise_counts), expected_counts))
expect_error(namespace_vz_rl_images(raw[[1L]], ""))
cat("PASS: 41 samples / 7 batches retain every FOV, coordinate, cell ID, count, and tested metadata column.\n")
