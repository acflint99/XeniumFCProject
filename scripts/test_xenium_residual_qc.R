#!/usr/bin/env Rscript
# Small synthetic checks only: no project objects or raw inputs are read.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_xenium_residual_qc.R")
source(here("scripts", "R", "xenium_residual_qc.R"))
expect_error <- function(expression) {
  if (!inherits(tryCatch(force(expression), error = identity), "error")) {
    stop("Expected validation error.")
  }
}
cells <- data.frame(
  cell_id = paste0("cell", 1:12),
  sample_id = rep(c("A", "B"), each = 6),
  cluster = c(rep("5", 2), rep("0", 4), rep("5", 4), rep("0", 2)),
  nCount_Xenium = c(rep(1, 2), rep(10, 4), rep(20, 4), rep(11, 2)),
  nucleus_area = c(NA_real_, 0, rep(5, 10)),
  nucleus_count = c(0, 0, rep(1, 10)), cell_area = 10
)
metrics <- c("nCount_Xenium", "nucleus_area")
original <- cells
result <- residual_qc_summarize(cells, metrics, min_cells_per_group = 2)
row <- result$balanced[result$balanced$cluster == "5" &
                        result$balanced$metric == "nCount_Xenium", ]
stopifnot(identical(cells, original), row$eligible_samples == 2,
          row$median_sample_difference == 0,
          row$samples_with_lower_target_median == 1,
          row$samples_with_higher_target_median == 1)
paired <- result$paired[result$paired$cluster == "5" &
                         result$paired$metric == "nCount_Xenium", ]
stopifnot(isTRUE(all.equal(paired$median_ratio_to_sample, c(.1, 1))),
          isTRUE(all.equal(sum(paired$fraction_of_cluster), 1)))
pooled <- result$pooled[result$pooled$cluster == "5" &
                         result$pooled$metric == "nCount_Xenium", ]
stopifnot(pooled$median == 20) # Cell pooling differs from equal sample weighting.
row <- result$pooled[result$pooled$cluster == "5" &
                      result$pooled$metric == "nucleus_area", ]
stopifnot(row$n_missing == 1, row$n_measured == 5)
small <- residual_qc_summarize(cells, metrics)
stopifnot(all(small$balanced$eligible_samples == 0),
          all(is.na(small$balanced$median_sample_difference)))
one_cluster <- cells; one_cluster$cluster <- "5"
single <- residual_qc_summarize(one_cluster, metrics, 2)
stopifnot(all(single$paired$background_cells == 0),
          all(is.na(single$paired$median_ratio)), !any(single$paired$eligible))
absent_area <- cells; absent_area$nucleus_area <- NA_real_
missing <- residual_qc_summarize(absent_area, metrics, 2)
stopifnot(all(is.na(missing$segmentation$fraction_nucleus_area_zero)))
duplicate <- cells; duplicate$cell_id[2] <- duplicate$cell_id[1]
expect_error(residual_qc_summarize(duplicate, metrics))
bad <- cells; bad$nCount_Xenium[1] <- Inf
expect_error(residual_qc_summarize(bad, metrics))
expect_error(residual_qc_numeric(c("1", "bad"), "counts"))
expect_error(residual_qc_numeric(c(1, NA_real_), "counts"))

path <- tempfile(fileext = ".csv")
raw <- data.frame(cell_id = c("b", "a"), cell_area = c(20.001, 10.001),
                  nucleus_area = c(0, NA_real_), nucleus_count = c(0, 1))
write.csv(raw, path, row.names = FALSE)
joined <- residual_qc_read_morphology(path, c("a", "b"), c(10, 20))
stopifnot(identical(joined$cell_id, c("a", "b")),
          is.na(joined$nucleus_area[1]), joined$nucleus_count[2] == 0)
expect_error(residual_qc_read_morphology(path, c("a", "missing"), c(10, 20)))
expect_error(residual_qc_read_morphology(path, c("a", "a"), c(10, 10)))
expect_error(residual_qc_read_morphology(path, c("a", "b"), c(9, 20)))
raw$cell_id[2] <- "b"; write.csv(raw, path, row.names = FALSE)
expect_error(residual_qc_read_morphology(path, "b", 20))
unlink(path)
cat("PASS: sample weighting, missing/zero values, empty backgrounds, IDs and morphology joins.\n")
