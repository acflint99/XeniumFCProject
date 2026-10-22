#!/usr/bin/env Rscript
# Small synthetic checks only; no production RDS or clustering job is run.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_regional_resolution_comparison.R")
source(here("scripts", "R", "regional_resolution_comparison.R"))
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
ids <- paste0("cell", 1:6)
from <- c("0", "0", "0", "0", "1", "1")
to <- c("a", "a", "b", "b", "b", "b")
x <- regional_resolution_contingency(from, to, ids)
stopifnot(sum(x$cells) == 6L, nrow(x) == 3L,
  x$fraction_of_from[x$from_cluster == "0" & x$to_cluster == "b"] == .5,
  x$fraction_of_to[x$from_cluster == "0" & x$to_cluster == "b"] == .5,
  x$fraction_of_from[x$from_cluster == "1"] == 1,
  all(tapply(x$fraction_of_from, x$from_cluster, sum) == 1),
  all(tapply(x$fraction_of_to, x$to_cluster, sum) == 1))
expect_error(regional_resolution_contingency(from, to[-1], ids))
expect_error(regional_resolution_contingency(from, to, rep("duplicate", 6)))
expect_error(regional_resolution_contingency(from, c(to[-1], NA), ids))
if (requireNamespace("Seurat", quietly = TRUE)) {
  future::plan(future::sequential)
  set.seed(42)
  embedding <- matrix(rnorm(400), 100, 4, dimnames = list(paste0("cell", 1:100), paste0("PC_", 1:4)))
  graph <- Seurat::FindNeighbors(embedding, k.param = 10, verbose = FALSE)$snn
  original <- graph
  anchors <- Seurat::FindClusters(graph, resolution = c(.5, .8), algorithm = 1,
    random.seed = 0, verbose = FALSE)
  for (res in c(.6, .7)) {
    partition <- Seurat::FindClusters(graph, resolution = res, algorithm = 1,
      random.seed = 0, verbose = FALSE)
    repeated <- Seurat::FindClusters(graph, resolution = res, algorithm = 1,
      random.seed = 0, verbose = FALSE)
    stopifnot(identical(rownames(partition), rownames(embedding)), identical(partition, repeated),
      identical(graph, original))
    transition <- regional_resolution_contingency(as.character(anchors[[2]]),
      as.character(partition[[1]]), rownames(embedding))
    stopifnot(sum(transition$cells) == 100)
  }
  message("Saved-graph comparison API and deterministic seeds passed.")
} else message("Synthetic Seurat comparison check skipped: package unavailable.")
cat("Regional resolution comparison synthetic checks passed.\n")
