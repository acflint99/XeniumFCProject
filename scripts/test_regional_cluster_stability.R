#!/usr/bin/env Rscript
# Lightweight synthetic checks; never read project analysis objects.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_regional_cluster_stability.R")
source(here("scripts", "R", "regional_cluster_stability.R"))
if (all(vapply(c("clustree", "ggraph", "ggplot2"), requireNamespace,
               logical(1), quietly = TRUE))) {
  regional_stability_check_plot()
} else message("Synthetic clustree rendering check skipped: clustree/ggraph dependencies unavailable.")
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
# Labels are arbitrary: exact recovery after permutation must equal one.
reference <- rep(c("a", "b"), each = 4)
exact <- regional_stability_jaccard(reference, rep(c("9", "3"), each = 4), c("a", "b", "absent"))
stopifnot(identical(exact$jaccard[1:2], c(1, 1)), is.na(exact$jaccard[3]),
          exact$retained_cells[3] == 0L)
# Both over-splitting and merging must lower scores, with correct denominators.
split <- regional_stability_jaccard(reference, rep(c("x", "y", "z", "w"), each = 2), c("a", "b"))
merge <- regional_stability_jaccard(reference, rep("x", 8), c("a", "b"))
stopifnot(all(split$jaccard == .5), all(merge$jaccard == .5))
expect_error(regional_stability_jaccard(reference, "x", c("a", "b")))
expect_error(regional_stability_jaccard(c("a", NA), c("x", "x"), "a"))
# Singleton biological samples must retain their actual index, never 1:index.
samples <- c(rep("A", 10), rep("B", 5), "singleton")
indices <- regional_stability_sample(samples, .8, 42)
stopifnot(identical(indices, regional_stability_sample(samples, .8, 42)),
  !anyDuplicated(indices), all(indices %in% seq_along(samples)),
  identical(as.integer(table(samples[indices])), c(8L, 4L, 1L)), 16L %in% indices)
expect_error(regional_stability_sample(samples, 1, 42))
expect_error(regional_stability_sample(c("A", ""), .8, 42))
# Exercise the same matrix -> SNN -> Louvain API and cell-order contract used
# by production review, on 100 artificial cells (no expression analysis).
if (requireNamespace("Seurat", quietly = TRUE)) {
  future::plan(future::sequential)
  set.seed(12)
  embedding <- matrix(rnorm(100 * 4), 100, 4,
                      dimnames = list(paste0("cell", seq_len(100)), paste0("PC_", 1:4)))
  graph <- Seurat::FindNeighbors(embedding, k.param = 10, verbose = FALSE)$snn
  clustered <- Seurat::FindClusters(graph, resolution = .5, algorithm = 1,
                                   random.seed = 42, verbose = FALSE)
  stopifnot(identical(rownames(clustered), rownames(embedding)), ncol(clustered) == 1L)
  # Confirm command records expose the settings needed for both reintegrated
  # and zero-removal production objects; do not guess graph parameters.
  counts <- matrix(rpois(500, 2), 5, 100,
    dimnames = list(paste0("gene", 1:5), rownames(embedding)))
  obj <- SeuratObject::CreateSeuratObject(counts, assay = "Xenium")
  colnames(embedding) <- paste0("HARMONY_", 1:4)
  obj[["harmony_clean"]] <- SeuratObject::CreateDimReducObject(embedding, key = "HARMONY_", assay = "Xenium")
  obj <- Seurat::FindNeighbors(obj, reduction = "harmony_clean", dims = 1:4, k.param = 10, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = c(.5, .8), verbose = FALSE)
  nn <- Filter(function(x) identical(x@params$reduction, "harmony_clean") &&
    "Xenium_snn" %in% x@params$graph.name, obj@commands[grepl("^FindNeighbors", names(obj@commands))])
  cc <- Filter(function(x) identical(x@params$graph.name, "Xenium_snn"),
    obj@commands[grepl("^FindClusters", names(obj@commands))])
  stopifnot(length(nn) == 1L, length(cc) == 1L,
    all(c("dims", "k.param", "nn.method", "n.trees", "annoy.metric", "nn.eps", "prune.SNN", "l2.norm") %in% names(nn[[1]]@params)),
    all(c("algorithm", "modularity.fxn", "n.start", "n.iter", "group.singletons") %in% names(cc[[1]]@params)))
  message("Synthetic Seurat graph and clustering API check passed.")
} else message("Seurat API check skipped: package unavailable.")
cat("Regional cluster stability synthetic checks passed.\n")
