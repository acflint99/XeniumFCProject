#!/usr/bin/env Rscript
# Tiny synthetic objects only; expected 20-90 seconds, no project RDS reads.
suppressPackageStartupMessages({library(here); library(Seurat); library(harmony); library(future)})
here::i_am("scripts/test_residual_post_qc_processing.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "residual_post_qc_processing.R"))
future::plan(future::sequential)
set.seed(42)
config <- load_pipeline_config()
# Isolated fixture settings; production configuration is never written.
config$regional_subsets$integration_dimensions_default <- 5L
config$regional_subsets$neighbor_k_post_qc <- 10L
counts <- matrix(rpois(120 * 120, 3), 120, 120,
                 dimnames = list(paste0("g", 1:120), paste0("cell", 1:120)))
x <- CreateSeuratObject(counts, assay = "Xenium")
x$orig.ident <- rep(c("S1", "S2"), each = 60)
x$sentinel_metadata <- seq_len(ncol(x))
x <- NormalizeData(x, verbose = FALSE)
x <- FindVariableFeatures(x, nfeatures = 100, verbose = FALSE)
x <- ScaleData(x, verbose = FALSE)
x <- RunPCA(x, npcs = 5, verbose = FALSE)
x <- RunHarmony(x, group.by.vars = "orig.ident", dims.use = 1:5, verbose = FALSE)
x <- RunUMAP(x, reduction = "harmony", dims = 1:5, n.neighbors = 10, reduction.name = "umap_harmony", verbose = FALSE)
x <- FindNeighbors(x, reduction = "harmony", dims = 1:5, k.param = 10, verbose = FALSE)
x <- FindClusters(x, resolution = c(.5, .8, 1), verbose = FALSE)
before <- x[[]]
reuse <- residual_postqc_representations(x, config, TRUE)
stopifnot(identical(reuse[[]][, names(before), drop = FALSE], before),
          identical(Embeddings(reuse, "harmony"), Embeddings(x, "harmony")),
          identical(Embeddings(reuse, "umap_harmony"), Embeddings(x, "umap_harmony")),
          identical(reuse[["Xenium_snn"]], x[["Xenium_snn"]]),
          identical(as.character(reuse$residual_postQC_Xenium_snn_res.0.8), as.character(x$Xenium_snn_res.0.8)))
y <- subset(x, cells = colnames(x)[-c(1:10, 61:70)])
retained <- y[[]]
z <- residual_postqc_representations(y, config, FALSE)
stopifnot(identical(colnames(z), colnames(y)),
          identical(z[[]][, names(retained), drop = FALSE], retained),
          all(is.finite(Embeddings(z, "harmony"))),
          identical(rownames(z[["Xenium_snn"]]), colnames(y)),
          all(paste0("residual_postQC_Xenium_snn_res.", c(.5, .8, 1)) %in% names(z[[]])))
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
expect_error(residual_postqc_representations(reuse, config, TRUE))
# Graph and clustering commands still describe the post-QC representations,
# even though original metadata columns have been preserved separately.
nn <- z@commands[grepl("^FindNeighbors", names(z@commands))]
cc <- z@commands[grepl("^FindClusters", names(z@commands))]
stopifnot(length(nn) == 1L, length(cc) == 1L,
          identical(nn[[1]]@params$reduction, "harmony"),
          identical(cc[[1]]@params$graph.name, "Xenium_snn"))
cat("PASS: no-removal reuse and after-removal refitting preserve cell IDs and all original metadata; post-QC partitions, embeddings, graph and saved commands align.\n")
