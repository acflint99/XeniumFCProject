#!/usr/bin/env Rscript
# Lightweight synthetic checks; no project RDS or biological analysis is loaded.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_vz_resolution_efficiency.R")
source(here("scripts", "R", "vz_resolution_efficiency.R"))
source(here("scripts", "R", "subcluster_review_markers.R"))
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
ids <- paste0("cell", 1:6)
columns <- paste0("VZ_compare_res.", c(.5, .6, .7, .8))
anchors <- data.frame(a = c("0", "0", "0", "1", "1", "1"),
                      b = c("0", "0", "1", "1", "2", "2"))
assignments <- data.frame(cell_id = ids, sample_id = rep(c("A", "B"), each = 3))
for (column in columns) assignments[[column]] <- anchors$b
assignments[[columns[1]]] <- anchors$a
key <- list(cell_ids = ids, sample_ids = assignments$sample_id, input_mtime = 1)
checkpoint <- list(key = key, assignments = assignments)
vz_comparison_validate_checkpoint(checkpoint, key, columns, anchors)
checkpoint_file <- tempfile(fileext = ".rds")
saveRDS(checkpoint, checkpoint_file, compress = FALSE)
vz_comparison_validate_checkpoint(readRDS(checkpoint_file), key, columns, anchors)
unlink(checkpoint_file)
stopifnot(identical(formals(write_focused_marker_dotplot)$pdf_device, quote(grDevices::cairo_pdf)))
changed_key <- key; changed_key$input_mtime <- 2
expect_error(vz_comparison_validate_checkpoint(checkpoint, changed_key, columns, anchors))
bad <- checkpoint; bad$assignments <- bad$assignments[6:1, ]
expect_error(vz_comparison_validate_checkpoint(bad, key, columns, anchors))
bad <- checkpoint; bad$assignments[[columns[1]]][1] <- "5"
expect_error(vz_comparison_validate_checkpoint(bad, key, columns, anchors))
bad <- checkpoint; bad$assignments[[columns[2]]][1] <- NA
expect_error(vz_comparison_validate_checkpoint(bad, key, columns, anchors))
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required for the slimming/marker equivalence test.")
future::plan(future::sequential)
set.seed(42)
counts <- matrix(rpois(40 * 120, 1), 40, 120,
  dimnames = list(paste0("gene", 1:40), paste0("cell", 1:120)))
for (i in 1:3) counts[(1:8) + (i - 1) * 8, (1:40) + (i - 1) * 40] <-
  counts[(1:8) + (i - 1) * 8, (1:40) + (i - 1) * 40] + 6
obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "Xenium")
obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj$cluster <- rep(as.character(0:2), each = 40)
SeuratObject::Idents(obj) <- "cluster"
coords <- matrix(rnorm(240), 120, 2,
  dimnames = list(colnames(obj), paste0("UMAP_", 1:2)))
obj[["umap_clean"]] <- SeuratObject::CreateDimReducObject(coords, key = "UMAP_", assay = "Xenium")
extra <- coords; colnames(extra) <- paste0("EXTRA_", 1:2)
obj[["extra"]] <- SeuratObject::CreateDimReducObject(extra, key = "EXTRA_", assay = "Xenium")
obj <- Seurat::FindNeighbors(obj, reduction = "umap_clean", dims = 1:2, k.param = 10, verbose = FALSE)
slim <- vz_comparison_slim(obj)
stopifnot(identical(colnames(obj), colnames(slim)), identical(obj[[]], slim[[]]),
  identical(SeuratObject::Embeddings(obj[["umap_clean"]]), SeuratObject::Embeddings(slim[["umap_clean"]])),
  identical(SeuratObject::LayerData(obj, assay = "Xenium", layer = "counts"),
            SeuratObject::LayerData(slim, assay = "Xenium", layer = "counts")),
  identical(SeuratObject::LayerData(obj, assay = "Xenium", layer = "data"),
            SeuratObject::LayerData(slim, assay = "Xenium", layer = "data")),
  length(slim@graphs) == 0, identical(names(slim@reductions), "umap_clean"),
  !"scale.data" %in% SeuratObject::Layers(slim[["Xenium"]]))
marker_test <- function(x) Seurat::FindAllMarkers(x, assay = "Xenium", only.pos = TRUE,
  min.pct = .25, logfc.threshold = .25, max.cells.per.ident = 30, random.seed = 42,
  test.use = "wilcox", verbose = FALSE)
before <- marker_test(obj); after <- marker_test(slim)
stopifnot(nrow(before) > 0, identical(before, after))
{
  pdf_path <- tempfile(fileext = ".pdf"); availability_path <- tempfile(fileext = ".csv")
  write_focused_marker_dotplot(slim, list("Synthetic" = c("gene1", "gene9", "gene17")),
    assay = "Xenium", group_by = "cluster", branch = "Synthetic", resolution_label = "test",
    availability_path = availability_path, tiff_path = NULL, pdf_path = pdf_path,
    pdf_device = grDevices::pdf)
  stopifnot(file.exists(pdf_path), file.info(pdf_path)$size > 0, file.exists(availability_path))
  unlink(c(pdf_path, availability_path))
  cat("PDF-only focused marker rendering passed.\n")
}
cat("Checkpoint validation, slimming preservation and exact marker equivalence passed.\n")
