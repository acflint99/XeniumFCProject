#!/usr/bin/env Rscript
# End-to-end synthetic fixture only: four cells, one sample. Expected 10-60 s.
# Never loads/writes production analysis objects. Temporary fixture retained for inspection.
origin <- normalizePath(".")
source(file.path(origin, "R", "atlas_qc.R"))
if (!requireNamespace("SeuratObject", quietly = TRUE) || !requireNamespace("yaml", quietly = TRUE)) stop("Use project dependencies; no installation attempted.")
fixture <- tempfile("atlas_qc_fixture_"); dir.create(fixture)
dir.create(file.path(fixture, "scripts", "R"), recursive = TRUE)
dir.create(file.path(fixture, "config")); file.create(file.path(fixture, ".here"))
for (f in c("xenium_atlas_qc.R", "xenium_atlas_plot_qc.R", "R/config.R", "R/atlas_qc.R")) {
  stopifnot(file.copy(file.path(origin, f), file.path(fixture, "scripts", f)))
}
config_dir <- if (file.exists(file.path(origin, "config", "config.yml"))) file.path(origin, "config") else file.path(dirname(origin), "config")
config <- yaml::read_yaml(file.path(config_dir, "config.yml"))
config$project$data_dir <- "data"; config$project$outputs_dir <- "outputs"
yaml::write_yaml(config, file.path(fixture, "config", "config.yml"))
samples <- atlas_read_csv(file.path(config_dir, "samples.csv"))[1, ]
samples$cell_stats_file <- "cerebellum_cells_stats.csv"
atlas_write(samples, file.path(fixture, "config", "samples.csv"))
for (f in c("atlas_sample_metadata.csv", "atlas_anatomical_coverage.csv", "atlas_segmentation_review.csv")) {
  x <- atlas_read_csv(file.path(origin, "config", f))
  if (nrow(x)) x <- x[x$sample_id == samples$sample_id, , drop = FALSE]
  atlas_write(x, file.path(fixture, "scripts", "config", f))
}
raw <- file.path(fixture, "data", samples$input_directory); dir.create(raw, recursive = TRUE)
m <- data.frame(nCount_Xenium = c(30, 10, 300, 40), nFeature_Xenium = c(3, 1, 3, 2),
  cell_area = c(30, 5, 40, 30), percent_control = c(0, 0, 0, 100*4/44), row.names = paste0("cell", 1:4))
object <- SeuratObject::CreateSeuratObject(counts = matrix(c(10,10,10,10,0,0,100,100,100,20,20,0), nrow = 3,
  dimnames = list(paste0("gene", 1:3), rownames(m))), assay = "Xenium")
object[["ControlProbe"]] <- SeuratObject::CreateAssayObject(counts = matrix(c(0,0,0,4), nrow = 1,
  dimnames = list("control", rownames(m))))
for (field in names(m)) object[[field]] <- m[[field]]
thresholds <- list(nCount_lower = 20, nCount_upper = 100, nFeature_lower = 2,
  nFeature_upper = 50, area_lower = 15, area_upper = 100, max_control_percent = 5)
write_rds <- function(x, folder, suffix) {
  p <- file.path(fixture, "outputs", "xenium", folder, "rds", paste0(samples$sample_id, suffix))
  dir.create(dirname(p), recursive = TRUE); saveRDS(x, p)
}
cropped <- object; cropped$percent_control <- NULL
write_rds(cropped, "preprocess/01_cropped", "_CB.rds")
qc <- subset(object, cells = "cell1"); qc@misc$QC_thresholds <- thresholds
write_rds(qc, "preprocess/02_qc", "_CB_QC.rds")
qc$comb_subcluster <- "example"; qc$PCW <- "PCW17"
write_rds(qc, "vz_rl/01_combined_labels", "_VZ_RL_QC_Subclusters.rds")
atlas_write(data.frame(`Cell ID` = rownames(m), cell_area = m$cell_area, check.names = FALSE), file.path(raw, "cerebellum_cells_stats.csv"))
atlas_write(data.frame(cell_id = rownames(m), cell_area = m$cell_area, nucleus_area = 10, nucleus_count = 1), file.path(raw, "cells.csv"))
atlas_write(data.frame(example_metric = "50%"), file.path(raw, "metrics_summary.csv"))
run <- function(driver, args, success = TRUE, pattern = NULL) {
  old <- setwd(fixture); on.exit(setwd(old))
  result <- suppressWarnings(system2(Sys.which("Rscript"), c("--vanilla", shQuote(file.path(fixture, "scripts", driver)), args), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status"); if (is.null(status)) status <- 0L
  if ((status == 0L) != success) stop(paste(result, collapse = "\n"))
  if (!is.null(pattern) && !any(grepl(pattern, result, fixed = TRUE))) stop("Expected diagnostic absent: ", pattern, "\n", paste(result, collapse = "\n"))
  invisible(result)
}
run("xenium_atlas_qc.R", c("--task", "1", "--dry-run"))
peer <- samples; peer$sample_id <- "fixture_peer"; peer$cell_stats_file <- "peer.csv"
write.csv(rbind(samples, peer), file.path(fixture, "config", "samples.csv"), row.names = FALSE)
atlas_write(data.frame(`Cell ID` = "cell1", check.names = FALSE), file.path(raw, "peer.csv"))
run("xenium_atlas_qc.R", c("--task", "1"), success = FALSE, pattern = "Overlapping crop membership")
write.csv(samples, file.path(fixture, "config", "samples.csv"), row.names = FALSE)
run("xenium_atlas_qc.R", c("--task", "1"))
run("xenium_atlas_qc.R", c("--task", "1"), success = FALSE)
run("xenium_atlas_qc.R", c("--compile"), success = FALSE)
run("xenium_atlas_qc.R", c("--compile", "--allow-incomplete"))
run("xenium_atlas_plot_qc.R", c("--pdf-only"))
out <- file.path(fixture, "outputs", "xenium", "atlas_qc")
stopifnot(all(atlas_read_csv(file.path(out, "retention.csv"))$n_cells == c(4, 1, 1)),
  file.exists(file.path(out, "plots", "atlas_cell_retention.pdf")),
  file.exists(file.path(out, "plots", "atlas_morphology_audit.png")))
run("xenium_atlas_plot_qc.R", c("--pdf-only"), success = FALSE)
cat("PASS: end-to-end four-cell extraction, raw-area/crop reconciliation, compile, PDF/PNG plots, strict incomplete refusal, and overwrite refusal.\nFixture: ", fixture, "\n", sep = "")
