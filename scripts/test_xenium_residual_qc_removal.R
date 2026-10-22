#!/usr/bin/env Rscript
# Tiny isolated fixtures, never project analysis objects. Expected runtime <1 minute.
suppressPackageStartupMessages({library(here); library(Seurat); library(yaml)})
here::i_am("scripts/test_xenium_residual_qc_removal.R")
original <- here("scripts")
source(file.path(original, "R", "config.R"))
source(file.path(original, "R", "consensus_labels.R"))
source(file.path(original, "R", "xenium_residual_qc_removal.R"))
expect_error <- function(expression) {
  if (!inherits(tryCatch(force(expression), error = identity), "error")) stop("Expected error.")
}
fixture <- tempfile("residual_removal_test_")
dir.create(file.path(fixture, "scripts", "R"), recursive = TRUE)
dir.create(file.path(fixture, "config"))
for (f in c("config.R", "consensus_labels.R", "xenium_residual_qc_removal.R", "xenium_replicate_marker_evidence.R")) {
  stopifnot(file.copy(file.path(original, "R", f), file.path(fixture, "scripts", "R", f)))
}
drivers <- c("xenium_residual_04b_apply_qc.R", "xenium_residual_04_review_subclusters.R",
             "xenium_residual_05_apply_reviewed_labels.R", "xenium_residual_06_map_to_samples.R")
for (f in drivers) stopifnot(file.copy(file.path(original, f), file.path(fixture, "scripts", f)))
config <- yaml::read_yaml(file.path(original, "config", "config.yml"))
yaml::write_yaml(config, file.path(fixture, "config", "config.yml"))
write.csv(data.frame(sample_id = c("S1", "S2", "EMPTY"), input_directory = c("S1", "S2", "EMPTY"),
                     input_layout = "single", cell_stats_file = "", notes = ""),
           file.path(fixture, "config", "samples.csv"), row.names = FALSE)
setwd(fixture)
here::i_am("scripts/xenium_residual_04b_apply_qc.R")
config <- load_pipeline_config()
p <- residual_qc_removal_paths(config)
raw <- paste0("cell", 1:30)
cell_ids <- c(paste0("S1_", raw), paste0("S2_", raw))
cluster <- rep(rep(c("0", "5", "9"), each = 10), 2)
counts <- matrix(2, nrow = 4, ncol = 60, dimnames = list(paste0("Gene", 1:4), cell_ids))
counts[1, cluster == "0"] <- 20; counts[4, cluster == "9"] <- 20
object <- CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "Xenium")
object <- NormalizeData(object, verbose = FALSE)
object$orig.ident <- rep(c("S1", "S2"), each = 30)
object$consensus_label <- "Unknown-01"
object$residual_reference_route <- "route_unresolved"
object$residual_subset_reason <- "retained_high_quality_other_or_unresolved"
object$Xenium_snn_res.0.8 <- cluster
object$unchanged_metadata <- seq_len(60)
coordinates <- matrix(seq_len(120), ncol = 2, dimnames = list(cell_ids, c("UMAP_1", "UMAP_2")))
object[["umap_harmony"]] <- CreateDimReducObject(embeddings = coordinates, key = "UMAP_", assay = "Xenium")
dir.create(dirname(p$source), recursive = TRUE)
dir.create(dirname(p$qc_cells), recursive = TRUE)
saveRDS(object, p$source)
info <- file.info(p$source)
write.csv(data.frame(cell_id = cell_ids, sample_id = object$orig.ident, cluster = cluster), p$qc_cells, row.names = FALSE)
write.csv(data.frame(input_path = normalizePath(p$source), input_size = info$size,
                     input_mtime = as.numeric(info$mtime), cluster_column = "Xenium_snn_res.0.8"),
           p$qc_provenance, row.names = FALSE)
stopifnot(!load_residual_qc_state(config)$active)
d <- data.frame(cluster = "5", cluster_column = "Xenium_snn_res.0.8", expected_input_cells = 60,
                 expected_removed_cells = 20, source_size = info$size, source_mtime = as.numeric(info$mtime),
                 review_notes = "Synthetic approved test only", approved = TRUE)
write.csv(d, p$decision, row.names = FALSE)
expect_error(load_residual_qc_state(config)) # Approved but not applied: fail closed.
run <- function(driver, args = character(), success = TRUE) {
  log <- tempfile()
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--max-ppsize=500000", file.path("scripts", driver), args), stdout = log, stderr = log)
  if ((status == 0L) != success) stop(paste(readLines(log), collapse = "\n"))
}
run(drivers[[1]], "--dry-run")
hash <- unname(tools::md5sum(p$source))
run(drivers[[1]])
state <- load_residual_qc_state(config)
filtered <- readRDS(p$filtered)
keep <- cluster != "5"
stopifnot(identical(unname(tools::md5sum(p$source)), hash),
          ncol(filtered) == 40, identical(colnames(filtered), cell_ids[keep]),
          identical(filtered$unchanged_metadata, object$unchanged_metadata[keep]),
          identical(Embeddings(filtered, "umap_harmony"), coordinates[keep, , drop = FALSE]),
          identical(as.character(filtered$Xenium_snn_res.0.8), cluster[keep]),
          identical(as.integer(state$manifest$removed_cells), c(10L, 10L, 0L)))
validate_residual_qc_object(filtered, state)
expect_error(validate_residual_qc_object(object, state)) # Cannot reintroduce removed cells.
exclusions <- residual_qc_exclusions_for_sample(raw, rep("Unknown-01", 30), "S1", state)
stopifnot(identical(exclusions, cluster[1:30] == "5"),
          !any(residual_qc_exclusions_for_sample(raw, rep("Unknown-01", 30), "OTHER", state)))
expect_error(residual_qc_exclusions_for_sample(raw[-11], rep("Unknown-01", 29), "S1", state))
expect_error(residual_qc_exclusions_for_sample(raw, rep("Known", 30), "S1", state))
run(drivers[[1]], success = FALSE) # Protected outputs.
run(drivers[[2]], "--dry-run")
run(drivers[[2]])
review_path <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels, config)
review <- read.csv(file.path(p$markers, "Xenium_Residual_Subcluster_Label_Review_Template.csv"),
                    stringsAsFactors = FALSE)
stopifnot(setequal(review$cluster, c("0", "9")))
review$final_label <- paste0("synthetic_", review$cluster)
for (name in c("annotation_confidence", "positive_marker_support", "anti_marker_support",
               "cross_sample_support", "review_notes")) review[[name]] <- "synthetic test only"
review$approved <- TRUE
write.csv(review, review_path, row.names = FALSE)
run(drivers[[3]], "--dry-run")
run(drivers[[3]])
master <- readRDS(file.path(here(config$project$outputs_dir), "xenium", "residual", "05_subclusters", "rds", "Xenium_Residual_Subclusters.rds"))
validate_residual_qc_object(master, state)
stopifnot(identical(colnames(master), colnames(filtered)))
# Exercise Stage 6 against tiny reviewed whole-tissue fixtures. Only the upstream
# input-manifest reader is stubbed; the production mapping driver runs unchanged.
writeLines("# Synthetic upstream fixture; no project routing data read.",
           file.path(fixture, "scripts", "R", "xenium_unknown_qc_review.R"))
writeLines(c("load_unknown_reviewed_input_manifest <- function(config) {",
             "read.csv(here::here('fixture_inputs.csv'), stringsAsFactors=FALSE)", "}"),
           file.path(fixture, "scripts", "R", "xenium_unknown_reference_routes.R"))
status_dir <- file.path(here(config$project$outputs_dir), "xenium", "residual", "01_subsets", "tables")
dir.create(status_dir, recursive = TRUE)
input_rows <- list()
for (id in c("S1", "S2", "EMPTY")) {
  n <- if (id == "EMPTY") 2L else 30L
  names <- if (id == "EMPTY") c("other1", "other2") else raw
  base_counts <- Matrix::Matrix(matrix(2, nrow = 4, ncol = n,
    dimnames = list(paste0("Gene", 1:4), names)), sparse = TRUE)
  base <- CreateSeuratObject(base_counts, assay = "Xenium")
  base$consensus_label <- if (id == "EMPTY") "RL" else "Unknown-01"
  base$PCW <- "PCW12"
  path <- file.path(fixture, paste0(id, "_Input.rds"))
  saveRDS(base, path)
  fi <- file.info(path)
  input_rows[[id]] <- data.frame(sample_id = id, mapping_input_path = normalizePath(path))
  write.csv(data.frame(sample_id = id, input_rds = normalizePath(path), input_size = fi$size,
                       input_mtime = as.numeric(fi$mtime), residual_cells = if (id == "EMPTY") 0L else n),
             file.path(status_dir, paste0(id, "_ResidualSubset_Status.csv")), row.names = FALSE)
}
write.csv(do.call(rbind, input_rows), "fixture_inputs.csv", row.names = FALSE)
run(drivers[[4]], "--dry-run")
run(drivers[[4]])
mapped_root <- file.path(here(config$project$outputs_dir), "xenium", "residual", "06_mapped")
mapping <- read.csv(file.path(mapped_root, "tables", "Xenium_Residual_Mapping_Manifest.csv"))
stopifnot(identical(mapping$sample_id, c("S1", "S2", "EMPTY")),
          identical(mapping$mapped_residual_cells, c(20L, 20L, 0L)),
          identical(mapping$removed_residual_qc_cells, c(10L, 10L, 0L)))
for (id in c("S1", "S2", "EMPTY")) {
  mapped <- readRDS(file.path(mapped_root, "rds", paste0(id, "_Residual_Subclusters.rds")))
  excluded <- residual_qc_exclusions_for_sample(colnames(mapped), mapped$consensus_label, id, state)
  stopifnot(identical(unname(mapped$residual_qc_status == "removed"), excluded),
            all(is.na(mapped$residual_subcluster[excluded])),
            identical(mapped@misc$residual_qc_decision_md5, state$decision_md5))
}
d$expected_removed_cells <- 19; write.csv(d, p$decision, row.names = FALSE)
expect_error(load_residual_qc_state(config)) # Changed decision cannot silently reuse outputs.
cat("PASS: approved removal, snapshot binding, exact IDs/metadata/embedding, empty samples, protected outputs, post-QC markers/labels, whole-tissue mapping and exclusion joins.\n")
