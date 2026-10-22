#!/usr/bin/env Rscript
# HPC extraction: 2-10 min/task, 28-140 min/41 tasks at concurrency 3.
# Compile: 1-5 min. Low confidence; RDS I/O dominates; no production objects written.
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
  scripts <- dirname(script)
  source(file.path(scripts, "R", "atlas_qc.R"), local = TRUE)
  source(file.path(scripts, "R", "config.R"), local = TRUE)
  config <- load_pipeline_config()
  samples <- load_sample_manifest(config)
  atlas_keys(samples, "sample_id", "Sample manifest")
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(); i <- 1L
  flags <- c("--list", "--dry-run", "--compile", "--allow-incomplete")
  valued <- c("--task", "--output-dir", "--metadata", "--anatomy", "--review")
  if ("--help" %in% args) {
    cat("--list | --dry-run | --task N | --compile [--allow-incomplete]\n",
      "[--output-dir PATH] [--metadata CSV] [--anatomy CSV] [--review CSV]\n",
      "--dry-run checks paths/small CSVs only; --allow-incomplete is compile-only.\n")
    return(invisible(NULL))
  }
  while (i <= length(args)) {
    key <- args[i]
    if (!key %in% c(flags, valued) || !is.null(opts[[key]])) stop("Unknown/repeated option: ", key)
    if (key %in% flags) opts[[key]] <- TRUE else {
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], "--")) stop("Missing option value: ", key)
      opts[[key]] <- args[i]
    }
    i <- i + 1L
  }
  ids <- samples$sample_id
  if (isTRUE(opts[["--list"]])) {
    if (length(opts) != 1L) stop("Use --list alone.")
    write.table(data.frame(task = seq_along(ids), sample_id = ids), stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  pick <- function(key, default) if (is.null(opts[[key]])) default else here::here(opts[[key]])
  out <- pick("--output-dir", here::here(config$project$outputs_dir, "xenium", "atlas_qc"))
  metadata_path <- pick("--metadata", file.path(scripts, "config", "atlas_sample_metadata.csv"))
  anatomy_path <- pick("--anatomy", file.path(scripts, "config", "atlas_anatomical_coverage.csv"))
  review_path <- pick("--review", file.path(scripts, "config", "atlas_segmentation_review.csv"))
  compile <- isTRUE(opts[["--compile"]]); dry <- isTRUE(opts[["--dry-run"]])
  if (compile && !is.null(opts[["--task"]])) stop("Choose compile or task, not both.")
  if (isTRUE(opts[["--allow-incomplete"]]) && !compile) stop("--allow-incomplete requires --compile.")
  if (!compile && is.null(opts[["--task"]]) && !dry) stop("Choose --task N, --compile, or --dry-run.")
  task <- if (is.null(opts[["--task"]])) seq_along(ids) else suppressWarnings(as.numeric(opts[["--task"]]))
  if (any(!is.finite(task)) || any(task < 1 | task > length(ids) | task != floor(task))) stop("Invalid task ID.")
  root <- here::here(config$project$outputs_dir, "xenium")
  stages <- c(cropped = "01_cropped", initial_qc = "02_qc")
  rds <- function(stage, suffix) file.path(root, "preprocess", stages[[stage]], "rds", paste0(ids, suffix))
  paths <- data.frame(sample_id = ids, cropped = rds("cropped", "_CB.rds"),
    initial_qc = rds("initial_qc", "_CB_QC.rds"),
    final = file.path(root, "vz_rl", "01_combined_labels", "rds", paste0(ids, "_VZ_RL_QC_Subclusters.rds")))
  bundles <- file.path(out, "samples", paste0(ids, "_atlas_qc"))
  suffixes <- c("retention", "quantiles", "failures", "populations", "controls", "provenance", "audit", "thresholds", "instrument_metrics")
  destinations <- file.path(out, paste0(c("coverage", "readiness", "evidence_inventory", "anatomy", "image_review", suffixes), ".csv"))
  expected <- unlist(lapply(bundles, function(b) paste0(b, "_", suffixes, ".csv")))
  required <- if (compile) c(metadata_path, anatomy_path, review_path, expected) else as.character(unlist(paths[task, -1, drop = FALSE]))
  if (compile && isTRUE(opts[["--allow-incomplete"]])) required <- c(metadata_path, anatomy_path, review_path)
  if (dry) {
    ok <- compact_dry_run("Atlas QC: paths only, no RDS loaded", inputs = required,
      outputs = if (compile) destinations else unlist(lapply(bundles[task], function(b) paste0(b, "_", suffixes, ".csv"))),
      checks = c(outputs_protected = !any(file.exists(if (compile) destinations else expected[unlist(lapply(task, function(t) ((t-1)*length(suffixes)+1):(t*length(suffixes))))]))))
    if (!ok) stop("Preflight failed. No analysis or writes performed.")
    return(invisible(NULL))
  }
  if (!compile) {
    dest <- paste0(bundles[task], "_", suffixes, ".csv")
    if (length(task) != 1L) stop("One task per process is required.")
    if (any(file.exists(dest))) stop("Extraction outputs exist; use a fresh --output-dir.")
    if (any(!file.exists(required))) stop("Missing RDS inputs: ", paste(required[!file.exists(required)], collapse = ", "))
    for (p in c("SeuratObject", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing project dependency: ", p)
    sample <- ids[task]
    raw <- resolve_sample_paths(samples[task, , drop = FALSE], config)$input_dir
    message("Atlas audit ", sample, ": extracting metadata sequentially; no RDS writes.")
    cropped <- atlas_metadata(paths$cropped[task], derive_qc = TRUE); gc()
    qc <- atlas_metadata(paths$initial_qc[task]); gc()
    flags <- atlas_qc_flags(cropped$metadata, qc$thresholds)
    keep <- rowSums(flags) == 0
    if (!setequal(rownames(qc$metadata), rownames(cropped$metadata)[keep])) stop("Stored QC thresholds do not reproduce retained cell IDs.")
    aligned <- cropped$metadata[rownames(qc$metadata), , drop = FALSE]
    for (metric in c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control")) {
      if (!metric %in% names(qc$metadata) || any(!is.finite(qc$metadata[[metric]])) ||
          any(abs(aligned[[metric]] - qc$metadata[[metric]]) > 1e-7)) stop("Initial QC metric differs from reconstructed crop: ", metric)
    }
    final <- atlas_metadata(paths$final[task]); gc()
    if (!all(rownames(final$metadata) %in% rownames(qc$metadata))) stop("Final cell IDs are not a subset of initial QC cells.")
    if (!"comb_subcluster" %in% names(final$metadata)) stop("Final object lacks comb_subcluster.")
    labels <- as.character(final$metadata$comb_subcluster)
    if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Missing final identities.")
    retention <- data.frame(sample_id = sample, stage = c("cropped", "initial_qc", "final"),
      n_cells = c(nrow(cropped$metadata), nrow(qc$metadata), nrow(final$metadata)))
    quantiles <- rbind(atlas_quantiles(cropped$metadata, sample, "cropped"),
      atlas_quantiles(qc$metadata, sample, "initial_qc"), atlas_quantiles(final$metadata, sample, "final"))
    patterns <- apply(flags, 1, function(z) if (any(z)) paste(names(flags)[z], collapse = "+") else "pass")
    failure_counts <- table(patterns)
    failures <- rbind(data.frame(sample_id = sample, kind = "individual_nonexclusive", reason = names(flags), n_cells = colSums(flags)),
      data.frame(sample_id = sample, kind = "exclusive_combination", reason = names(failure_counts), n_cells = as.integer(failure_counts)))
    tab <- table(labels)
    populations <- data.frame(sample_id = sample, label = names(tab), n_cells = as.integer(tab))
    control_names <- c("BlankCodeword", "ControlCodeword", "ControlProbe", "GenomicControl")
    controls <- data.frame(sample_id = sample, assay = control_names, imported = control_names %in% cropped$assays,
      interpretation = "Absence is unavailable evidence, not measured zero; present assay does not prove raw-import completeness.")
    audit <- data.frame(sample_id = sample, all_cropped_areas_equal_20 = all(cropped$metadata$cell_area == 20),
      missing_control_assays = paste(control_names[!controls$imported], collapse = ";"),
      object_PCW = if ("PCW" %in% names(final$metadata)) paste(unique(as.character(final$metadata$PCW)), collapse = ";") else NA_character_,
      note = "Constant area is a review flag, not proof dummy fallback occurred. Raw crop-list and area provenance require review.")
    thresholds <- data.frame(sample_id = sample, parameter = names(qc$thresholds), value = unlist(qc$thresholds))
    raw_paths <- resolve_sample_paths(samples[task, , drop = FALSE], config)
    if (!file.exists(raw_paths$cell_stats_path)) stop("Missing cerebellar crop list: ", raw_paths$cell_stats_path)
    crop_list <- read.csv(raw_paths$cell_stats_path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
    atlas_keys(crop_list, "Cell ID", "Cerebellar crop list")
    if (!setequal(crop_list$`Cell ID`, rownames(cropped$metadata))) stop("Crop-list IDs do not exactly match cropped object; investigate partial matches.")
    peer_paths <- character()
    peers <- which(samples$input_directory == samples$input_directory[task] & samples$sample_id != sample)
    for (peer in peers) {
      peer_path <- resolve_sample_paths(samples[peer, , drop = FALSE], config)$cell_stats_path
      if (!file.exists(peer_path)) stop("Missing same-acquisition crop list: ", peer_path)
      peer_list <- read.csv(peer_path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
      atlas_keys(peer_list, "Cell ID", "Same-acquisition crop list")
      if (length(intersect(crop_list$`Cell ID`, peer_list$`Cell ID`))) stop("Overlapping crop membership within raw acquisition: ", sample, " and ", samples$sample_id[peer])
      peer_paths <- c(peer_paths, peer_path)
    }
    raw_cells_path <- c(file.path(raw, "cells.csv.gz"), file.path(raw, "cells.csv"))
    raw_cells_path <- raw_cells_path[file.exists(raw_cells_path)]
    audit$raw_cells_available <- length(raw_cells_path) > 0L
    audit$raw_area_mismatch_cells <- NA_integer_
    audit$raw_nucleus_area_available_fraction <- NA_real_
    audit$raw_multiple_nuclei_fraction <- NA_real_
    audit$raw_no_nuclei_fraction <- NA_real_
    if (length(raw_cells_path)) {
      morphology <- atlas_read_csv(raw_cells_path[1])
      atlas_keys(morphology, "cell_id", "Raw cells table")
      index <- match(rownames(cropped$metadata), morphology$cell_id)
      if (anyNA(index)) stop("Cropped cells missing from raw morphology.")
      if (!"cell_area" %in% names(morphology)) stop("Raw cells table lacks cell_area.")
      m <- morphology[index, , drop = FALSE]
      if (!is.numeric(m$cell_area) || any(!is.finite(m$cell_area) | m$cell_area <= 0)) stop("Invalid raw cell areas.")
      audit$raw_area_mismatch_cells <- sum(abs(m$cell_area - cropped$metadata$cell_area) > .011)
      if ("nucleus_area" %in% names(m)) audit$raw_nucleus_area_available_fraction <- mean(is.finite(m$nucleus_area))
      if ("nucleus_count" %in% names(m)) {
        if (!is.numeric(m$nucleus_count) || any(!is.finite(m$nucleus_count) | m$nucleus_count < 0 | m$nucleus_count != floor(m$nucleus_count))) stop("Invalid raw nucleus counts.")
        audit$raw_multiple_nuclei_fraction <- mean(m$nucleus_count > 1)
        audit$raw_no_nuclei_fraction <- mean(m$nucleus_count == 0)
      }
      rm(morphology, m); gc()
    }
    provenance <- atlas_provenance(required); provenance$sample_id <- sample
    provenance <- rbind(provenance, transform(atlas_provenance(c(raw_paths$cell_stats_path, peer_paths, raw_cells_path,
      resolve_config_path(config$manifests$samples, config), attr(config, "config_path"))), sample_id = sample))
    metric_path <- file.path(raw, "metrics_summary.csv")
    instrument <- data.frame(sample_id = sample, input_directory = samples$input_directory[task],
      source_path = metric_path, scope = "input_directory_not_split_sample", metric = NA_character_, value = NA_character_,
      available = FALSE)
    if (file.exists(metric_path)) {
      raw_metrics <- atlas_read_csv(metric_path)
      if (!nrow(raw_metrics)) stop("Empty instrument metric file: ", metric_path)
      if (nrow(raw_metrics) == 1L) {
        instrument <- data.frame(sample_id = sample, input_directory = samples$input_directory[task],
          source_path = metric_path, scope = "input_directory_not_split_sample", metric = names(raw_metrics),
          value = as.character(unlist(raw_metrics[1, ], use.names = FALSE)), available = TRUE)
      } else if (ncol(raw_metrics) == 2L) {
        instrument <- data.frame(sample_id = sample, input_directory = samples$input_directory[task],
          source_path = metric_path, scope = "input_directory_not_split_sample", metric = as.character(raw_metrics[[1]]),
          value = as.character(raw_metrics[[2]]), available = TRUE)
      } else stop("Unsupported metrics_summary.csv schema; preserve and inspect: ", metric_path)
    }
    provenance <- rbind(provenance, transform(atlas_provenance(c(metric_path, file.path(raw, "experiment.xenium"))), sample_id = sample))
    result <- list(retention, quantiles, failures, populations, controls, provenance, audit, thresholds, instrument)
    # RDS MD5 is deliberately omitted: it would reread every large file.
    for (j in seq_along(dest)) atlas_write(result[[j]], dest[j])
    return(invisible(NULL))
  }
  if (any(file.exists(destinations))) stop("Compiled outputs exist; choose a fresh --output-dir.")
  meta <- atlas_check_metadata(atlas_read_csv(metadata_path), ids)
  anatomy <- atlas_read_csv(anatomy_path)
  atlas_keys(anatomy, c("sample_id", "compartment"), "Anatomy coverage")
  if (!all(c("status", "area_mm2", "area_method", "source") %in% names(anatomy))) stop("Invalid anatomy schema.")
  if (!all(anatomy$sample_id %in% ids) || !setequal(unique(anatomy$sample_id), ids)) stop("Anatomy table sample coverage mismatch.")
  if (anyNA(anatomy$status) || any(!anatomy$status %in% c("present", "absent", "uncertain", "not_assessed"))) stop("Invalid anatomy status.")
  areas <- suppressWarnings(as.numeric(anatomy$area_mm2))
  if (any(!is.na(anatomy$area_mm2) & (!is.finite(areas) | areas < 0))) stop("Invalid tissue areas.")
  if (any(!is.na(areas) & (is.na(anatomy$area_method) | is.na(anatomy$source)))) stop("Measured areas need method and source.")
  review <- atlas_read_csv(review_path)
  if (!all(c("sample_id", "field_id", "compartment", "selection_basis", "image_path", "reviewer", "review_date", "status", "findings") %in% names(review))) stop("Invalid image-review schema.")
  if (nrow(review)) {
    atlas_keys(review, c("sample_id", "field_id"), "Image review")
    if (any(!review$sample_id %in% ids) || anyNA(review$status) || any(!review$status %in% c("pending", "reviewed"))) stop("Invalid image-review sample/status.")
    complete_review <- review$status == "reviewed"
    for (column in c("compartment", "selection_basis", "image_path", "reviewer", "review_date", "findings")) {
      if (any(complete_review & (is.na(review[[column]]) | !nzchar(trimws(review[[column]]))))) stop("Reviewed fields need ", column)
    }
  }
  bundle_complete <- vapply(bundles, function(b) all(file.exists(paste0(b, "_", suffixes, ".csv"))), logical(1))
  verified <- !is.na(meta$donor_verified) & meta$donor_verified == "yes"
  anatomy_assessed <- vapply(ids, function(id) all(anatomy$status[anatomy$sample_id == id] %in% c("present", "absent")), logical(1))
  image_reviewed <- ids %in% review$sample_id[!is.na(review$status) & review$status == "reviewed"]
  readiness <- data.frame(sample_id = ids, donor_verified = verified, age_recorded = !is.na(meta$PCW),
    section_recorded = !is.na(meta$section_id), run_recorded = !is.na(meta$run_id),
    anatomy_assessed = anatomy_assessed, image_review_recorded = image_reviewed, extraction_complete = bundle_complete)
  for (field in c("section_plane", "slide_id", "panel_id", "xoa_version", "segmentation_method")) {
    readiness[[paste0(field, "_recorded")]] <- !is.na(meta[[field]])
  }
  readiness$report_complete <- apply(readiness[-1], 1, all)
  if (!all(readiness$report_complete) && !isTRUE(opts[["--allow-incomplete"]])) {
    stop("Incomplete atlas evidence. Fill templates/extract all samples, or explicitly use --allow-incomplete for a DRAFT gap report.")
  }
  manifest_coverage <- samples
  collisions <- intersect(setdiff(names(meta), "sample_id"), names(manifest_coverage))
  names(manifest_coverage)[match(collisions, names(manifest_coverage))] <- paste0("manifest_", collisions)
  coverage <- cbind(manifest_coverage, meta[setdiff(names(meta), "sample_id")])
  coverage$report_status <- if (all(readiness$report_complete)) "COMPLETE_REPORT_NOT_SCIENTIFIC_APPROVAL" else "DRAFT_INCOMPLETE"
  # Inventory existing small evidence without silently accepting absent outputs.
  inventory_paths <- c(file.path(root, "preprocess", "02_qc", "reports", paste0(ids, "_QC_thresholds.txt")),
    file.path(here::here(config$project$outputs_dir), "validation", "resolution5_through_consensus", c("sample_validation.csv", "consensus_label_counts.csv")))
  inventory <- data.frame(path = inventory_paths, exists = file.exists(inventory_paths), bytes = file.info(inventory_paths)$size)
  tables <- lapply(seq_along(suffixes), function(j) {
    if (!any(bundle_complete)) return(data.frame(sample_id = character()))
    present <- paste0(bundles[bundle_complete], "_", suffixes[j], ".csv")
    if (!length(present)) return(data.frame(sample_id = character()))
    values <- lapply(present, atlas_read_csv)
    for (k in seq_along(values)) if (anyNA(values[[k]]$sample_id) || !nrow(values[[k]]) || !all(values[[k]]$sample_id == ids[bundle_complete][k])) stop("Extraction sample ID mismatch: ", present[k])
    do.call(rbind, values)
  })
  audit_table <- tables[[match("audit", suffixes)]]
  if (any(bundle_complete)) atlas_validate_extracted(tables[[1]], tables[[4]], tables[[3]], ids[bundle_complete])
  if (nrow(audit_table)) {
    observed_age <- suppressWarnings(as.numeric(sub("^PCW\\s*", "", audit_table$object_PCW, ignore.case = TRUE)))
    expected_age <- meta$PCW[match(audit_table$sample_id, meta$sample_id)]
    if (any(!is.finite(observed_age)) || any(!is.na(expected_age) & observed_age != expected_age)) {
      stop("Final-object PCW missing, invalid, or inconsistent with atlas metadata; resolve before compiling.")
    }
  }
  provenance_table <- tables[[match("provenance", suffixes)]]
  if (nrow(provenance_table)) {
    recorded <- provenance_table$exists
    now <- file.info(provenance_table$path)
    if (any(recorded & (!file.exists(provenance_table$path) | is.na(now$size) | now$size != provenance_table$bytes |
      format(now$mtime, tz = "UTC", usetz = TRUE) != provenance_table$modified_utc))) stop("Extraction source missing/changed; regenerate in a fresh report directory.")
    hashed <- recorded & !is.na(provenance_table$md5)
    if (any(unname(tools::md5sum(provenance_table$path[hashed])) != provenance_table$md5[hashed])) stop("Extraction source checksum changed.")
  }
  readiness$raw_morphology_available <- FALSE
  readiness$raw_area_consistent <- FALSE
  readiness$instrument_metrics_recorded <- FALSE
  if (nrow(audit_table)) {
    matched <- match(ids, audit_table$sample_id)
    readiness$raw_morphology_available <- !is.na(matched) & !is.na(audit_table$raw_cells_available[matched]) & audit_table$raw_cells_available[matched]
    readiness$raw_area_consistent <- !is.na(matched) & !is.na(audit_table$raw_area_mismatch_cells[matched]) & audit_table$raw_area_mismatch_cells[matched] == 0
  }
  instrument_table <- tables[[match("instrument_metrics", suffixes)]]
  if (nrow(instrument_table)) readiness$instrument_metrics_recorded <- ids %in% instrument_table$sample_id[!is.na(instrument_table$available) & instrument_table$available]
  readiness$report_complete <- apply(readiness[setdiff(names(readiness), c("sample_id", "report_complete"))], 1, all)
  if (!all(readiness$report_complete) && !isTRUE(opts[["--allow-incomplete"]])) stop("Raw morphology/instrument evidence incomplete or area audit inconsistent; inspect extraction tables using a DRAFT report.")
  coverage$report_status <- if (all(readiness$report_complete)) "COMPLETE_REPORT_NOT_SCIENTIFIC_APPROVAL" else "DRAFT_INCOMPLETE"
  result <- c(list(coverage, readiness, inventory, anatomy, review), tables)
  for (j in seq_along(destinations)) atlas_write(result[[j]], destinations[j])
  atlas_write(atlas_provenance(c(metadata_path, anatomy_path, review_path, resolve_config_path(config$manifests$samples, config))), file.path(out, "report_sources.csv"))
  message("Atlas report: ", sum(readiness$report_complete), "/", length(ids), " entries complete. Missing evidence remains explicit.")
}
if (sys.nframe() == 0L) main()
