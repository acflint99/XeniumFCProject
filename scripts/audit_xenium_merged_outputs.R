#!/usr/bin/env Rscript
# Read-only inventory of production cohort-wide outputs through regional routing.
# Does not load RDS files, regenerate reports, or write outputs.
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "R", "config.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, c("--dry-run", "--list")))) {
  stop("Usage: Rscript scripts/audit_xenium_merged_outputs.R [--dry-run|--list]")
}
cfg <- load_pipeline_config()
sample_ids <- load_sample_manifest(cfg)$sample_id
production <- get_production_output_paths(cfg)
root <- here::here(cfg$project$outputs_dir)
annotation <- production$annotation_root
entries <- list()
add <- function(group, directory, files, expected_pages = NA_integer_) {
  for (file in files) entries[[length(entries) + 1L]] <<- data.frame(
    group = group, path = file.path(directory, file), expected_pages = expected_pages,
    stringsAsFactors = FALSE
  )
}
add("Consensus PDF", production$consensus_plot_dir, c(
  "Xenium_Consensus_Res5_Weighted2of3_All_Samples_Report.pdf",
  "Xenium_Consensus_Res5_Weighted2of3_All_Samples_KnownOnly_Report.pdf"
), length(sample_ids))
add("Aldinger PDF", file.path(production$label_transfer_dir, "aldinger", "plots"),
    "Xenium_Aldinger_Res5_All_Samples_Plot_Report.pdf", length(sample_ids))
add("Unknown QC PDF", file.path(annotation, "04_unknown_qc_review", "plots"),
    "Xenium_Unknown_QC_All_Samples.pdf")
add("Unknown QC", file.path(annotation, "04_unknown_qc_review", "tables"), c(
  "Xenium_Unknown_QC_Evidence.csv", "Xenium_Unknown_QC_Review_Template.csv",
  "Xenium_Unknown_QC_Audit_manifest.csv"
))
add("Broad calibration", file.path(root, "validation", "broad_reference_workflow", "03_calibration", "tables"), c(
  "Xenium_Broad_All_Cluster_Scores.csv", "Xenium_Broad_Reference_Threshold_Sensitivity.csv",
  "Xenium_Broad_Current_Numeric_Threshold_Consensus.csv", "Xenium_Broad_Current_Numeric_Threshold_Summary.csv",
  "Xenium_Broad_Reference_Threshold_Review_Template.csv", "Xenium_Broad_Calibration_Provenance.csv"
))
add("Broad consensus", file.path(annotation, "05_broad_reference_consensus", "tables"), c(
  "Xenium_Broad_Reference_Consensus.csv", "Xenium_Broad_Reference_Consensus_Summary.csv",
  "Xenium_Broad_Known_Consensus_Conflicts.csv", "Xenium_Broad_Reference_Consensus_Provenance.csv"
))
add("Reviewed QC inputs", file.path(annotation, "04_unknown_qc_reviewed", "tables"),
    "Xenium_Unknown_QC_routing_manifest.csv")
add("Reference routing", file.path(annotation, "04_unknown_reference_routing", "tables"), c(
  "Xenium_Unknown_Reference_Routes.csv", "Xenium_Unknown_Reviewed_Input_manifest.csv",
  "Xenium_Reference_Route_Validation.csv", "Xenium_Unknown_Reference_Route_Summary.csv",
  "Xenium_Unknown_Reference_Score_Marker_Evidence.csv"
))
add("Unrouted Unknowns", file.path(annotation, "05_unrouted_unknown_report", "tables"), c(
  "Xenium_Unrouted_Unknown_Populations.csv", "Xenium_Unrouted_Unknown_BySample.csv",
  "Xenium_Unrouted_Unknown_Summary.csv"
))
add("Consensus validation", file.path(root, "validation", "resolution5_through_consensus"), c(
  "reference_validation.csv", "sample_validation.csv", "consensus_label_counts.csv", "validation_summary.txt"
))
add("Consensus proportions", file.path(production$consensus_plot_dir, "proportions"), c(
  "XeniumConsensus_Res5_Weighted2of3_ClusterPropPlot.pdf",
  "XeniumConsensus_Res5_Weighted2of3_ClusterPropPlot.tif"
))
inventory <- do.call(rbind, entries)
if ("--list" %in% args) {
  write.table(inventory, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no")
}
pdfinfo <- Sys.which("pdfinfo")
inspect <- function(path, expected_pages) {
  result <- c(status = "MISSING", coverage = "", pages = "", modified = "")
  if (!file.exists(path)) return(result)
  info <- file.info(path)
  result[["modified"]] <- format(info$mtime, "%Y-%m-%d %H:%M:%S %Z")
  result[["status"]] <- if (is.na(info$size) || info$size <= 0) "EMPTY" else "PRESENT"
  if (result[["status"]] == "EMPTY") return(result)
  if (grepl("\\.csv$", path)) {
    data <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = identity)
    if (inherits(data, "error")) {
      result[["status"]] <- "UNREADABLE"
    } else {
      key <- intersect(c("sample_id", "sample", "sample_name"), names(data))
      if (length(key)) {
        ids <- unique(as.character(data[[key[[1]]]]))
        ids <- ids[!is.na(ids) & nzchar(ids)]
        result[["coverage"]] <- paste0(sum(sample_ids %in% ids), "/", length(sample_ids), " samples")
        if (length(setdiff(ids, sample_ids))) result[["status"]] <- "UNEXPECTED_SAMPLE_IDS"
      }
    }
  }
  if (grepl("\\.pdf$", path) && nzchar(pdfinfo)) {
    output <- suppressWarnings(system2(pdfinfo, shQuote(path), stdout = TRUE, stderr = TRUE))
    page_line <- grep("^Pages:[[:space:]]+[0-9]+", output, value = TRUE)
    if (length(page_line) == 1L) {
      pages <- as.integer(sub("^Pages:[[:space:]]+", "", page_line))
      result[["pages"]] <- as.character(pages)
      if (!is.na(expected_pages) && pages != expected_pages) result[["status"]] <- "PAGE_COUNT_MISMATCH"
    } else result[["status"]] <- "PDF_NOT_VERIFIED"
  }
  result
}
details <- t(vapply(seq_len(nrow(inventory)), function(i) {
  inspect(inventory$path[[i]], inventory$expected_pages[[i]])
}, c(status = "", coverage = "", pages = "", modified = "")))
result <- cbind(inventory, as.data.frame(details, stringsAsFactors = FALSE))
under_root <- startsWith(result$path, paste0(root, "/"))
result$path[under_root] <- substring(result$path[under_root], nchar(root) + 2L)
cat("Read-only audit: ", length(sample_ids), " manifest samples; no RDS objects loaded.\n", sep = "")
cat("Coverage is descriptive: conflict/unrouted tables may legitimately cover only a subset.\n")
cat("Page counts alone do not establish current content; this audit does not refresh outputs.\n")
write.table(result, row.names = FALSE, quote = FALSE, sep = "\t")
