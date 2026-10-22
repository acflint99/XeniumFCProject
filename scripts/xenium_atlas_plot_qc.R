#!/usr/bin/env Rscript
# Small compiled tables only: expected 1-5 min, no RDS loads.
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
  scripts <- dirname(script)
  source(file.path(scripts, "R", "atlas_qc.R"), local = TRUE)
  source(file.path(scripts, "R", "config.R"), local = TRUE)
  config <- load_pipeline_config()
  args <- commandArgs(trailingOnly = TRUE)
  out <- here::here(config$project$outputs_dir, "xenium", "atlas_qc")
  pdf_only <- FALSE; dry <- FALSE; i <- 1L
  while (i <= length(args)) {
    if (args[i] == "--output-dir") {
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], "--")) stop("Missing output directory.")
      out <- here::here(args[i])
    } else if (args[i] == "--pdf-only") pdf_only <- TRUE else if (args[i] == "--dry-run") dry <- TRUE else {
      stop("Usage: xenium_atlas_plot_qc.R [--output-dir PATH] [--pdf-only] [--dry-run]")
    }
    i <- i + 1L
  }
  files <- file.path(out, paste0(c("coverage", "readiness", "anatomy", "retention", "quantiles", "failures", "populations", "controls", "audit"), ".csv"))
  plot_dir <- file.path(out, "plots")
  # A new plot directory is required; prevents partial or concurrent replacement.
  if (dry) {
    ok <- compact_dry_run("Atlas plots: small tables only; no RDS", inputs = files,
      checks = c(new_plot_directory = !dir.exists(plot_dir), Cairo_available = pdf_only || requireNamespace("Cairo", quietly = TRUE)))
    if (!ok) stop("Plot preflight failed.")
    return(invisible(NULL))
  }
  if (any(!file.exists(files))) stop("Compile atlas tables first.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing ggplot2 in project renv.")
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))
  if (!pdf_only && !requireNamespace("Cairo", quietly = TRUE)) stop("Cairo is required for production TIFFs; --pdf-only supports lightweight local review.")
  if (!dir.create(plot_dir, showWarnings = FALSE)) stop("Plot directory exists or cannot be created; choose a fresh report output directory.")
  coverage <- atlas_read_csv(files[1]); ready <- atlas_read_csv(files[2]); anatomy <- atlas_read_csv(files[3])
  atlas_keys(coverage, "sample_id", "Compiled coverage")
  ready <- atlas_align(ready, coverage$sample_id, "Readiness")
  draft <- any(coverage$report_status == "DRAFT_INCOMPLETE")
  caption <- if (draft) "DRAFT: incomplete evidence; section entries are not independent donors" else "QC evidence summary; completeness does not certify scientific validity"
  plots <- list(); manifest <- data.frame(panel = character(), status = character(), reason = character())
  record <- function(name, status, reason = "") {
    manifest <<- rbind(manifest, data.frame(panel = name, status = status, reason = reason))
  }
  add <- function(name, plot, width = 8, height = 5) {
    plot <- plot +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(caption = caption)
    path <- file.path(plot_dir, paste0(name, ".pdf"))
    pdf_device <- if (pdf_only) function(filename, width, height, ...) {
      grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE)
    } else grDevices::cairo_pdf
    ggplot2::ggsave(path, plot, device = pdf_device, width = width, height = height, limitsize = FALSE)
    if (!file.exists(path) || file.info(path)$size <= 0) stop("PDF export failed: ", path)
    # PNG previews allow inspection without loading full-resolution TIFFs.
    ggplot2::ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 130, limitsize = FALSE)
    if (!pdf_only) {
      Cairo::CairoTIFF(file.path(plot_dir, paste0(name, ".tiff")), width = width, height = height, units = "in", res = 600)
      print(plot); grDevices::dev.off()
    }
    plots[[name]] <<- plot
    record(name, "generated")
  }
  y_order <- rev(coverage$sample_id)
  checks <- setdiff(names(ready), c("sample_id", "report_complete"))
  missingness <- do.call(rbind, lapply(checks, function(k) data.frame(sample_id = ready$sample_id,
    evidence = k, recorded = ready[[k]])))
  missingness$sample_id <- factor(missingness$sample_id, levels = y_order)
  add("atlas_evidence_completeness", ggplot2::ggplot(missingness, ggplot2::aes(evidence, sample_id, fill = recorded)) +
    ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_manual(values = c(`FALSE` = "#E6E6E6", `TRUE` = "#247BA0")) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Recorded", title = "Atlas evidence coverage") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)), 12, max(6, nrow(coverage)*.21))
  anatomy$sample_id <- factor(anatomy$sample_id, levels = y_order)
  add("atlas_anatomical_coverage", ggplot2::ggplot(anatomy, ggplot2::aes(compartment, sample_id, fill = status)) +
    ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_manual(values = c(present = "#247BA0", absent = "#E9C46A", uncertain = "#B79CED", not_assessed = "#E6E6E6")) +
    ggplot2::labs(x = NULL, y = NULL, title = "Reviewed anatomical coverage") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)), 10, max(6, nrow(coverage)*.21))
  # Input directory is a technical grouping only, never relabeled as donor or run.
  inputs <- as.data.frame(table(coverage$input_directory), stringsAsFactors = FALSE)
  names(inputs) <- c("input_directory", "sample_entries")
  add("atlas_input_directory_coverage", ggplot2::ggplot(inputs, ggplot2::aes(sample_entries, reorder(input_directory, sample_entries))) +
    ggplot2::geom_col(fill = "#247BA0") + ggplot2::labs(x = "Manifest entries per raw-input directory", y = NULL,
      title = "Technical input coverage", subtitle = "Input directories are not verified physical slides or biological donors"), 9, max(5, nrow(inputs)*.22))
  verified <- !is.na(coverage$donor_verified) & coverage$donor_verified == "yes" & !is.na(coverage$donor_id)
  known_age <- verified & !is.na(coverage$PCW)
  if ("section_order_lateral_to_medial" %in% names(coverage)) {
    sections <- coverage[!is.na(coverage$section_order_lateral_to_medial) & !is.na(coverage$donor_id), ]
    if (nrow(sections)) add("atlas_section_position_coverage", ggplot2::ggplot(sections,
      ggplot2::aes(section_order_lateral_to_medial, donor_id, shape = donor_verified)) +
      ggplot2::geom_point(size = 2, color = "#247BA0") + ggplot2::labs(x = "Section order: lateral toward medial", y = "Cerebellum group",
        title = "Within-cerebellum section coverage", subtitle = "Ordinal position, not distance; unnumbered sections omitted", shape = "Grouping status"), 10, 5)
  }
  working <- !is.na(coverage$donor_id) & !is.na(coverage$PCW)
  if (any(working)) {
    working_age <- unique(coverage[working, c("donor_id", "PCW")])
    working_counts <- aggregate(rep(1L, nrow(working_age)), list(PCW = working_age$PCW), sum)
    names(working_counts)[2] <- "working_donors"
    add("atlas_working_donors_by_age", ggplot2::ggplot(working_counts, ggplot2::aes(PCW, working_donors)) +
      ggplot2::geom_col(width = .4, fill = "#B79CED") + ggplot2::labs(x = "Age (PCW)", y = "Working cerebellum groups",
        title = "Provisional developmental coverage", subtitle = "Includes assumed FB/GZFB merges; PCW from dated metadata snapshot until HPC reconciliation"))
  }
  if (any(known_age)) {
    donor_age <- unique(coverage[known_age, c("donor_id", "PCW")])
    age_counts <- aggregate(rep(1L, nrow(donor_age)), list(PCW = donor_age$PCW), sum)
    names(age_counts)[2] <- "independent_donors"
    add("atlas_donors_by_age", ggplot2::ggplot(age_counts, ggplot2::aes(PCW, independent_donors)) +
      ggplot2::geom_col(width = .4, fill = "#247BA0") + ggplot2::scale_y_continuous(breaks = scales::breaks_width(1)) +
      ggplot2::labs(x = "Age (PCW)", y = "Confirmed cerebellum groups", title = "Developmental coverage",
        subtitle = "Excludes provisional FB/GZFB groups; PCW awaits HPC reconciliation"))
    run <- coverage[known_age & !is.na(coverage$run_id), , drop = FALSE]
    if (nrow(run)) {
      add("atlas_age_run_confounding", ggplot2::ggplot(run, ggplot2::aes(PCW, run_id)) +
        ggplot2::geom_count(color = "#247BA0") + ggplot2::labs(x = "Age (PCW)", y = "Verified run ID", size = "Sections", title = "Age versus technical run"))
    } else record("atlas_age_run_confounding", "blocked", "No verified run IDs with donor/age data")
  } else record("atlas_donors_by_age", "blocked", "No verified donor/PCW pairs; no identities inferred from names")
  retention <- atlas_read_csv(files[4]); quantiles <- atlas_read_csv(files[5]); failures <- atlas_read_csv(files[6])
  populations <- atlas_read_csv(files[7]); controls <- atlas_read_csv(files[8])
  if (nrow(retention)) {
    atlas_keys(retention, c("sample_id", "stage"), "Retention")
    baseline <- retention[retention$stage == "cropped", ]
    if (any(!is.finite(baseline$n_cells) | baseline$n_cells <= 0)) stop("Invalid retention baseline.")
    retention$percent_of_cropped <- 100 * retention$n_cells / baseline$n_cells[match(retention$sample_id, baseline$sample_id)]
    retention$stage <- factor(retention$stage, c("cropped", "initial_qc", "final"))
    add("atlas_cell_retention", ggplot2::ggplot(retention, ggplot2::aes(stage, percent_of_cropped, group = sample_id)) +
      ggplot2::geom_line(alpha = .35, color = "#247BA0") + ggplot2::geom_point(size = .8) +
      ggplot2::labs(x = NULL, y = "% of cropped cerebellar cells", title = "Cell retention through the atlas", subtitle = "Regional branches are not successive filtering stages"))
    atlas_write(retention, file.path(plot_dir, "cell_retention_source.csv"))
    area <- anatomy[as.character(anatomy$compartment) == "whole_cerebellar_tissue" & anatomy$status == "present" &
      !is.na(anatomy$area_mm2) & anatomy$area_mm2 > 0, c("sample_id", "area_mm2")]
    density <- merge(retention[retention$stage == "final", ], area, by = "sample_id")
    if (nrow(density)) {
      density$cells_per_mm2 <- density$n_cells / density$area_mm2
      add("atlas_cerebellar_cell_density", ggplot2::ggplot(density, ggplot2::aes(cells_per_mm2, sample_id)) +
        ggplot2::geom_point(color = "#247BA0") + ggplot2::labs(x = "Final retained cells per measured cerebellar mm2", y = NULL,
          title = "Sampled cerebellar cell density", subtitle = "Requires reviewed crop area; includes filtering and detection effects"), 10, max(5, nrow(density)*.22))
      atlas_write(density, file.path(plot_dir, "cell_density_source.csv"))
    } else record("atlas_cerebellar_cell_density", "blocked", "No reviewed positive cerebellar crop areas")
  } else record("atlas_cell_retention", "blocked", "No complete extraction bundles")
  if (nrow(quantiles)) {
    quantiles$stage <- factor(quantiles$stage, c("cropped", "initial_qc", "final"))
    add("atlas_qc_distributions", ggplot2::ggplot(quantiles, ggplot2::aes(stage, median, group = sample_id)) +
      ggplot2::geom_line(alpha = .3, color = "#247BA0") + ggplot2::geom_point(size = .7) +
      ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::labs(x = NULL, y = "Per-section median", title = "QC across processing stages", subtitle = "Cell-level 5th/25th/75th/95th percentiles are retained in quantiles.csv"), 10, 7)
  } else record("atlas_qc_distributions", "blocked", "No extracted metric summaries")
  if (nrow(failures)) {
    f <- failures[failures$kind == "individual_nonexclusive", ]
    n <- retention[retention$stage == "cropped", ]
    f$percent <- 100 * f$n_cells / n$n_cells[match(f$sample_id, n$sample_id)]
    add("atlas_qc_failure_reasons", ggplot2::ggplot(f, ggplot2::aes(reason, percent)) +
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position = ggplot2::position_jitter(width = .12, seed = 42), alpha = .5, size = .8) +
      ggplot2::labs(x = NULL, y = "% of cropped cells per section", title = "Initial QC failure reasons", subtitle = "Reasons overlap; exclusive combinations are exported separately") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)), 9, 5)
  } else record("atlas_qc_failure_reasons", "blocked", "No reconstructed QC flags")
  if (nrow(controls)) {
    add("atlas_control_availability", ggplot2::ggplot(controls, ggplot2::aes(assay, factor(sample_id, levels = y_order), fill = imported)) +
      ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_manual(values = c(`FALSE` = "#E6E6E6", `TRUE` = "#247BA0")) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Imported", title = "Control-assay availability", subtitle = "Missing assays do not establish zero background"), 9, max(6, nrow(coverage)*.21))
  } else record("atlas_control_availability", "blocked", "No assay availability records")
  audit <- atlas_read_csv(files[9])
  if (nrow(audit)) {
    metrics <- c("raw_area_mismatch_cells", "raw_multiple_nuclei_fraction", "raw_no_nuclei_fraction")
    a <- do.call(rbind, lapply(metrics, function(k) data.frame(sample_id = audit$sample_id, metric = k, value = audit[[k]])))
    a <- a[is.finite(a$value), ]
    if (nrow(a)) add("atlas_morphology_audit", ggplot2::ggplot(a, ggplot2::aes(value, sample_id)) +
      ggplot2::geom_point(color = "#247BA0") + ggplot2::facet_wrap(~metric, scales = "free_x") +
      ggplot2::labs(x = "Raw-table morphology diagnostic", y = NULL, title = "Cerebellar crop morphology audit", subtitle = "Diagnostic flags do not automatically remove cells"), 12, max(5, nrow(audit)*.22))
    else record("atlas_morphology_audit", "blocked", "Raw morphology values unavailable")
  } else record("atlas_morphology_audit", "blocked", "No raw morphology extraction")
  if (nrow(populations)) {
    atlas_keys(populations, c("sample_id", "label"), "Population counts")
    if (any(!is.finite(populations$n_cells) | populations$n_cells < 0)) stop("Invalid population counts.")
    p <- merge(populations, coverage[, c("sample_id", "donor_id", "donor_verified", "PCW")], by = "sample_id", all.x = TRUE, sort = FALSE)
    complete_samples <- ready$sample_id[ready$extraction_complete]
    eligible_donors <- unique(coverage$donor_id[verified])
    eligible_donors <- eligible_donors[vapply(eligible_donors, function(d) all(coverage$sample_id[which(coverage$donor_id == d)] %in% complete_samples), logical(1))]
    p <- p[!is.na(p$donor_verified) & p$donor_verified == "yes" & !is.na(p$donor_id) & p$donor_id %in% eligible_donors, ]
    if (nrow(p)) {
      donor_counts <- aggregate(n_cells ~ donor_id + label, p, sum)
      grid <- expand.grid(donor_id = unique(p$donor_id), label = unique(populations$label), stringsAsFactors = FALSE)
      donor_counts <- merge(grid, donor_counts, by = c("donor_id", "label"), all.x = TRUE)
      donor_counts$n_cells[is.na(donor_counts$n_cells)] <- 0
      donor_counts$sampled_tissue_fraction <- donor_counts$n_cells / ave(donor_counts$n_cells, donor_counts$donor_id, FUN = sum)
      donor_counts$PCW <- coverage$PCW[match(donor_counts$donor_id, coverage$donor_id)]
      support <- aggregate(n_cells ~ label, donor_counts, function(x) sum(x > 0))
      names(support)[2] <- "donors_detected"
      add("atlas_population_donor_support", ggplot2::ggplot(support, ggplot2::aes(donors_detected, reorder(label, donors_detected))) +
        ggplot2::geom_col(fill = "#247BA0") + ggplot2::labs(x = "Verified donors with at least 1 observed cell", y = NULL, title = "Population detection across donors", subtitle = "Detection is not independent validation of identity"), 9, max(5, nrow(support)*.25))
      atlas_write(donor_counts, file.path(plot_dir, "donor_population_counts.csv"))
      # Zero-filling applies only to donors with a complete extraction bundle.
      if (length(unique(donor_counts$donor_id)) >= 3L) {
        sensitivity <- do.call(rbind, lapply(unique(donor_counts$donor_id), function(d) {
          x <- donor_counts[donor_counts$donor_id != d, ]
          a <- aggregate(sampled_tissue_fraction ~ label, x, mean)
          a$omitted_donor <- d; a
        }))
        atlas_write(sensitivity, file.path(plot_dir, "leave_one_donor_out_composition.csv"))
        add("atlas_composition_sensitivity", ggplot2::ggplot(sensitivity, ggplot2::aes(sampled_tissue_fraction, reorder(label, sampled_tissue_fraction))) +
          ggplot2::geom_point(alpha = .5, size = 1) + ggplot2::labs(x = "Mean donor fraction after omitting one donor", y = NULL,
            title = "Descriptive composition sensitivity", subtitle = "Sections pooled within donor; unequal anatomical sampling remains a limitation"), 10, max(5, nrow(support)*.25))
      } else record("atlas_composition_sensitivity", "blocked", "Fewer than three verified donors with extracted populations")
    } else record("atlas_population_donor_support", "blocked", "No verified donor mapping for extracted cells")
  } else record("atlas_population_donor_support", "blocked", "No final-population extraction")
  atlas_write(manifest, file.path(plot_dir, "plot_manifest.csv"))
  atlas_write(atlas_provenance(files), file.path(plot_dir, "plot_sources.csv"))
  message("Generated ", length(plots), " panels. Blocked panels are listed in plot_manifest.csv.")
}
if (sys.nframe() == 0L) main()
