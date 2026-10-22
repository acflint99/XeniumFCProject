#!/usr/bin/env Rscript
# Lightweight regression checks only; no Seurat objects or biological data.
# Estimated runtime: under 10 seconds, based on small in-memory marker tables.
source(file.path("R", "nonspatial_merge_plots.R"))
levels <- c("RL", "Granule", "Unresolved")
top <- data.frame(cluster = c("RL", "Granule"), gene = c("OTX2", "PAX6"))
notes <- character()
summary <- withCallingHandlers(
  nonspatial_heatmap_marker_summary(top, levels),
  message = function(m) { notes <<- c(notes, conditionMessage(m)); invokeRestart("muffleMessage") }
)
stopifnot(
  identical(summary$identity, levels),
  identical(summary$markers_selected, c(1L, 1L, 0L)),
  summary$status[summary$identity == "Unresolved"] == "no_qualifying_positive_markers_cells_retained",
  any(grepl("retaining its cells", notes, fixed = TRUE))
)
complete <- nonspatial_heatmap_marker_summary(
  rbind(top, data.frame(cluster = "Unresolved", gene = "GENE")), levels
)
stopifnot(all(complete$markers_selected == 1L))
# Missing a named biological group still fails, as do entirely empty results.
missing_known <- tryCatch(nonspatial_heatmap_marker_summary(top[1, , drop = FALSE], levels), error = identity)
empty <- tryCatch(nonspatial_heatmap_marker_summary(top[FALSE, ], levels), error = identity)
stopifnot(inherits(missing_known, "error"), grepl("Granule", conditionMessage(missing_known)), inherits(empty, "error"))
cat("PASS: absent Unresolved markers are reported, all groups remain represented, and unexpected missing marker groups still fail.\n")
