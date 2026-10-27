#!/usr/bin/env Rscript
# Small label/palette checks only; expected under 10 seconds. No analysis RDS.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/test_vascular_labels.R")
source(here("scripts", "R", "science_reviewed_mapping.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))
palette <- new.env(parent = globalenv())
source(here("scripts", "color_palette.R"), local = palette)
source_labels <- factor(c("Endothelial", "Pericytes", "Vascular", "Granule"))
before <- source_labels
out <- harmonize_vascular_reference_labels(source_labels)
stopifnot(identical(source_labels, before),
          identical(out, c("Vascular", "Vascular", "Vascular", "Granule")),
          identical(map_corrected_broad_labels(c("Vascular", "Pericytes", "Unknown-01")), c("Vascular", "Vascular", "Unknown")))
for (x in list(palette$cluster_colors, palette$markers, palette$rl_palette, palette$rl_other_markers)) {
  stopifnot(!anyDuplicated(names(x)), !"Endothelial" %in% names(x))
}
stopifnot(!anyDuplicated(palette$celltype_order), !"Endothelial" %in% palette$celltype_order,
          !"Pericytes" %in% palette$celltype_order,
          "Pericytes" %in% palette$master_subcluster_order,
          all(c("CLDN5", "PECAM1", "PDGFRB", "ABCC9", "KCNJ8") %in% palette$markers$Vascular),
          identical(palette$cluster_colors[["Vascular"]], "#FF1A1A"))
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
expect_error(harmonize_vascular_reference_labels(c("Vascular", NA_character_)))
expect_error(harmonize_vascular_reference_labels(""))
# Old downstream names fail rather than silently creating a separate broad group.
expect_error(map_corrected_broad_labels("Endothelial"))
cat("PASS: reference harmonization, source preservation, Vascular broad identity, retained Pericytes fine identity, unique palettes and marker groups.\n")
