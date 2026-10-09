# Focused lineage-marker panels and plotting utilities for raw VZ/RL
# subcluster review. Marker availability is checked against the Xenium assay so
# the review remains usable if the targeted panel changes.

rl_focused_review_markers <- list(
  "RL VZ" = c("LMX1A", "WLS", "SOX2", "PAX6", "CALCB", "CRYAB", "COL2A1"),
  "RL SVZ / UBC progenitors" = c(
    "EOMES", "OTX2", "RSPO3", "MKI67", "TOP2A", "PCNA"
  ),
  "UBC differentiation" = c(
    "TRPC3", "GRIA2", "CBLN1", "CALB2", "ZFHX4", "CNTNAP5"
  ),
  "GCP proliferation" = c(
    "ATOH1", "SFRP1", "CXCR4", "MYCN", "CCND2", "MCM2"
  ),
  "GC differentiation / migration" = c(
    "NEUROD1", "NHLH1", "PRPH", "CNTN2", "DCX", "TUBB3", "RELN"
  ),
  "Mature GC" = c(
    "SLC17A7", "GRM4", "GABRA6", "CBLN3", "RBFOX3", "NEFM"
  ),
  "eCN" = c("MEIS2", "TBR1", "SLC17A6", "TAC1", "NEFL", "NTNG2"),
  "BG / radial glia" = c(
    "SLC1A3", "GDF10", "TNC", "ATP1A2", "HOPX", "FABP7", "SOX9"
  ),
  "PC / GABAergic" = c(
    "CALB1", "RORA", "PCP4", "CA8", "FOXP2", "GAD1", "GAD2"
  ),
  "OPC / ependymal" = c(
    "OLIG1", "OLIG2", "PDGFRA", "SOX10", "FOXJ1", "PIFO", "TTR"
  ),
  "Immune / endothelial" = c(
    "P2RY12", "AIF1", "TYROBP", "C1QB", "CLDN5", "PECAM1", "KDR", "EPAS1"
  )
)

vz_focused_review_markers <- list(
  "VZ progenitors" = c(
    "PTF1A", "PRDM13", "ASCL1", "HES1", "DLL1", "HES6", "SOX1", "SOX2"
  ),
  "Cell cycle" = c(
    "MKI67", "TOP2A", "PCNA", "RRM2", "TYMS",
    "CDC20", "CCNB1", "CCNB2", "BUB1", "UBE2C"
  ),
  "Purkinje core" = c(
    "RORA", "PCP4", "CA8", "CALB1", "ITPR1", "SKOR2", "DAB1", "FOXP2"
  ),
  "Purkinje patterning" = c(
    "FOXP1", "FOXP4", "EBF1", "EBF2", "EN1",
    "NDNF", "TRPC3", "PCDH10", "EPHA4", "TMEFF2"
  ),
  "GABA progenitors" = c(
    "PAX2", "NEUROG2", "TFAP2A", "TFAP2B", "GAD1", "GAD2"
  ),
  "Golgi / MLI" = c(
    "SLC6A5", "NXPH2", "GRIK3", "SLC6A1", "SP9",
    "PVALB", "SST", "ERBB4", "GRIK1", "KIT", "KCNA2", "SORCS3", "NXPH1"
  ),
  "iCN / GABA maturation" = c(
    "LHX1", "LHX5", "SOX14", "DMBX1", "MEIS1",
    "TAC1", "SLC24A4", "CHRM2"
  ),
  "Radial / cycling glia" = c(
    "NES", "FABP7", "PTPRZ1", "EGFR", "BCAN", "SOX9", "NOTCH1", "NOTCH2"
  ),
  "Bergmann glia" = c(
    "SLC1A3", "GDF10", "TNC", "ATP1A2", "HOPX",
    "CACNG5", "FAM107A", "FREM2", "ETV5", "MEGF10"
  ),
  "Astrocyte / ependymal" = c(
    "AQP4", "SLC6A11", "GFAP", "ALDH1L1", "FOXJ1",
    "PIFO", "TTR", "KCNJ10", "GJA1", "FGFR3"
  ),
  "Pre-OPC / OPC" = c(
    "OLIG1", "OLIG2", "SOX8", "PDGFRA", "SOX10",
    "GPR17", "CSPG4", "TNR", "PROM1"
  ),
  "RL VZ" = c(
    "LMX1A", "WLS", "PAX6", "CALCB", "CRYAB",
    "COL2A1", "OTX2", "RSPO3", "WNT3A", "TPBG"
  ),
  "Immature UBC" = c(
    "EOMES", "ZIC2", "ZIC4", "CACNA2D1", "CCND2"
  ),
  "Mature UBC" = c(
    "CBLN1", "RELN", "TRH", "SSTR2", "NEUROD1", "MAB21L1"
  ),
  "eCN" = c(
    "SLC17A6", "C1QL3", "MEIS2", "NTNG2", "NEFL", "CREG2"
  ),
  "Microglia" = c(
    "P2RY12", "AIF1", "TYROBP", "C1QB", "SLCO2B1", "OLFML3", "ITM2A"
  ),
  "Endothelial contrast" = c(
    "CLDN5", "PECAM1", "KDR", "EPAS1"
  ),
  "Granule contrast" = c("ATOH1", "SFRP1", "CXCR4", "MYCN", "NHLH1", "SLC17A7")
)

add_palette_markers_to_review <- function(
    focused_marker_sets, palette_marker_sets, palette_object_name) {
  validate_focused_marker_sets(focused_marker_sets)
  validate_focused_marker_sets(palette_marker_sets)
  if (length(palette_object_name) != 1L || is.na(palette_object_name) ||
      !nzchar(palette_object_name)) {
    stop("palette_object_name must be one non-empty string.")
  }

  focused_markers <- unlist(focused_marker_sets, use.names = FALSE)
  palette_additions <- lapply(
    palette_marker_sets,
    function(markers) markers[!markers %in% focused_markers]
  )
  palette_additions <- palette_additions[lengths(palette_additions) > 0L]
  if (length(palette_additions)) {
    names(palette_additions) <- paste0(
      "color_palette.R ", palette_object_name, ": ",
      names(palette_additions)
    )
  }

  combined_marker_sets <- c(focused_marker_sets, palette_additions)
  validate_focused_marker_sets(combined_marker_sets)
  combined_marker_sets
}

validate_focused_marker_sets <- function(marker_sets) {
  if (!is.list(marker_sets) || !length(marker_sets) ||
      is.null(names(marker_sets)) || any(!nzchar(names(marker_sets)))) {
    stop("Focused marker sets must be a named, non-empty list.")
  }
  markers <- unlist(marker_sets, use.names = FALSE)
  if (!length(markers) || anyNA(markers) || any(!nzchar(markers))) {
    stop("Focused marker sets contain missing or blank gene names.")
  }
  duplicated_markers <- unique(markers[duplicated(markers)])
  if (length(duplicated_markers)) {
    stop(
      "Focused marker genes must occur in exactly one group. Duplicated: ",
      paste(duplicated_markers, collapse = ", ")
    )
  }
  invisible(TRUE)
}

focused_marker_set_signature <- function(marker_sets) {
  validate_focused_marker_sets(marker_sets)
  paste(
    vapply(
      names(marker_sets),
      function(group) paste0(group, "=", paste(marker_sets[[group]], collapse = "|")),
      character(1)
    ),
    collapse = ";"
  )
}

write_focused_marker_dotplot <- function(
    object, marker_sets, assay, group_by, branch, resolution_label,
    availability_path, tiff_path, pdf_path, width = 21, height = 10) {
  validate_focused_marker_sets(marker_sets)
  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Focused marker DotPlot assay is absent: ", assay)
  }
  if (!group_by %in% colnames(object[[]])) {
    stop("Focused marker DotPlot grouping column is absent: ", group_by)
  }

  marker_groups <- rep(names(marker_sets), lengths(marker_sets))
  requested_markers <- unlist(marker_sets, use.names = FALSE)
  assay_features <- rownames(object[[assay]])
  availability <- data.frame(
    branch = branch,
    resolution = resolution_label,
    cluster_column = group_by,
    marker_group = marker_groups,
    gene = requested_markers,
    available = requested_markers %in% assay_features,
    stringsAsFactors = FALSE
  )
  write.csv(availability, availability_path, row.names = FALSE)

  available_marker_sets <- lapply(
    marker_sets,
    function(markers) markers[markers %in% assay_features]
  )
  empty_groups <- names(available_marker_sets)[lengths(available_marker_sets) == 0L]
  if (length(empty_groups)) {
    warning(
      branch, " focused marker group(s) have no genes in the Xenium assay: ",
      paste(empty_groups, collapse = ", ")
    )
  }
  available_marker_sets <- available_marker_sets[lengths(available_marker_sets) > 0L]
  if (!length(available_marker_sets)) {
    stop("None of the focused ", branch, " markers are present in the Xenium assay.")
  }
  available_markers <- unlist(available_marker_sets, use.names = FALSE)
  marker_group_lookup <- stats::setNames(
    rep(names(available_marker_sets), lengths(available_marker_sets)),
    available_markers
  )

  missing_markers <- availability$gene[!availability$available]
  if (length(missing_markers)) {
    message(
      branch, " focused DotPlot omits ", length(missing_markers),
      " marker(s) absent from the Xenium assay; see ", availability_path
    )
  }

  plot <- Seurat::DotPlot(
    object,
    # Passing a named list to DotPlot triggers Seurat's legacy facet_grid()
    # interface, which is incompatible with current ggplot2. Add the marker
    # group facets explicitly below instead.
    features = available_markers,
    assay = assay,
    group.by = group_by,
    cols = c("grey90", "#B2182B"),
    dot.scale = 6,
    cluster.idents = FALSE
  )
  plot$data$marker_group <- factor(
    unname(marker_group_lookup[as.character(plot$data$features.plot)]),
    levels = names(available_marker_sets)
  )
  if (anyNA(plot$data$marker_group)) {
    stop("Focused marker group lookup failed while constructing the DotPlot.")
  }
  plot <- plot +
    ggplot2::facet_grid(
      cols = ggplot2::vars(marker_group),
      scales = "free_x",
      space = "free_x"
    ) +
    Seurat::RotatedAxis() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, face = "italic"),
      axis.text.y = ggplot2::element_text(size = 9, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 9, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::ggtitle(
      paste0(branch, " focused lineage markers - resolution ", resolution_label)
    )

  Cairo::CairoTIFF(
    filename = tiff_path,
    width = width,
    height = height,
    units = "in",
    res = 600
  )
  print(plot)
  grDevices::dev.off()
  ggplot2::ggsave(
    pdf_path, plot, device = grDevices::cairo_pdf,
    width = width, height = height
  )
  invisible(availability)
}
