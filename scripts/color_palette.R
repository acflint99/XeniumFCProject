# =====================================================
# Master Cell Type Color Palette
# Cerebellar Xenium Project
# =====================================================
library(scales) # Required for alpha()
source(here::here("scripts", "R", "consensus_labels.R"))

####broad cluster label colors & order####
# Named vector: names MUST match cluster labels exactly
cluster_colors <- c(
  Purkinje    = "#1E90FF",  # vivid blue
  GABA        = "#39FF14",  # neon green
  RL          = "#FF00FF",  # bright magenta
  UBC         = "#FF7F00",  # bright neon orange
  Granule     = "#00CED1",  # neon cyan / teal
  VZ          = "#00FFD5",  # neon teal
  Glia        = "#7FBFFF",  # pastel blue
  OPC         = "#FFD700",  # neon yellow
  Meninges    = "#9370DB",  # medium purple
  Endothelial = "#FF1A1A",  # deep red
  Immune      = "#E066FF",   # neon lavender
  Unknown = "grey80"
)

# Optional: function to validate palette
validate_palette <- function(labels) {
  missing <- setdiff(labels[!is_unknown_consensus_label(labels)], names(cluster_colors))
  if (length(missing) > 0) {
    stop(
      paste(
        "Missing colors for:",
        paste(missing, collapse = ", ")
      )
    )
  }
}

# Detailed Unknown clusters receive deterministic grey shades. The numbering
# order controls the shade order, so plots remain stable across reruns.
resolve_cluster_colors <- function(labels) {
  labels <- as.character(labels)
  validate_palette(labels)
  colors <- cluster_colors[labels]
  unknown_labels <- labels[is_unknown_consensus_label(labels)]
  if (length(unknown_labels)) {
    unknown_numbers <- suppressWarnings(
      as.integer(sub("^Unknown-", "", unknown_labels))
    )
    unknown_numbers[is.na(unknown_numbers)] <- seq_len(sum(is.na(unknown_numbers)))
    grey_values <- grDevices::gray(
      pmin(0.35 + 0.08 * (unknown_numbers - 1L), 0.83)
    )
    colors[match(unknown_labels, labels)] <- grey_values
  }
  stats::setNames(unname(colors), labels)
}

# Regional VZ/RL objects preserve the original consensus_label for joins and
# append the sample ID only in regional_consensus_label. These helpers keep all
# known biological identities on the master palette while assigning arbitrary,
# deterministic colors to the sample-qualified Unknown identities.
is_sample_qualified_unknown_consensus_label <- function(labels) {
  labels <- trimws(as.character(labels))
  !is.na(labels) & grepl("^Unknown(?:-[0-9]{2,})?_.+$", labels, perl = TRUE)
}

order_regional_consensus_levels <- function(labels, known_order) {
  labels <- unique(trimws(as.character(labels)))
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Regional consensus labels must not be blank or missing.")
  }
  unknown <- is_unknown_consensus_label(labels) |
    is_sample_qualified_unknown_consensus_label(labels)
  known_labels <- labels[!unknown]
  unknown_labels <- sort(labels[unknown])
  c(
    intersect(known_order, known_labels),
    setdiff(known_labels, known_order),
    unknown_labels
  )
}

resolve_regional_consensus_colors <- function(labels) {
  labels <- as.character(labels)
  qualified_unknown <- is_sample_qualified_unknown_consensus_label(labels)
  standard_labels <- labels[!qualified_unknown]
  colors <- stats::setNames(rep(NA_character_, length(labels)), labels)
  if (length(standard_labels)) {
    colors[standard_labels] <- resolve_cluster_colors(standard_labels)
  }
  qualified_labels <- labels[qualified_unknown]
  if (length(qualified_labels)) {
    colors[qualified_labels] <- grDevices::hcl.colors(
      length(qualified_labels), "Dynamic"
    )
  }
  colors
}

resolve_sample_qualified_unknown_colors <- function(labels) {
  labels <- as.character(labels)
  colors <- stats::setNames(
    grDevices::hcl.colors(length(labels), "Dynamic"), labels
  )
  if ("Non-Unknown" %in% labels) colors[["Non-Unknown"]] <- "#D9D9D9"
  colors
}

#lock factor order
celltype_order <- c(
  "RL", "UBC", "Granule",
  "VZ", "Purkinje", "GABA", "Glia", "OPC", "Meninges",
  "Endothelial", "Immune", "Unknown"
)

# Shared scales and figure styling for all current workflow DotPlots.
source(here::here("scripts", "R", "dotplot_helpers.R"))

standardize_broad_dotplot <- function(plot, label_order = celltype_order) {
  standardize_workflow_dotplot(plot) +
    ggplot2::scale_y_discrete(limits = rev(label_order), drop = FALSE)
}

markers <- list(
  "RL" = c("MKI67", "LTBP1", "OTX2"),
  "UBC" = c("EOMES"),
  "Granule" = c("ATOH1", "PAX6", "NEUROD1", "RELN"),
  "VZ" = c("PRDM13", "ASCL1"),
  "Purkinje" = c("FOXP2", "CALB1", "DAB1"),
  "GABA" = c("PAX2", "GAD1", "GAD2"),
  "Glia" = c("SOX9", "TNC"),
  "OPC" = c("PDGFRA", "OLIG1"),
  "Meninges" = c("FOXC1", "SLC7A11"),
  "Endothelial" = c("CLDN5", "PECAM1"),
  "Immune" = c("P2RY12")
)


####VZ subcluster label colors & order####
vz_markers <- list(
  "VZPs" = c(
    "PTF1A", "PRDM13", "ASCL1", "DLL1", "HES6", "SOX1", "SOX2"
  ),
  "Cell cycle" = c(
    "MKI67", "TOP2A", "PCNA", "RRM2", "TYMS",
    "CDC20", "CCNB1", "CCNB2", "BUB1", "UBE2C"
  ),
  "Purkinje core" = c(
    "RORA", "PCP4", "CA8", "CALB1", "ITPR1", "SKOR2", "DAB1", "FOXP2"
  ),
  "Purkinje subtype" = c(
    "FOXP1", "FOXP4", "EBF1", "EBF2", "EN1",
    "NDNF", "TRPC3", "PCDH10", "EPHA4", "TMEFF2"
  ),
  "GABA progenitors" = c(
    "PAX2", "NEUROG2", "TFAP2A", "TFAP2B", "GAD1", "GAD2"
  ),
  "Golgi / MLI" = c(
    "SLC6A5", "NXPH2", "GRIK3", "SLC6A1", "SP9",
    "PVALB", "SST", "ERBB4", "GRIK1", "KIT"
  ),
  "iCN" = c("SOX14", "DMBX1", "MEIS1", "TAC1", "SLC24A4", "CHRM2"),
  "RG / cycling glia" = c(
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
  "Immature UBC" = c("EOMES", "ZIC2", "ZIC4", "CACNA2D1", "CCND2"),
  "Mature UBC" = c("CBLN1", "RELN", "TRH", "SSTR2", "NEUROD1", "MAB21L1"),
  "eCN" = c("SLC17A6", "C1QL3", "MEIS2", "NTNG2", "NEFL", "CREG2"),
  "Microglia" = c(
    "P2RY12", "AIF1", "TYROBP", "C1QB", "SLCO2B1", "OLFML3", "ITM2A"
  ),
  "Endothelial contrast" = c("CLDN5", "PECAM1", "KDR", "EPAS1")
)

# 1. Define the subcluster order (shared across scripts)
vz_subcluster_order <- c(
  "VZPs", "Cycling VZPs",
  "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma",
  "Pre-OPCs", "OPCs",
  "RL VZ", "Immature UBCs", "Mature UBCs", "eCN",
  "Cycling Cells", "Microglia"
)

# Compact production panels for readable lineage-specific DotPlots. These are
# deliberately separate from the broader vz_markers review/merged panel.
vz_lineage_subclusters <- c(
  "VZPs", "Cycling VZPs",
  "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma",
  "Pre-OPCs", "OPCs"
)
vz_other_subclusters <- c(
  "RL VZ", "Immature UBCs", "Mature UBCs", "eCN",
  "Cycling Cells", "Microglia"
)

vz_compact_lineage_markers <- list(
  "VZ" = c("PRDM13", "ASCL1", "SOX2"),
  "Purkinje" = c(
    "DAB1", "FOXP4", "CALB1", "NDNF", "EBF1", "PCDH10", "EPHA4"
  ),
  "GABA" = c("PAX2", "TFAP2A", "KIT", "NXPH2", "SP9", "SOX14"),
  "Glia" = c("NES", "FABP7", "SLC1A3", "TNC", "BCAN", "AQP4", "FOXJ1"),
  "OPCs" = c("EGFR", "OLIG1", "OLIG2", "PDGFRA", "SOX10"),
  "Cycling" = c("MKI67", "TOP2A", "PCNA", "TYMS")
)

vz_other_markers <- list(
  "RL VZ" = c("CALCB", "OTX2", "RSPO3"),
  "UBCs" = c("EOMES", "CCND2", "CBLN1", "RELN"),
  "eCN" = c("MEIS2", "SLC17A6"),
  "Cycling Cells" = c("MKI67", "TOP2A"),
  "Microglia" = c("P2RY12")
)

# Shared grouped panel for the Stage-6 all-subcluster DotPlot and every Stage-8
# sample DotPlot. MKI67 and TOP2A are already present under Cycling, so the
# duplicate Cycling Cells marker group is omitted without removing either gene.
vz_all_compact_markers <- c(
  vz_compact_lineage_markers,
  vz_other_markers[c("RL VZ", "UBCs", "eCN", "Microglia")]
)

validate_compact_marker_panel <- function(
    marker_sets, expected_labels, panel_name, min_markers, max_markers) {
  if (!identical(names(marker_sets), expected_labels)) {
    stop(panel_name, " labels must exactly match their configured cluster order.")
  }
  marker_counts <- lengths(marker_sets)
  if (any(marker_counts < min_markers | marker_counts > max_markers)) {
    stop(
      panel_name, " must provide ", min_markers, "-", max_markers,
      " marker(s) per cluster."
    )
  }
  markers <- unlist(marker_sets, use.names = FALSE)
  if (anyNA(markers) || any(!nzchar(markers)) || anyDuplicated(markers)) {
    stop(panel_name, " contains blank, missing, or duplicated marker genes.")
  }
  invisible(TRUE)
}

validate_compact_marker_panel(
  vz_compact_lineage_markers,
  c("VZ", "Purkinje", "GABA", "Glia", "OPCs", "Cycling"),
  "VZ compact lineage panel", 3L, 7L
)
validate_compact_marker_panel(
  vz_other_markers,
  c("RL VZ", "UBCs", "eCN", "Cycling Cells", "Microglia"),
  "VZ Other panel", 1L, 4L
)
validate_compact_marker_panel(
  vz_all_compact_markers,
  c(
    "VZ", "Purkinje", "GABA", "Glia", "OPCs", "Cycling",
    "RL VZ", "UBCs", "eCN", "Microglia"
  ),
  "VZ all-subcluster compact panel", 1L, 7L
)
if (!setequal(
  unlist(vz_all_compact_markers, use.names = FALSE),
  unique(c(
    unlist(vz_compact_lineage_markers, use.names = FALSE),
    unlist(vz_other_markers, use.names = FALSE)
  ))
)) {
  stop("VZ all-subcluster compact panel must contain the exact 42-gene union.")
}
if (!setequal(c(vz_lineage_subclusters, vz_other_subclusters), vz_subcluster_order)) {
  stop("VZ lineage and Other DotPlot groups do not cover the full VZ taxonomy.")
}

# 2. Define the Master Palette
# You can use hex codes or R color names.
vz_palette <- c(
  # --- Purkinje Lineage (High Contrast Blues/Purples) ---
  "VZPs"                = "#BBFFFF", # Neon Teal (Original VZ) - Very Bright
  "Cycling VZPs"        = "#008C8C", # Dark teal within the VZP gradient
  "Maturing PCs"       = "#1E90FF", # Vivid Blue (Original PC) - Medium
  "Early-born PCs" = "#CC66FF", # Medium Purple (Original UBC) - Medium Dark
  "Late-born PCs"         = "#00CED1", # Indigo - Darkest (Pops against Early/Migrating)
  "Patterning PCs"         = "#0047AB", # Indigo - Darkest (Pops against Early/Migrating)
  
  # --- GABAergic Lineage (Revised for Maximum Contrast) ---
  "GABA Progenitors"   = "#CCFF00", 
  "Golgi Cells"        = "#00FF7F", 
  "MLIs"  = "#32CD32", 
  "iCN"                = "#006400", 
  
  # --- Glial Lineage (blue through violet gradient) ---
  "RG Progenitors"     = "#7FBFFF", # Pastel Blue (Original Glia)
  "Cycling Glia"       = "#A25FD1", # Blue-purple bridge from RG to BG
  "BG"                 = "#8A2BE2", # Neon Lavender (Original Immune)
  "Astrocytes/Ependyma"         = "#E6B3FF", # Light lavender endpoint
  
  # --- OPC, RL/UBC, and off-lineage populations ---
  "Pre-OPCs"           = "#C99700", # Dark mustard within the OPC gradient
  "OPCs"               = "#FFD700", # Gold/Neon Yellow
  "RL VZ"              = "#C2185B", # Dark magenta; contrasts with light-pink RL SVZ
  "Immature UBCs"      = "#FFA600", # Shared with RL palette
  "Mature UBCs"        = "#8B0000", # Shared with RL palette
  "Cycling Cells" = "#FFFF00",
  "eCN"                = "#FF4500", # Orange Red (Moved from UBC)
  "Microglia"          = "#E066FF"  # Shared with RL palette
)
vz_palette["NA"] <- alpha("gray20", 0.3)

missing_vz_colors <- setdiff(vz_subcluster_order, names(vz_palette))
if (length(missing_vz_colors)) {
  stop(
    "VZ subcluster labels lack palette colors: ",
    paste(missing_vz_colors, collapse = ", ")
  )
}

# 3. Define a background/null color for highlight plots
bg_color <- alpha("gray20", 0.3)

####RL subcluster label colors & order####
rl_subcluster_order <- c(
  "RL SVZ", "Immature UBCs", "Mature UBCs",
  "Prolif GCPs", "Maturing GCPs",
  "Differentiating GCs", "Migrating GCs", "Mature GCs",
  "BG", "Maturing PCs",
  "eCN", "GABA Progenitors", "Microglia", "Endothelial"
)

rl_palette <- c(
  # --- Rhombic Lip (High-Contrast Cyans/Blues) ---
  "RL SVZ"               = "#FF4FA3", # Vivid magenta; distinct from Immature UBC orange
  
  # --- The Transition (The "Standout" Bridge) ---
    "Prolif GCPs"    = "#FFD700", # Neon Green
  "Maturing GCPs"       = "#00FF7F", # Deep Emerald
  "Differentiating GCs"      = "#006400", # Green-Yellow (Bright pop to distinguish from RL)
  "Migrating GCs" = "#00CED1",
  "Mature GCs" = "#0047AB",
  
  "Immature UBCs" = "#FFA600", # Vivid orange; distinct from RL SVZ and eCN
  "Mature UBCs"        = "#8B0000", # Dark Red
  
  #"RG Progenitors"     = "#7FBFFF", # Light Blue
  "BG"         = "#8A2BE2", # Magenta
  "Maturing PCs" = "#1E90FF",

  # --- Reviewed off-lineage populations retained for provenance ---
  "eCN" = "#FF4500",
  "GABA Progenitors" = "#CCFF00",
  "Microglia" = "#E066FF",
  "Endothelial" = "#FF1A1A"

)
rl_palette["NA"] <- alpha("gray20", 0.3)

missing_rl_colors <- setdiff(rl_subcluster_order, names(rl_palette))
if (length(missing_rl_colors)) {
  stop(
    "RL subcluster labels lack palette colors: ",
    paste(missing_rl_colors, collapse = ", ")
  )
}

rl_markers <- list(
  "RL" = c(
    "MKI67", "RSPO3", "OTX2"
  ),
  "UBC" = c(
    "EOMES", "CACNA2D1", "SEMA3A", "PTPRM", "TRPC3"
  ),
  "Granule" = c(
    "TOP2A", "ATOH1", "CXCR4", "SFRP1",
    "NHLH2", "CNTN2", "NTF3", "GABRA6"
  ),
  "Other" = c(
    "SLC1A3", "CALB1", "SLC17A6", "PAX2", "P2RY12", "CLDN5"
  )
)

rl_lineage_subclusters <- c(
  "RL SVZ", "Immature UBCs", "Mature UBCs",
  "Prolif GCPs", "Maturing GCPs", "Differentiating GCs",
  "Migrating GCs", "Mature GCs"
)
rl_other_subclusters <- c(
  "BG", "Maturing PCs", "eCN", "GABA Progenitors", "Microglia", "Endothelial"
)

rl_compact_lineage_markers <- list(
  "RL" = c("MKI67", "RSPO3", "OTX2"),
  "UBC" = c("EOMES", "CACNA2D1", "TRPC3"),
  "Granule" = c("ATOH1", "NHLH2", "GABRA6")
)

rl_other_markers <- list(
  "BG" = c("HOPX", "FAM107A", "FREM2"),
  "Maturing PCs" = c("CALB1", "ITPR1", "DAB1"),
  "eCN" = c("SLC17A6", "MEIS2", "NTNG2"),
  "GABA Progenitors" = c("PAX2", "GAD1", "GAD2"),
  "Microglia" = c("P2RY12", "SLCO2B1", "OLFML3"),
  "Endothelial" = c("CLDN5", "PECAM1", "KDR")
)

validate_compact_marker_panel(
  rl_compact_lineage_markers, c("RL", "UBC", "Granule"),
  "RL compact lineage panel", 3L, 3L
)
validate_compact_marker_panel(
  rl_other_markers, rl_other_subclusters,
  "RL Other panel", 2L, 3L
)
if (!setequal(c(rl_lineage_subclusters, rl_other_subclusters), rl_subcluster_order)) {
  stop("RL lineage and Other DotPlot groups do not cover the full RL taxonomy.")
}

# Residual-only identities that are not present in the reviewed VZ or RL maps.
# Keep them in the shared lookup so downstream combined plots cannot silently
# convert an approved residual label to NA.
residual_palette <- c(
  "Pericytes" = "#8B4513"
)

# Merge all palettes into one master lookup for the combined column
# using 'c()' on named vectors merges them
subcluster_palette <- c(vz_palette, rl_palette, residual_palette, cluster_colors)

# Remove any duplicate names (e.g., if "Purkinje" is in both) 
# to keep the vector clean
subcluster_palette <- subcluster_palette[!duplicated(names(subcluster_palette))]

master_subcluster_order <- c(
  "VZPs", "Cycling VZPs",
  "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma",
  "Pre-OPCs", "OPCs",
  "RL VZ", "RL SVZ", "Immature UBCs", "Mature UBCs", "eCN",
  "Prolif GCPs", "Maturing GCPs",
  "Differentiating GCs", "Migrating GCs", "Mature GCs",
  "Cycling Cells", "Pericytes", "Meninges", "Microglia", "Endothelial",
  "Immune", "Unknown"
)
