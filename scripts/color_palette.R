# =====================================================
# Master Cell Type Color Palette
# Cerebellar Xenium Project
# =====================================================
library(scales) # Required for alpha()
source(here::here("scripts", "R", "consensus_labels.R"))

####broad cluster label colors & order####
# Named vector: names MUST match cluster labels exactly
# Approved neon families: Granule green, GABA teal, RL warm pink, VZ cool pink.
# Shared identities must use the same colors in the VZ and RL palettes below.
cluster_colors <- c(
  Purkinje    = "#1E90FF",  # vivid blue
  GABA        = "#00FFD5",  # neon turquoise; matches GABA Progenitors
  RL          = "#E00065",  # raspberry pink; matches RL SVZ
  UBC         = "#FFC48C",  # pastel peach; matches Immature UBCs
  Granule     = "#39FF14",  # neon green
  VZ          = "#FF00D4",  # cool neon pink
  Glia        = "#8A2BE2",  # electric violet; matches BG
  OPC         = "#FFD700",  # neon yellow
  Meninges    = "#9370DB",  # medium purple
  Endothelial = "#FF1A1A",  # vivid red; retained unchanged
  Immune      = "#E066FF",   # neon lavender
  eCN         = "#FF4500",  # approved new broad identity
  Pericytes   = "#8B4513",  # approved new broad identity
  "Cycling Cells" = "#FFFF00",
  Unresolved  = "#707070",
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
  "RL", "UBC", "Granule", "eCN",
  "VZ", "Purkinje", "GABA", "Glia", "OPC", "Meninges",
  "Endothelial", "Pericytes", "Immune",
  "Cycling Cells", "Unresolved", "Unknown"
)

# Shared scales and figure styling for all current workflow DotPlots.
source(here::here("scripts", "R", "dotplot_helpers.R"))
# Detect an incomplete HPC script sync before any caller loads large objects.
required_dotplot_helpers <- c(
  "workflow_dotplot", "workflow_dotplot_dimensions", "workflow_dotplot_grob"
)
missing_dotplot_helpers <- required_dotplot_helpers[!vapply(
  required_dotplot_helpers,
  function(name) exists(name, mode = "function", inherits = TRUE),
  logical(1)
)]
if (length(missing_dotplot_helpers)) {
  stop(
    "DotPlot helper file is incomplete or outdated; missing: ",
    paste(missing_dotplot_helpers, collapse = ", "),
    ". Sync scripts/R/dotplot_helpers.R together with scripts/color_palette.R ",
    "and the calling plotting script before retrying."
  )
}


if (!"strip_colors" %in% names(formals(workflow_dotplot_grob))) {
  stop("DotPlot helper is outdated. Sync scripts/R/dotplot_helpers.R with scripts/color_palette.R.")
}

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
  "Immune" = c("P2RY12"),
  "eCN" = c("SLC17A6", "MEIS2", "NEUROD2"),
  "Pericytes" = c("PDGFRB", "ABCC9", "KCNJ8"),
  "Cycling Cells" = c("TOP2A", "MKI67", "CENPF")
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
# Current approved VZ labels only. Legacy shared colors remain available for
# other branches; no cell identities or biological routing are changed here.
vz_subcluster_order <- c(
  "VZPs",
  "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "GABA Interneurons", "iCN",
  "RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma",
  "Pre-OPCs", "OPCs", "RL VZ", "Differentiating GCs", "eCN", "Unresolved"
)

# Unresolved cells remain visible in global/Other plots without a forced lineage.
vz_other_subclusters <- c("RL VZ", "Differentiating GCs", "eCN", "Unresolved")
vz_lineage_subclusters <- setdiff(vz_subcluster_order, vz_other_subclusters)
vz_dotplot_lineage_subclusters <- vz_lineage_subclusters
vz_dotplot_other_subclusters <- vz_other_subclusters

# Display panels aligned with the VZ taxonomy updated 2026-10-19. These do not
# change raw review panels, marker tests, labels, expression scaling, or QC.
vz_compact_lineage_markers <- list(
  "VZ" = c("PRDM13", "ASCL1"),
  "Purkinje" = c(
    "DAB1", "FOXP1", "FOXP4", "EBF2", "EBF1", "EPHA4"
  ),
  "GABA" = c(
    "TFAP2B", "NXPH2", "SLC6A5", "GRIK3", "SP9", "SOX14", "DMBX1"
  ),
  "Glia" = c(
    "SOX2", "TOP2A", "TNC", "AQP4", "FOXJ1"
  ),
  "OPCs" = c("EGFR", "OLIG1", "OLIG2", "PDGFRA")
)
vz_other_markers <- list(
  "RL VZ" = c("CALCB", "OTX2", "RSPO3"),
  "Granule" = c("NEUROD1", "CBLN1", "RELN"),
  "eCN" = c("MEIS2", "SLC17A6")
)
vz_all_compact_markers <- c(vz_compact_lineage_markers, vz_other_markers)

# Keep only panels represented by observed identities in each merged/sample
# plot. Never drop cells or metadata when hiding an unused display panel.
vz_dotplot_group_identities <- list(
  "VZ" = "VZPs",
  "Purkinje" = c("Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs"),
  "GABA" = c("GABA Progenitors", "Golgi Cells", "MLIs", "GABA Interneurons", "iCN"),
  "Glia" = c("RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma"),
  "OPCs" = c("Pre-OPCs", "OPCs"), "RL VZ" = "RL VZ",
  "Granule" = "Differentiating GCs", "eCN" = "eCN"
)
active_vz_dotplot_markers <- function(identities, scope = c("all", "lineage", "other")) {
  scope <- match.arg(scope)
  panels <- switch(scope, all = vz_all_compact_markers,
                   lineage = vz_compact_lineage_markers, other = vz_other_markers)
  keep <- vapply(names(panels), function(group) {
    any(vz_dotplot_group_identities[[group]] %in% identities)
  }, logical(1))
  panels[keep]
}

validate_compact_marker_panel <- function(
    marker_sets, expected_labels, panel_name, min_markers, max_markers) {
  if (!identical(names(marker_sets), expected_labels)) {
    stop(panel_name, " labels must exactly match their configured cluster order.")
  }
  marker_counts <- lengths(marker_sets)
  if (any(marker_counts < min_markers | marker_counts > max_markers)) {
    stop(panel_name, " must provide ", min_markers, "-", max_markers, " marker(s) per cluster.")
  }
  markers <- unlist(marker_sets, use.names = FALSE)
  if (anyNA(markers) || any(!nzchar(markers)) || anyDuplicated(markers)) {
    stop(panel_name, " contains blank, missing, or duplicated marker genes.")
  }
  invisible(TRUE)
}
validate_compact_marker_panel(
  vz_compact_lineage_markers, c("VZ", "Purkinje", "GABA", "Glia", "OPCs"),
  "VZ compact lineage panel", 2L, 10L
)
validate_compact_marker_panel(
  vz_other_markers, c("RL VZ", "Granule", "eCN"), "VZ Other panel", 2L, 5L
)
validate_compact_marker_panel(
  vz_all_compact_markers, names(vz_dotplot_group_identities),
  "VZ all-subcluster compact panel", 2L, 10L
)
if (!setequal(c(vz_lineage_subclusters, vz_other_subclusters), vz_subcluster_order) ||
    !setequal(c(vz_dotplot_lineage_subclusters, vz_dotplot_other_subclusters), vz_subcluster_order)) {
  stop("VZ plotting partitions do not cover the approved VZ taxonomy.")
}

# 2. Define the Master Palette
# You can use hex codes or R color names.
vz_palette <- c(
  # --- VZ progenitors: cool pinks; Purkinje cells: blues and purple ---
  "VZPs"                = "#FF00D4", # Cool neon pink, matching broad VZ
  "Cycling VZPs"        = "#B000A8", # Deeper purple-magenta
  "Maturing PCs"        = "#1E90FF", # Vivid blue
  "Early-born PCs"      = "#CC66FF", # Existing purple retained
  "Late-born PCs"       = "#00BFFF", # Sky blue; swapped with proposed Golgi color
  "Patterning PCs"      = "#0047AB", # Deep blue
  
  # --- GABAergic lineage: turquoise through teal ---
  "GABA Progenitors"   = "#00FFD5",
  "GABA Interneurons"  = "#00AFA0",
  "Golgi Cells"        = "#A6D6F5",
  "MLIs"               = "#00FF7F",
  "iCN"                = "#005A70",
  
  # --- Glial Lineage (blue through violet gradient) ---
  "RG Progenitors"     = "#7FBFFF", # Pastel Blue (Original Glia)
  "Cycling Glia"       = "#6868FF", # Blue-violet, distinct from Early-born PC purple
  "BG"                 = "#8A2BE2", # Electric violet
  "Astrocytes/Ependyma"         = "#E6B3FF", # Light lavender endpoint
  
  # --- OPC, RL/granule, and off-lineage populations ---
  "Pre-OPCs"           = "#C99700", # Dark mustard within the OPC gradient
  "OPCs"               = "#FFD700", # Gold/Neon Yellow
  "RL VZ"              = "#FFB3D9", # Light pink; contrasts with raspberry RL SVZ
  "Immature UBCs"      = "#FFC48C", # Pastel peach; shared with RL palette
  "Differentiating GCs" = "#32CD32", # Shared with the RL granule palette
  "Neuroblasts"         = "#FF66B2", # Working neuronal-state label; lineage unresolved
  "Cycling Cells" = "#FFFF00",
  "eCN"                = "#FF4500", # Orange Red (Moved from UBC)
  "Unresolved"         = "#808080", # Retained without a forced lineage
  "Immune"             = "#E066FF",
  "Microglia"          = "#E066FF"  # Legacy stored-label alias
)
vz_palette["NA"] <- alpha("gray20", 0.3)

vz_dotplot_strip_colors <- c(
  VZ = unname(cluster_colors["VZ"]), Neuroblasts = unname(vz_palette["Neuroblasts"]),
  Purkinje = unname(cluster_colors["Purkinje"]), GABA = unname(cluster_colors["GABA"]),
  Glia = unname(cluster_colors["Glia"]), OPCs = unname(cluster_colors["OPC"]),
  "RL VZ" = unname(vz_palette["RL VZ"]),
  Granule = unname(cluster_colors["Granule"]), eCN = unname(vz_palette["eCN"])
)

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
  "BG Progenitors", "Maturing PCs",
  "eCN", "GABA Progenitors", "Microglia", "Endothelial"
)

rl_palette <- c(
  # --- Rhombic lip: raspberry and pastel pinks ---
  "RL SVZ"             = "#E00065", # Raspberry pink; matches broad RL
  
  # --- Granule lineage: yellow-green through emerald ---
  "Prolif GCPs"        = "#CCFF00",
  "Maturing GCPs"      = "#39FF14",
  "Differentiating GCs" = "#32CD32",
  "Migrating GCs"      = "#00CC99",
  "Mature GCs"         = "#008C46",
  
  "Immature UBCs"      = "#FFC48C", # Pastel peach
  "Mature UBCs"        = "#FF7F00", # Vivid orange
  
  #"RG Progenitors"     = "#7FBFFF", # Light Blue
  "BG"         = "#8A2BE2", # Electric violet
  "BG Progenitors" = "#8A2BE2", # Approved RL progenitor label; retain BG legacy color
  "Maturing PCs" = "#1E90FF",

  # --- Reviewed off-lineage populations retained for provenance ---
  "eCN" = "#FF4500",
  "GABA Progenitors" = "#00FFD5",
  "Immune" = "#E066FF",
  "Microglia" = "#E066FF", # Legacy stored-label alias
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
  "BG Progenitors", "Maturing PCs", "eCN", "GABA Progenitors", "Microglia", "Endothelial"
)

rl_compact_lineage_markers <- list(
  "RL" = c("MKI67", "RSPO3", "OTX2"),
  "UBC" = c("EOMES", "CACNA2D1", "TRPC3"),
  "Granule" = c(
    "SFRP1", "HEY1", "NTF3", "CNTN2", "GAP43", "GABRA6"
  )
)

rl_other_markers <- list(
  "BG Progenitors" = c("SLC1A3", "SOX2", "BCAN"),
  "Maturing PCs" = c("CALB1", "FOXP4", "DAB1"),
  "eCN" = c("SLC17A6", "MEIS2", "NTNG2"),
  "GABA Progenitors" = c("PAX2", "GAD1", "GAD2"),
  "Microglia" = c("P2RY12", "SLCO2B1", "OLFML3"),
  "Endothelial" = c("CLDN5", "PECAM1", "KDR")
)

validate_compact_marker_panel(
  rl_compact_lineage_markers, c("RL", "UBC", "Granule"),
  "RL compact lineage panel", 3L, 6L
)
validate_compact_marker_panel(
  rl_other_markers, rl_other_subclusters,
  "RL Other panel", 2L, 3L
)
if (!setequal(c(rl_lineage_subclusters, rl_other_subclusters), rl_subcluster_order)) {
  stop("RL lineage and Other DotPlot groups do not cover the full RL taxonomy.")
}

# One production panel contract for merged and per-sample RL DotPlots. Keep
# raw review rl_markers unchanged so saved review signatures remain valid.
rl_production_marker_panel <- function(identities) {
  identities <- unique(as.character(identities))
  if (!length(identities) || anyNA(identities) || any(!nzchar(identities)) ||
      any(!identities %in% rl_subcluster_order)) {
    stop("RL production marker panel requires nonblank registered identities.")
  }
  broad_sets <- list(
    RL = rl_markers$RL, UBC = rl_markers$UBC,
    Granule = rl_compact_lineage_markers$Granule
  )
  broad_present <- c(
    RL = "RL SVZ" %in% identities,
    UBC = any(c("Immature UBCs", "Mature UBCs") %in% identities),
    Granule = any(c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs",
                    "Migrating GCs", "Mature GCs") %in% identities)
  )
  c(broad_sets[broad_present], rl_other_markers[intersect(rl_other_subclusters, identities)])
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
  "VZPs", "Neuroblasts", "Cycling VZPs",
  "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "GABA Interneurons", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "Cycling Glia", "BG", "BG Progenitors", "Astrocytes/Ependyma",
  "Pre-OPCs", "OPCs",
  "RL VZ", "RL SVZ", "Immature UBCs", "Mature UBCs", "eCN",
  "Prolif GCPs", "Maturing GCPs",
  "Differentiating GCs", "Migrating GCs", "Mature GCs",
  "Cycling Cells", "Pericytes", "Meninges", "Microglia", "Endothelial",
  "Immune", "Unresolved", "Unknown"
)
