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
  RL          = "#FF7F00",  # bright neon orange
  UBC         = "#9370DB",  # medium purple
  Granule     = "#FF00FF",  # bright magenta
  VZ          = "#00FFD5",  # neon teal
  Glia        = "#7FBFFF",  # pastel blue
  OPC         = "#FFD700",  # neon yellow
  Meninges    = "#00CED1",  # neon cyan / teal
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

#lock factor order
celltype_order <- c(
  "RL", "UBC", "Granule",
  "VZ", "Purkinje", "GABA", "Glia", "OPC", "Meninges",
  "Endothelial", "Immune", "Unknown"
)

# Shared broad-label DotPlot scales. Keeping these fixed makes color, dot size,
# label order, and legend orientation directly comparable across samples,
# references, voting methods, and final consensus plots.
broad_dotplot_col_min <- -2.5
broad_dotplot_col_max <- 2.5
broad_dotplot_dot_min <- 0
broad_dotplot_dot_max <- 100
broad_dotplot_dot_scale <- 6

standardize_broad_dotplot <- function(plot, label_order = celltype_order) {
  plot +
    ggplot2::scale_y_discrete(limits = rev(label_order), drop = FALSE) +
    ggplot2::scale_color_gradient(
      low = "lightgrey",
      high = "red",
      limits = c(broad_dotplot_col_min, broad_dotplot_col_max),
      breaks = c(broad_dotplot_col_min, 0, broad_dotplot_col_max),
      oob = scales::squish,
      name = "Average expression\n(scaled)"
    ) +
    ggplot2::scale_radius(
      range = c(0, broad_dotplot_dot_scale),
      limits = c(broad_dotplot_dot_min, broad_dotplot_dot_max),
      breaks = c(0, 25, 50, 75, 100),
      name = "Percent expressed"
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, direction = "vertical"),
      color = ggplot2::guide_colorbar(order = 2, direction = "vertical")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      legend.direction = "vertical"
    )
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
  "VZP" = c("PRDM13", "ASCL1"),
  "Purkinje"    = c("FOXP1", "ITPR1","EBF2", "NDNF", "CALB1","TRPC3", "PCDH10", "EN1", "EBF1"),
  "GABAergic"   = c("NEUROG2", "NXPH2", "SP9", "SOX14"),
  "Glia/OPC"    = c("SOX2", "TNC", "AQP4", "FOXJ1", "OLIG1"),
  "Other"  = c("BUB1", "SLC17A6", "EOMES")
)

# 1. Define the subcluster order (shared across scripts)
vz_subcluster_order <- c(
  "VZPs", "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "BG", "Astrocytes/Ependyma",
  "OPCs", "Cycling Cells", "eCN", "GCPs"
)

# 2. Define the Master Palette
# You can use hex codes or R color names.
vz_palette <- c(
  # --- Purkinje Lineage (High Contrast Blues/Purples) ---
  "VZPs"                = "#BBFFFF", # Neon Teal (Original VZ) - Very Bright
  "Maturing PCs"       = "#1E90FF", # Vivid Blue (Original PC) - Medium
  "Early-born PCs" = "#CC66FF", # Medium Purple (Original UBC) - Medium Dark
  "Late-born PCs"         = "#00CED1", # Indigo - Darkest (Pops against Early/Migrating)
  "Patterning PCs"         = "#0047AB", # Indigo - Darkest (Pops against Early/Migrating)
  
  # --- GABAergic Lineage (Revised for Maximum Contrast) ---
  "GABA Progenitors"   = "#CCFF00", 
  "Golgi Cells"        = "#00FF7F", 
  "MLIs"  = "#32CD32", 
  "iCN"                = "#006400", 
  
  # --- Glia Lineage (Pinks/Cyans) ---
  "RG Progenitors"     = "#7FBFFF", # Pastel Blue (Original Glia)
  "BG"                 = "#8A2BE2", # Neon Lavender (Original Immune)
  "Astrocytes/Ependyma"         = "#E6B3FF", # Magenta (Pure contrast vs BG/RG)
  
  # --- Other (Warm Tones / High Contrast) ---
  "OPCs"               = "#FFD700", # Gold/Neon Yellow
  "Cycling Cells" = "#FFFF00",
  "eCN"                = "#FF4500", # Orange Red (Moved from UBC)
  "GCPs"    = "#FF1493"  # Deep Pink (Changed to free up orange/red space)
)
vz_palette["NA"] <- alpha("gray20", 0.3)

# 3. Define a background/null color for highlight plots
bg_color <- alpha("gray20", 0.3)

####RL subcluster label colors & order####
rl_subcluster_order <- c(
  "RL VZ", "RL SVZ",
  "Intermediate Progenitors", "Immature UBCs", "Mature UBCs",
  "Prolif GCPs", "Maturing GCPs",
  "Differentiating GCs", "Migrating GCs", "Mature GCs",
  "BG", "Ependymal Cells", "Maturing PCs" #"RG Progenitors"
)

rl_palette <- c(
  # --- Rhombic Lip (High-Contrast Cyans/Blues) ---
  "RL VZ"       = "#BBFFFF", # Azure Blue (Transitioning out)
  "RL SVZ"               = "#FFB6C1", # Medium Blue (Dark & Rich anchor - better than Navy)
  
  # --- The Transition (The "Standout" Bridge) ---
    "Prolif GCPs"    = "#FFD700", # Neon Green
  "Maturing GCPs"       = "#00FF7F", # Deep Emerald
  "Differentiating GCs"      = "#006400", # Green-Yellow (Bright pop to distinguish from RL)
  "Migrating GCs" = "#00CED1",
  "Mature GCs" = "#0047AB",
  
  "Intermediate Progenitors" = "#FF1493",
  "Immature UBCs" = "#FF4500", # Orange Red
  "Mature UBCs"        = "#8B0000", # Dark Red
  
  #"RG Progenitors"     = "#7FBFFF", # Light Blue
  "BG"         = "#8A2BE2", # Magenta
  "Ependymal Cells"     = "#E6B3FF",
  "Maturing PCs" = "#1E90FF"

)
rl_palette["NA"] <- alpha("gray20", 0.3)

rl_markers <- list(
  "RL" = c("HES1", "RSPO3"),
  "UBC"   = c("EOMES", "BCAN", "ZFHX4", "TRPC3"),
  "Granule"    = c("CENPA", "PRPH", "CNTN2", "NTF3", "BMP2", "GABRA6"),
  "Other"    = c("SLC1A3", "FOXJ1", "FOXP4")
)

# Merge all palettes into one master lookup for the combined column
# using 'c()' on named vectors merges them
subcluster_palette <- c(vz_palette, rl_palette, cluster_colors)

# Remove any duplicate names (e.g., if "Purkinje" is in both) 
# to keep the vector clean
subcluster_palette <- subcluster_palette[!duplicated(names(subcluster_palette))]

master_subcluster_order <- c(
  "VZPs", "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN",
  "RG Progenitors", "BG", "Astrocytes/Ependyma", "Ependymal Cells", "OPCs",
  "RL VZ", "RL SVZ",
  "Intermediate Progenitors", "Immature UBCs", "Mature UBCs", "eCN",
  "Prolif GCPs", "Maturing GCPs", "GCPs",
  "Differentiating GCs", "Migrating GCs", "Mature GCs",
  "Cycling Cells", "Meninges",  "Endothelial", "Immune", "Unknown"
)
