# =====================================================
# Master Cell Type Color Palette
# Cerebellar Xenium Project
# =====================================================

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
  Immune      = "#E066FF"   # neon lavender
)

# Optional: function to validate palette
validate_palette <- function(labels) {
  missing <- setdiff(labels, names(cluster_colors))
  if (length(missing) > 0) {
    stop(
      paste(
        "Missing colors for:",
        paste(missing, collapse = ", ")
      )
    )
  }
}

#lock factor order
celltype_order <- c(
  "RL", "UBC", "Granule",
  "VZ","Purkinje", "GABA", "Glia", "OPC", "Meninges",
  "Endothelial", "Immune"
)