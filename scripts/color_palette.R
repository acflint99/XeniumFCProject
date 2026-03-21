# =====================================================
# Master Cell Type Color Palette
# Cerebellar Xenium Project
# =====================================================
library(scales) # Required for alpha()

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
  "VZ", "Purkinje", "GABA", "Glia", "OPC", "Meninges",
  "Endothelial", "Immune"
)


####VZ subcluster label colors & order####
vz_markers <- list(
  "VZP" = c("ASCL1"),
  "Purkinje"    = c("KITLG", "FOXP2", "EBF1", "EBF2", "CALB1", "TRPC3"),
  "GABAergic"   = c("NEUROG2", "THBS1", "SP9", "SOX14"),
  "Glia/OPC"    = c("SOX2", "TOP2A", "TNC", "AQP4", "FOXJ1", "OLIG1"),
  "Excitatory"  = c("SLC17A6", "EOMES")
)

# 1. Define the subcluster order (shared across scripts)
vz_subcluster_order <- c(
  "VZP", "Migrating PCs", "Immature PCs", "Differentiated PCs", "Mature PCs", 
  "GABA Progenitors", "Golgi Cells", "GABA Interneurons", "MLI", 
  "RG Progenitors", "Prolif RG", "BG", "Astrocytes", "Ependymal Cells", "OPCs", "eCN", "UBC Progenitors"
)

# 2. Define the Master Palette
# You can use hex codes or R color names.
vz_palette <- c(
  # --- Purkinje Lineage (High Contrast Blues/Purples) ---
  "VZP"                = "#00FFD5", # Neon Teal (Original VZ) - Very Bright
  "Migrating PCs"      = "#BBFFFF", # Pale Cyan - Extremely Bright (Contrast vs Migrating)
  "Immature PCs"       = "#1E90FF", # Vivid Blue (Original PC) - Medium
  "Differentiated PCs" = "#CC66FF", # Medium Purple (Original UBC) - Medium Dark
  "Mature PCs"         = "#8A2BE2", # Indigo - Darkest (Pops against Early/Migrating)
  
  # --- GABAergic Lineage (Revised for Maximum Contrast) ---
  "GABA Progenitors"   = "#FFFF00", 
  "Golgi Cells"        = "#CCFF00", 
  "GABA Interneurons"  = "#39FF14", 
  "MLI"                = "#228B22", 
  
  # --- Glia Lineage (Pinks/Cyans) ---
  "RG Progenitors"     = "#7FBFFF", # Pastel Blue (Original Glia)
  "Prolif RG"          = "#00CED1", # Neon Cyan (Original Meninges)
  "BG"                 = "#F0B3FF", # Neon Lavender (Original Immune)
  "Astrocytes"         = "#FF00FF", # Magenta (Pure contrast vs BG/RG)
  "Ependymal Cells"    = "#E6B3FF", # 
  
  # --- Other (Warm Tones / High Contrast) ---
  "OPCs"               = "#FFD700", # Gold/Neon Yellow
  "eCN"                = "#FF4500", # Orange Red (Moved from UBC)
  "UBC Progenitors"    = "#FF1493"  # Deep Pink (Changed to free up orange/red space)
)
vz_palette["NA"] <- alpha("gray20", 0.3)

# 3. Define a background/null color for highlight plots
bg_color <- alpha("gray20", 0.3)

####RL subcluster label colors & order####
rl_subcluster_order <- c(
  "RLVZ", "RL Transition", "RLSVZ", "Cycling Cells",
  "Cycling GCPs/EGL", "Differentiating GCPs", "Migrating GCs", "Mature GCs/IGL",
  "UBC Progenitors", "Transitioning UBCs", "Mature UBCs",
  "RG Progenitors", "Astrocytes", "Purkinje Cells"
)

rl_palette <- c(
  # --- Rhombic Lip (High-Contrast Cyans/Blues) ---
  "RLVZ"                = "#E0FFFF", # Dark Turquoise (The "Teal" Bridge)
  "RL Transition"       = "#00CED1", # Azure Blue (Transitioning out)
  "RLSVZ"               = "#0080FF", # Medium Blue (Dark & Rich anchor - better than Navy)
  
  # --- The Transition (The "Standout" Bridge) ---
  "Cycling Cells"       = "#FFFF00", # Pure Yellow (Back in, but bracketed by Blue/Green)
  
  "Cycling GCPs/EGL"    = "#ADFF2F", # Neon Green
  "Differentiating GCPs"= "#00FF7F", # Spring Green
  "Migrating GCs"       = "#32CD32", # Deep Emerald
  "Mature GCs/IGL"      = "#006400", # Green-Yellow (Bright pop to distinguish from RL)
  
  "UBC Progenitors"    = "#FFB6C1", # Deep Pink
  "Transitioning UBCs" = "#FF4500", # Orange Red
  "Mature UBCs"        = "#8B0000", # Dark Red
  
  "RG Progenitors"     = "#7FBFFF", # Light Blue
  "Astrocytes"         = "#FF00FF", # Magenta
  "Purkinje Cells"     = "#8A2BE2"  # Indigo (No gradient)
)

rl_markers <- list(
  "RL" = c("CALCB", "RSPO3", "EOMES", "AURKB"),
  "Granule"    = c("MCM2", "BMP2", "NTF3", "KITLG"),
  "UBC"   = c("SOX5", "ZFHX4", "TRPC3"),
  "Other"    = c("SLC1A3", "BCAN", "PRPH")
)

####PC subcluster label colors & order####
PC_subcluster_order <- c(
  "1 - Progenitor PCs", "5 - Progenitor PCs", "0 - Early PCs", "3 - Early PCs",
  "6 - Transition PCs", "8 - Transition PCs", "2 - Mature PCs", "4 - Mature PCs",
  "7 - GCPs"
)

PC_palette <- c(
  # --- Rhombic Lip (High-Contrast Cyans/Blues) ---
  "1 - Progenitor PCs"  = "#E0FFFF", # Dark Turquoise (The "Teal" Bridge)
  "5 - Progenitor PCs"  = "#00CED1", # Azure Blue (Transitioning out)
  "0 - Early PCs"       = "#0080FF", # Medium Blue (Dark & Rich anchor - better than Navy)
  "3 - Early PCs"       = "#00FFD5", # Pure Yellow (Back in, but bracketed by Blue/Green)
  
  "6 - Transition PCs"    = "#FF4500", # Neon Green
  "8 - Transition PCs"   = "#8B0000", # Spring Green
  "2 - Mature PCs"       = "#FF00FF", # Deep Emerald
  "4 - Mature PCs"      = "#8A2BE2", # Green-Yellow (Bright pop to distinguish from RL)
  
  "7 - GCPs"           = "#00FF7F" # Deep Pink

)

PC_markers <- list(
  "Progenitor" = c("NES", "SOX9"),
  "Early"    = c("POSTN", "PTGER3", "CNTNAP4", "FOXP1", "FOXP4"),
  "Transition"   = c("GRIK3", "ZEB2", "SP9", "MAB21L1"),
  "Mature"    = c("EBF2", "RORB", "GAD1", "EN1", "EN2"),
  "Luo" = c("ITPR1", "VSTM2L", "NEFL", "PCDH10")
)

####VZ subcluster label colors & order post BANKSY####
vz_B_markers <- list(
  "Purkinje"    = c("FOXP4", "FOXP2", "BCL11A", "EBF1", "EBF2", "EN1", "CALB1", "HPCA"),
  "GABAergic"   = c("NEUROG2", "SOX14", "DMBX1", "MEIS2", "SLC17A6"),
  "Glia/OPC"    = c("BCAN", "TNC", "SLC1A3","AQP4", "FOXJ1", "OLIG1"),
  "Cycling"  = c("MKI67", "AURKB")
)

# 1. Define the subcluster order (shared across scripts)
vz_B_subcluster_order <- c(
  "PC Progenitors", "Maturing PCs", "Migratory PCs", "Intermediate PCs", "Rostral PCs", "Mature PCs", 
  "GABAergic Progenitors", "GABAergic DCN", "Glutamatergic DCN", 
  "Glial Progenitors", "Bergmann Glia", "Astrocytes", "Ependymal Cells", "OPCs", "Cycling Cells"
)

# 2. Define the Master Palette
# You can use hex codes or R color names.
vz_B_palette <- c(
  # --- Purkinje Lineage (High Contrast Blues/Purples) ---
  "PC Progenitors"      = "#00FFD5", # Pale Cyan - Extremely Bright (Contrast vs Migrating)
  "Maturing PCs" = "#BBFFFF",
  "Migratory PCs"       = "#00CED1", # Vivid Blue (Original PC) - Medium
  "Intermediate PCs" = "#1E90FF", # Medium Purple (Original UBC) - Medium Dark
  "Rostral PCs"         = "#CC66FF", # Indigo - Darkest (Pops against Early/Migrating)
  "Mature PCs"  ="#E6B3FF",
  
  # --- GABAergic Lineage (Revised for Maximum Contrast) ---
  "GABAergic Progenitors"   = "#FFFF00", 
  "GABAergic DCN"        = "#CCFF00", 
  "Glutamatergic DCN"                = "#228B22", 
  
  # --- Glia Lineage (Pinks/Cyans) ---
  "Glial Progenitors" = "#8A2BE2",
  "Bergmann Glia"                 = "#F0B3FF", # Neon Lavender (Original Immune)
  "Astrocytes"         = "#FF00FF", # Magenta (Pure contrast vs BG/RG)
  "Ependymal Cells"    = "#FF4500", # 
  
  # --- Other (Warm Tones / High Contrast) ---
  "OPCs"               = "#8B0000", # Gold/Neon Yellow
  "Cycling Cells" = "#FFD700"
)
vz_B_palette["NA"] <- alpha("gray20", 0.3)
