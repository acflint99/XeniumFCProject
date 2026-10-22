# Approved 2026-10-19 fine-to-broad mapping for the combined VZ/RL/residual
# cohort. The original consensus_label is retained separately for provenance.
corrected_broad_from_subcluster <- c(
  "VZPs" = "VZ",
  "Maturing PCs" = "Purkinje", "Early-born PCs" = "Purkinje",
  "Late-born PCs" = "Purkinje", "Patterning PCs" = "Purkinje",
  "GABA Progenitors" = "GABA", "GABA Interneurons" = "GABA",
  "Golgi Cells" = "GABA", "MLIs" = "GABA", "iCN" = "GABA",
  "RG Progenitors" = "Glia", "Cycling Glia" = "Glia",
  "BG" = "Glia", "BG Progenitors" = "Glia",
  "Astrocytes/Ependyma" = "Glia",
  "Pre-OPCs" = "OPC", "OPCs" = "OPC",
  "RL VZ" = "RL", "RL SVZ" = "RL",
  "Immature UBCs" = "UBC", "Mature UBCs" = "UBC",
  "eCN" = "eCN",
  "Prolif GCPs" = "Granule", "Maturing GCPs" = "Granule",
  "Differentiating GCs" = "Granule", "Migrating GCs" = "Granule",
  "Mature GCs" = "Granule",
  "Cycling Cells" = "Cycling Cells",
  "Pericytes" = "Pericytes", "Microglia" = "Immune",
  "Meninges" = "Meninges", "Endothelial" = "Endothelial", "Vascular" = "Vascular",
  "Immune" = "Immune", "Unresolved" = "Unresolved",
  "RL" = "RL", "UBC" = "UBC", "Granule" = "Granule",
  "VZ" = "VZ", "Purkinje" = "Purkinje", "GABA" = "GABA",
  "Glia" = "Glia", "OPC" = "OPC", "Unknown" = "Unknown"
)

map_corrected_broad_labels <- function(subclusters) {
  labels <- as.character(subclusters)
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
    stop("Combined subcluster labels contain missing or blank values.")
  }
  unknown_numbered <- grepl("^Unknown-[0-9]{2,}$", labels)
  lookup <- labels
  lookup[unknown_numbered] <- "Unknown"
  missing <- setdiff(unique(lookup), names(corrected_broad_from_subcluster))
  if (length(missing)) {
    stop("No approved corrected broad label for: ", paste(missing, collapse = ", "))
  }
  unname(corrected_broad_from_subcluster[lookup])
}
