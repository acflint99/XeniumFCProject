library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

qc_xenium <- function(xenium_obj, sample_name) {
  
  # Record starting cell count
  initial_cells <- ncol(xenium_obj)
  
  # -------------------------------
  # 1. Ensure correct assay
  # -------------------------------
  DefaultAssay(xenium_obj) <- "Xenium"
  
  # -------------------------------
  # 2. Calculate QC metrics
  # -------------------------------
  xenium_obj$nCount_Xenium <- colSums(GetAssayData(xenium_obj, layer = "counts"))
  xenium_obj$nFeature_Xenium <- colSums(GetAssayData(xenium_obj, layer = "counts") > 0)
  
  ## Safely calculate total control probes (prevents missing assay crashes)
  safe_pull_assay <- function(assay_name) {
    if (assay_name %in% Assays(xenium_obj)) {
      return(colSums(GetAssayData(xenium_obj, assay = assay_name, layer = "counts")))
    } else {
      return(rep(0, ncol(xenium_obj)))
    }
  }
  
  blank_counts <- safe_pull_assay("BlankCodeword")
  ctrl_codeword_counts <- safe_pull_assay("ControlCodeword")
  ctrl_probe_counts <- safe_pull_assay("ControlProbe")
  genomic_counts <- safe_pull_assay("GenomicControl")
  
  total_controls <- blank_counts + ctrl_codeword_counts + ctrl_probe_counts + genomic_counts
  
  # Calculate percentage and safely handle 0/0 division (NaN)
  xenium_obj$percent_control <- (total_controls / (xenium_obj$nCount_Xenium + total_controls)) * 100
  xenium_obj$percent_control[is.na(xenium_obj$percent_control)] <- 0
  
  # -------------------------------
  # 3. Prepare outputs folder and PDF
  # -------------------------------
  qc_plot_dir <- here("outputs", "XeniumQCPlots")
  dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf_file <- here("outputs", "XeniumQCPlots", paste0(sample_name, "_QCplots.pdf"))
  pdf(pdf_file, width = 10, height = 7)
  
  # -------------------------------
  # 4. Determine thresholds (MAD Strategy)
  # -------------------------------
  mad_multiplier <- 3
  min_counts_floor <- 20
  min_features_floor <- 10
  
  # Calculate MAD for nCount
  nCount_med <- median(xenium_obj$nCount_Xenium, na.rm = TRUE)
  nCount_mad <- mad(xenium_obj$nCount_Xenium, na.rm = TRUE)
  nCount_lower <- max(min_counts_floor, nCount_med - (mad_multiplier * nCount_mad))
  nCount_upper <- quantile(xenium_obj$nCount_Xenium, 0.99, na.rm = TRUE)
  
  # Calculate MAD for nFeature
  nFeature_med <- median(xenium_obj$nFeature_Xenium, na.rm = TRUE)
  nFeature_mad <- mad(xenium_obj$nFeature_Xenium, na.rm = TRUE)
  nFeature_lower <- max(min_features_floor, nFeature_med - (mad_multiplier * nFeature_mad))
  nFeature_upper <- quantile(xenium_obj$nFeature_Xenium, 0.99, na.rm = TRUE)
  
  ## Area and Control Thresholds
  area_lower <- 15
  area_upper <- quantile(xenium_obj$cell_area, 0.99, na.rm = TRUE)
  max_control_percent <- 5
  
  # -------------------------------
  # 4b. Store thresholds in metadata
  # ------------------------------- 
  xenium_obj@misc$QC_thresholds <- list(
    nCount_lower = nCount_lower,
    nCount_upper = nCount_upper,
    nFeature_lower = nFeature_lower,
    nFeature_upper = nFeature_upper,
    area_lower = area_lower,
    area_upper = area_upper,
    max_control_percent = max_control_percent
  )
  
  # -------------------------------
  # 5. Violin plots with thresholds
  # -------------------------------
  vln_data <- FetchData(xenium_obj, vars = c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control"))
  
  vln_long <- vln_data %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(
      cols = c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control"),
      names_to = "feature",
      values_to = "value"
    )
  
  thresholds <- tibble::tibble(
    feature = c("nCount_Xenium", "nCount_Xenium", "nFeature_Xenium", "nFeature_Xenium", "cell_area", "cell_area", "percent_control"),
    value = c(nCount_lower, nCount_upper, nFeature_lower, nFeature_upper, area_lower, area_upper, max_control_percent),
    Threshold = c("Lower", "Upper", "Lower", "Upper", "Lower", "Upper", "Upper")
  )
  
  vln <- ggplot(vln_long, aes(x = 1, y = value)) +
    geom_violin(fill = "lightgray") +
    geom_hline(data = thresholds, aes(yintercept = value, color = Threshold), linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~feature, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Lower" = "red", "Upper" = "blue")) +
    ggtitle(paste0(sample_name, " QC distributions")) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  print(vln)
  
  # -------------------------------
  # 6. Feature scatter with thresholds
  # -------------------------------
  fs <- FeatureScatter(xenium_obj, "nCount_Xenium", "nFeature_Xenium") +
    geom_hline(yintercept = nFeature_lower, linetype = "dashed", color = "red") +
    geom_hline(yintercept = nFeature_upper, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = nCount_lower, linetype = "dashed", color = "red") +
    geom_vline(xintercept = nCount_upper, linetype = "dashed", color = "blue") +
    ggtitle(paste0(sample_name, " nCount vs nFeature")) +
    theme(legend.position = "none")
  
  print(fs)
  
  # -------------------------------
  # 7. Histograms
  # -------------------------------
  hist(xenium_obj$nCount_Xenium, breaks = 50, main = paste0(sample_name, ": Transcripts per Cell"), xlab = "nCount_Xenium")
  abline(v = c(nCount_lower, nCount_upper), col = c("red", "blue"), lty = 2)
  
  hist(xenium_obj$nFeature_Xenium, breaks = 50, main = paste0(sample_name, ": Genes per Cell"), xlab = "nFeature_Xenium")
  abline(v = c(nFeature_lower, nFeature_upper), col = c("red", "blue"), lty = 2)
  
  # -------------------------------
  # 8. Flag cells outside thresholds
  # -------------------------------
  xenium_obj$lowQC <- xenium_obj$nCount_Xenium < nCount_lower | 
    xenium_obj$nFeature_Xenium < nFeature_lower |
    xenium_obj$cell_area < area_lower
  
  xenium_obj$highQC <- xenium_obj$nCount_Xenium > nCount_upper | 
    xenium_obj$nFeature_Xenium > nFeature_upper |
    xenium_obj$cell_area > area_upper |
    xenium_obj$percent_control > max_control_percent
  
  # -------------------------------
  # 9. Spatial overview of failed cells
  # -------------------------------
  if ("lowQC" %in% colnames(xenium_obj@meta.data) & "highQC" %in% colnames(xenium_obj@meta.data)) {
    coords <- GetTissueCoordinates(xenium_obj)
    df <- cbind(coords, xenium_obj@meta.data)
    df$lowQC <- as.logical(df$lowQC)
    df$highQC <- as.logical(df$highQC)
    
    failed_flag <- df$lowQC | df$highQC
    df_failed <- df[failed_flag, ]
    df_passed <- df[!failed_flag, ]
    
    gg_qc <- ggplot() +
      geom_point(data = df_passed, aes(x = x, y = y), color = "grey80", size = 0.3) +
      geom_point(data = df_failed, aes(x = x, y = y), color = "red", size = 0.3) +
      coord_fixed() +
      theme_void() +
      ggtitle(paste0(sample_name, " Spatial overview (QC failed cells in red)"))
    
    print(gg_qc)
  }
  
  # -------------------------------
  # 10. Spatial Feature Plots for Area & Controls
  # -------------------------------
  p_area <- ImageFeaturePlot(xenium_obj, features = "cell_area", max.cutoff = "q99") +
    ggtitle(paste0(sample_name, ": Spatial Distribution of Cell Area")) +
    theme(legend.position = "right")
  
  # Safe plotting fallback if control probes are entirely zero
  if (max(xenium_obj$percent_control, na.rm = TRUE) > 0) {
    p_control <- ImageFeaturePlot(xenium_obj, features = "percent_control", max.cutoff = "q99") +
      ggtitle(paste0(sample_name, ": Spatial Distribution of % Control Probes")) +
      theme(legend.position = "right")
  } else {
    p_control <- ImageFeaturePlot(xenium_obj, features = "percent_control") +
      ggtitle(paste0(sample_name, ": ZERO Control Probes Detected")) +
      theme(legend.position = "right")
  }
  
  print(p_area + p_control)
  
  dev.off()  # close PDF
  
  # -------------------------------
  # 11. Filter QC cells
  # -------------------------------
  xenium_obj <- subset(
    xenium_obj,
    subset = nCount_Xenium >= nCount_lower &
      nCount_Xenium <= nCount_upper &
      nFeature_Xenium >= nFeature_lower &
      nFeature_Xenium <= nFeature_upper &
      cell_area >= area_lower &
      cell_area <= area_upper &
      percent_control <= max_control_percent
  )
  
  final_cells <- ncol(xenium_obj)
  cells_removed <- initial_cells - final_cells
  
  # -------------------------------
  # 12. Export thresholds & cell counts to .txt
  # -------------------------------
  txt_file <- here("outputs", "XeniumQCPlots", paste0(sample_name, "_QC_thresholds.txt"))
  
  threshold_text <- c(
    paste0("QC Report for ", sample_name, ":"),
    "--------------------------------------------------",
    paste0("  Initial cells: ", initial_cells),
    paste0("  Final cells:   ", final_cells),
    paste0("  Cells removed: ", cells_removed, " (", round((cells_removed/initial_cells)*100, 2), "%)"),
    "",
    "Applied QC filters:",
    "--------------------------------------------------",
    paste0("  nCount_Xenium: lower = ", round(nCount_lower, 2), ", upper = ", round(nCount_upper, 2)),
    paste0("    (Calculated via 3 MADs. Median: ", round(nCount_med, 2), " | MAD: ", round(nCount_mad, 2), ")"),
    paste0("  nFeature_Xenium: lower = ", round(nFeature_lower, 2), ", upper = ", round(nFeature_upper, 2)),
    paste0("    (Calculated via 3 MADs. Median: ", round(nFeature_med, 2), " | MAD: ", round(nFeature_mad, 2), ")"),
    paste0("  cell_area: lower = ", area_lower, ", upper = ", round(area_upper, 2)),
    paste0("  percent_control: upper = ", max_control_percent, "%")
  )
  
  writeLines(threshold_text, con = txt_file)
  
  # -------------------------------
  # 13. Save QC'd object
  # -------------------------------
  rds_out_dir <- here("outputs", "XeniumRDS")
  dir.create(rds_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_file <- here("outputs", "XeniumRDS", paste0(sample_name, "_CB_QC.rds"))
  saveRDS(xenium_obj, file = output_file, compress = FALSE)
  
  cat("QC complete for sample:", sample_name, "\n")
  cat("Initial cells:", initial_cells, "| Final cells:", final_cells, "\n")
  cat("Saved RDS file to:", output_file, "\n")
  cat("QC plots saved to:", pdf_file, "\n")
  cat("QC thresholds and cell counts saved to:", txt_file, "\n")
  
  return(xenium_obj)
}