# ===============================
# Function: QC Xenium object with annotated thresholds PDF
# ===============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

qc_xenium <- function(xenium_obj, sample_name) {
  
  # -------------------------------
  # 1. Ensure correct assay
  # -------------------------------
  DefaultAssay(xenium_obj) <- "Xenium"
  
  # -------------------------------
  # 2. Calculate QC metrics
  # -------------------------------
  xenium_obj$nCount_Xenium <- colSums(GetAssayData(xenium_obj, layer = "counts"))
  xenium_obj$nFeature_Xenium <- colSums(GetAssayData(xenium_obj, layer = "counts") > 0)
  
  # -------------------------------
  # 3. Prepare outputs folder and PDF
  # -------------------------------
  pdf_file <- here("outputs", paste0(sample_name, "_QCplots.pdf"))
  pdf(pdf_file, width = 10, height = 7)
  
  # -------------------------------
  # 4. Determine thresholds
  # -------------------------------
  nCount_lower <- 50
  nCount_upper <- quantile(xenium_obj$nCount_Xenium, 0.99)
  nFeature_lower <- 50
  nFeature_upper <- quantile(xenium_obj$nFeature_Xenium, 0.99)
  
  # -------------------------------
  # 4b. Report thresholds
  # -------------------------------
  cat(paste0("Applied QC filters for ", sample_name, ":\n"))
  cat(paste0("  nCount_Xenium: lower = ", nCount_lower, ", upper = ", nCount_upper, "\n"))
  cat(paste0("  nFeature_Xenium: lower = ", nFeature_lower, ", upper = ", nFeature_upper, "\n"))
  
  # -------------------------------
  # 4c. Store thresholds in metadata
  # ------------------------------- 
  xenium_obj@misc$QC_thresholds <- list(
    nCount_lower = nCount_lower,
    nCount_upper = nCount_upper,
    nFeature_lower = nFeature_lower,
    nFeature_upper = nFeature_upper
  )
  
  # -------------------------------
  # 5. Violin plots with thresholds
  # -------------------------------
  # Extract data used in VlnPlot
  vln_data <- FetchData(xenium_obj, vars = c("nCount_Xenium", "nFeature_Xenium"))
  
  # Convert to long format for ggplot
  vln_long <- vln_data %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(
      cols = c("nCount_Xenium", "nFeature_Xenium"),
      names_to = "feature",
      values_to = "value"
    )
  
  # Thresholds in long format
  thresholds <- tibble::tibble(
    feature = c("nCount_Xenium", "nCount_Xenium", "nFeature_Xenium", "nFeature_Xenium"),
    value = c(nCount_lower, nCount_upper, nFeature_lower, nFeature_upper),
    Threshold = c("Lower", "Upper", "Lower", "Upper")
  )
  
  # Plot using ggplot2 directly
  vln <- ggplot(vln_long, aes(x = 1, y = value)) +
    geom_violin(fill = "lightgray") +
    geom_hline(data = thresholds, aes(yintercept = value, color = Threshold), linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~feature, scales = "free_y") +
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
  # 7. Histograms with thresholds
  # -------------------------------
  hist(
    xenium_obj$nCount_Xenium,
    breaks = 50,
    main = paste0(sample_name, ": Transcripts per Cell"),
    xlab = "nCount_Xenium"
  )
  abline(v = c(nCount_lower, nCount_upper), col = c("red", "blue"), lty = 2)
  
  hist(
    xenium_obj$nFeature_Xenium,
    breaks = 50,
    main = paste0(sample_name, ": Genes per Cell"),
    xlab = "nFeature_Xenium"
  )
  abline(v = c(nFeature_lower, nFeature_upper), col = c("red", "blue"), lty = 2)
  
  # -------------------------------
  # 8. Flag cells outside thresholds (for spatial plotting)
  # -------------------------------
  xenium_obj$lowQC <- xenium_obj$nCount_Xenium < nCount_lower | xenium_obj$nFeature_Xenium < nFeature_lower
  xenium_obj$highQC <- xenium_obj$nCount_Xenium > nCount_upper | xenium_obj$nFeature_Xenium > nFeature_upper

  
  # -------------------------------
  # 10. Spatial QC plots
  # -------------------------------
  if ("lowQC" %in% colnames(xenium_obj@meta.data) & "highQC" %in% colnames(xenium_obj@meta.data)) {
    # Get coordinates
    coords <- GetTissueCoordinates(xenium_obj)
    df <- cbind(coords, xenium_obj@meta.data)
    
    # Make sure QC columns are logical
    df$lowQC <- as.logical(df$lowQC)
    df$highQC <- as.logical(df$highQC)
    
    # Identify failed QC cells
    failed_flag <- df$lowQC | df$highQC
    
    # Split for plotting so red points are on top
    df_failed <- df[failed_flag, ]
    df_passed <- df[!failed_flag, ]
    
    # Plot
    gg_qc <- ggplot() +
      geom_point(data = df_passed, aes(x = x, y = y), color = "grey80", size = 0.3) +
      geom_point(data = df_failed, aes(x = x, y = y), color = "red", size = 0.3) +
      coord_fixed() +
      theme_void() +
      ggtitle(paste0(sample_name, " Spatial overview (QC failed cells in red)"))
    
    print(gg_qc)
    
  } else {
    print(ImageFeaturePlot(xenium_obj, features = c("nCount_Xenium", "nFeature_Xenium")))
  }
  
  dev.off()  # close PDF
  
  # -------------------------------
  # 9. Filter QC cells
  # -------------------------------
  xenium_obj <- subset(
    xenium_obj,
    subset = nCount_Xenium > nCount_lower &
      nCount_Xenium < nCount_upper &
      nFeature_Xenium > nFeature_lower &
      nFeature_Xenium < nFeature_upper
  )
  
  # -------------------------------
  # 11. Save QC'd object
  # -------------------------------
  output_file <- here("outputs", paste0(sample_name, "_CB_QC.rds"))
  saveRDS(xenium_obj, file = output_file, compress = FALSE)
  
  cat("QC complete for sample:", sample_name, "\n")
  cat("Final number of cells:", ncol(xenium_obj), "\n")
  cat("Saved RDS file to:", output_file, "\n")
  cat("QC plots saved to:", pdf_file, "\n")
  
  return(xenium_obj)
}