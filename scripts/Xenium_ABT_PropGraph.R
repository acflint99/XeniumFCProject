rm(list = ls())

library(Seurat)
library(purrr)
library(readxl)
library(here)

# 1. Load Master Palette & Settings
source(here("scripts", "color_palette.R")) 

# SPECIFY CLUSTERS TO INCLUDE: 
# Edit this list to include only the clusters you want to see in the plot
target_clusters <- c("RL", "UBC", "Granule", "Purkinje", "GABA") 

# 2. Load Sample Mapping
sample_map <- read_excel(here("metadata", "samples_meta.xlsx")) 

# 1. Create a new directory for the updated files
out_dir <- here("outputs", "Xenium_AldingerABT_Res1.5_PCW_RDS")
if (!dir.exists(out_dir)) dir.create(out_dir)
plot_dir <- here("outputs", "Xenium_AldingerABT_Res1.5_PCW_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# 3. List and Process RDS Files
data_path <- here("outputs", "Xenium_AldingerABT_Res1.5_RDS")
sample_files <- list.files(path = data_path, pattern = "\\.rds$", full.names = TRUE)

all_data <- map_dfr(sample_files, function(f) {
  fname <- basename(f)
  s_id <- gsub(".rds", "", fname)
  
  # IMPROVED REGEX: 
  # 1. Remove the Aldinger suffix
  temp_id <- gsub("_Aldinger.*$", "", s_id)
  
  # 2. Remove trailing _1, _2, etc. (e.g., GZFB_12_X_G_1 -> GZFB_12_X_G)
  # This only removes an underscore + digit IF it's at the very end
  base_id <- gsub("_\\d+$", "", temp_id)
  
  # Map PCW from Excel
  age_string <- sample_map$PCW[match(base_id, sample_map$sample)]
  
  # Diagnostic message to help you see the match
  message("File: ", s_id, " -> Matched ID: ", base_id, " -> Age: ", age_string)
  
  # Stop if age is still NA to avoid processing the whole heavy object for nothing
  if(is.na(age_string)) {
    warning("Could not find age for: ", base_id, ". check Excel 'sample' column.")
    return(NULL)
  }
  
  obj <- readRDS(f)
  
  # --- PERSISTENT MODIFICATION ---
  # Add the PCW to the object's metadata permanently
  obj$PCW <- age_string
  
  # Overwrite the original file with the updated version
  saveRDS(obj, file = file.path(out_dir, fname))
  # -------------------------------
  
  df <- obj@meta.data %>%
    mutate(PCW_label = age_string) %>%
    select(PCW_label, cluster_weighted)
  
  rm(obj); gc()
  return(df)
})

# 4. Filter and Aggregate
plot_df <- all_data %>%
  # NEW: Only keep the clusters you listed in target_clusters
  filter(cluster_weighted %in% target_clusters) %>%
  
  # Clean PCW for numeric sorting
  mutate(PCW_num = as.numeric(gsub("PCW", "", PCW_label))) %>%
  
  # Calculate proportions based ONLY on the filtered set
  group_by(PCW_num, cluster_weighted) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(PCW_num) %>%
  mutate(relative_percent = (count / sum(count)) * 100) %>%
  ungroup()

# 5. Set Factor Levels for Plotting
plot_df$PCW_num <- factor(plot_df$PCW_num, levels = sort(unique(plot_df$PCW_num)))

# This ensures the legend and stack order match your master celltype_order
plot_df$cluster_weighted <- factor(plot_df$cluster_weighted, 
                                   levels = intersect(celltype_order, target_clusters))

# 6. Plot with your custom colors
p <- ggplot(plot_df, aes(x = PCW_num, y = relative_percent, fill = cluster_weighted)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
  scale_fill_manual(values = cluster_colors) + 
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    x = "Age (PCW)", 
    y = "Relative percent (of selected clusters)", 
    fill = "Cell Type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black")
  )

Cairo::CairoTIFF(
  filename = here(plot_dir, "XenAldingerABT_res1.5_ClusterPropPlot.tif"),
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()
ggplot2::ggsave(
  filename = here(plot_dir, "XenAldingerABT_res1.5_ClusterPropPlot.pdf"),
  plot = p, device = grDevices::cairo_pdf, width = 10, height = 6
)
