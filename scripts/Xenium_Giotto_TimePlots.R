# ==============================================================================
# SCRIPT 1: TEMPORAL AXIS META-ANALYSIS (USING PCW & SIGNIFICANCE FILTERING)
# ==============================================================================
# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

library(dplyr)
library(ggplot2)
library(stringr)
library(here)

output_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_GiottoBroad_RDS")

# SET THIS ONCE: "knn" or "delaunay"
target_network <- "delaunay" 

# Create a sub-folder for this specific network to keep outputs organized
plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_GiottoBroad_Plots", target_network)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. Define Temporal Metadata 
meta_df <- data.frame(
  Sample = c("GZFB_12_X_G_3", "GZFB5_X_G", "FB328_1_X_G", "FB330_1_X_G", "GZFB4_X_G",  "FB124_X_G"), #
  PCW = c("10PCW", "11PCW", "12PCW", "13PCW", "16PCW",  "17PCW") # Using strings to match the gsub logic #
)

# Ensure PCW is numeric so ggplot scales the x-axis proportionally
meta_df <- meta_df %>%
  mutate(PCW_num = as.numeric(gsub("PCW", "", PCW, ignore.case = TRUE)))

# 2. Load and Aggregate Proximity Stats (The "Master" Unfiltered Table)
prox_files <- list.files(output_dir, pattern = "_proximity_stats\\.rds$", full.names = TRUE)
prox_files <- prox_files[basename(prox_files) %in% paste0(meta_df$Sample, "_proximity_stats.rds")]

temporal_prox <- bind_rows(lapply(prox_files, function(f) {
  s_name <- str_remove(basename(f), "_proximity_stats\\.rds")
  res <- readRDS(f)
  
  if(!is.null(res)) {
    # FIX 1: Use the correct slot name (enrichm_res)
    df <- as.data.frame(res$enrichm_res) 
    
    if(nrow(df) > 0) {
      df$Sample <- s_name
      
      # FIX 2: Split "CellA--CellB" into two columns for easy filtering
      df <- df %>%
        mutate(
          cell_1 = sapply(str_split(unified_int, "--"), `[`, 1),
          cell_2 = sapply(str_split(unified_int, "--"), `[`, 2)
        )
      
      return(df)
    }
  }
  return(NULL)
})) %>% left_join(meta_df, by = "Sample")

# ==============================================================================
# 3. UPDATED: Load and Aggregate Ligand-Receptor Stats (Dual Network)
# ==============================================================================
# Update the pattern to match your new "dual" suffix
lr_files <- list.files(output_dir, pattern = "_LR_interactions_dual\\.rds$", full.names = TRUE)
lr_files <- lr_files[basename(lr_files) %in% paste0(meta_df$Sample, "_LR_interactions_dual.rds")]

temporal_lr <- bind_rows(lapply(lr_files, function(f) {
  s_name <- str_remove(basename(f), "_LR_interactions_dual\\.rds")
  res_list <- readRDS(f)
  
  if(!is.null(res_list[[target_network]])) {
    df <- as.data.frame(res_list[[target_network]])
    df$Sample <- s_name
    df$Network <- target_network # Keep track of which network this came from
    return(df)
  }
  return(NULL)
})) %>% left_join(meta_df, by = "Sample")


# ==============================================================================
# ANALYSIS A: CELL-CELL PROXIMITY (HETERO ONLY)
# ==============================================================================

# Filter for strictly significant, hetero spatial attraction
sig_proximity <- temporal_prox %>%
  filter(enrichm > 0) %>%            
  filter(p.adj_higher < 0.05) %>%    
  filter(type_int == "hetero")       # FIX: Keep only interactions between different cell types

# Find the most consistent neighbors across time
top_prox_pairs <- sig_proximity %>%
  group_by(cell_1, cell_2) %>%
  summarise(times_significant = n(), avg_score = mean(enrichm), .groups = "drop") %>% 
  arrange(desc(times_significant), desc(avg_score))

print("Top 10 Spatially Colocalized Cell Pairs (Hetero Only):")
head(top_prox_pairs, 10)

# --- PLOT: PROXIMITY GRADIENT OVER PCW ---
# Automatically grab the #1 most consistent pair to plot:
target_cell_1 <- top_prox_pairs$cell_1[1] 
target_cell_2 <- top_prox_pairs$cell_2[1]

# Query the UNFILTERED 'temporal_prox' to ensure we see the score drop to zero 
plot_prox_data <- temporal_prox %>% 
  filter(cell_1 == target_cell_1 & cell_2 == target_cell_2)

# Plot using the corrected 'enrichm' column
ggplot(plot_prox_data, aes(x = PCW_num, y = enrichm)) +
  geom_line(color = "firebrick", size = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = unique(meta_df$PCW_num)) + 
  labs(
    title = paste("Proximity Dynamics:", target_cell_1, "and", target_cell_2), 
    x = "Post-Conceptional Weeks (PCW)", 
    y = "Proximity Enrichment Score (Log2FC)"
  ) +
  theme_minimal()


# ==============================================================================
# ANALYSIS B: LIGAND-RECEPTOR SIGNALING (CORRECTED COLUMNS)
# ==============================================================================

# Filter for strong, significant communication
# Standard 'Bio-sanity' Filters for Spatial Signaling
sig_lr_data <- temporal_lr %>%
  filter(p.adj < 0.05) %>%
  filter(PI > 0) %>%
  # Ensure the Sender actually has the Ligand
  filter(lig_expr > 0.1) %>% 
  # Ensure the Receiver actually has the Receptor
  filter(rec_expr > 0.1)

# Find the strongest, most consistent signaling pathways
top_lr_pairs <- sig_lr_data %>%
  group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
  summarise(times_significant = n(), max_PI = max(PI), .groups = "drop") %>%
  arrange(desc(times_significant), desc(max_PI))

print("Top 10 Strongest Ligand-Receptor Interactions:")
head(top_lr_pairs, 10)

# --- PLOT: LR HEATMAP OVER PCW ---
# Let's take the Top 15 pathways to visualize in a temporal heatmap
top_15_targets <- head(top_lr_pairs, 15)

# Filter the UNFILTERED master table for our Top 15 targets
plot_lr_data <- temporal_lr %>%
  inner_join(top_15_targets %>% select(ligand, receptor, lig_cell_type, rec_cell_type), 
             by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
  # Create a clean label for the Y-axis: Ligand(Sender) -> Receptor(Receiver)
  mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", receptor, "(", rec_cell_type, ")"))

# Create the heatmap
ggplot(plot_lr_data, aes(x = factor(PCW_num), y = Interaction_Label, fill = PI)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(name = "Comm Score (PI)") +
  labs(
    title = "Top Signaling Pathways Over Development",
    x = "Post-Conceptional Weeks (PCW)", 
    y = "Ligand (Sender) -> Receptor (Receiver)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

# ==============================================================================
# ANALYSIS C: PURKINJE-CENTRIC SIGNALING
# ==============================================================================

# 1. Filter the master LR table for Purkinje involvement (as sender OR receiver)
# Note: Ensure "Purkinje" matches your exact cluster name (e.g., "Purkinje_cell")
purkinje_lr_clean <- temporal_lr %>%
  filter(lig_cell_type == "Purkinje" | rec_cell_type == "Purkinje") %>%
  filter(p.adj < 0.05) %>%
  filter(PI > 0) %>%
  filter(lig_expr > 0.1) %>% # Sender must actually have the ligand
  filter(rec_expr > 0.1)   # Receiver must actually have the receptor

# 2. Identify the most significant Purkinje-specific interactions
top_purkinje_pairs <- purkinje_lr_clean %>%
  group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
  summarise(times_significant = n(), max_PI = max(PI), .groups = "drop") %>%
  arrange(desc(times_significant), desc(max_PI))

print("Top Purkinje-related Signaling Interactions:")
print(head(top_purkinje_pairs, 15))

# 3. Create a Purkinje-specific Heatmap
# We'll take the top 20 interactions involving Purkinje cells
top_20_purkinje <- head(top_purkinje_pairs, 20)

plot_purkinje_data <- purkinje_lr_clean %>%
  inner_join(top_20_purkinje %>% select(ligand, receptor, lig_cell_type, rec_cell_type), 
             by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
  mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", 
                                    receptor, "(", rec_cell_type, ")"))

ggplot(plot_purkinje_data, aes(x = factor(PCW_num), y = Interaction_Label, fill = PI)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(name = "Comm Score (PI)", option = "magma") + # Using 'magma' to visually distinguish from the general plot
  labs(
    title = "Purkinje Cell Signaling Dynamics (10-17 PCW)",
    subtitle = "Interactions where Purkinje is either Sender or Receiver",
    x = "Post-Conceptional Weeks (PCW)", 
    y = "Interaction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# ==============================================================================
# ANALYSIS D: DIRECTIONAL PURKINJE HEATMAPS (SENDER VS RECEIVER)
# ==============================================================================

# 1. Separate the Cleaned Data
purkinje_as_sender <- purkinje_lr_clean %>% filter(lig_cell_type == "Purkinje")
purkinje_as_receiver <- purkinje_lr_clean %>% filter(rec_cell_type == "Purkinje")

# 2. Function to create the Label and Factor the Y-axis for better plotting
prep_plot_data <- function(df, top_n = 20) {
  top_pairs <- df %>%
    group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
    summarise(times_sig = n(), max_pi = max(PI), .groups = "drop") %>%
    arrange(desc(times_sig), desc(max_pi)) %>%
    head(top_n)
  
  df_plot <- df %>%
    inner_join(top_pairs %>% select(ligand, receptor, lig_cell_type, rec_cell_type), 
               by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
    mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", 
                                      receptor, "(", rec_cell_type, ")"))
  return(df_plot)
}

plot_sender_data <- prep_plot_data(purkinje_as_sender)
plot_receiver_data <- prep_plot_data(purkinje_as_receiver)

# --- PLOT 1: PURKINJE AS SENDER ---
p_sender <- ggplot(plot_sender_data, aes(x = factor(PCW_num), y = Interaction_Label, fill = PI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Comm Score", option = "magma") +
  labs(title = "Purkinje Cell Output (Sender)", x = "PCW", y = "Target") +
  theme_minimal()

# --- PLOT 2: PURKINJE AS RECEIVER ---
p_receiver <- ggplot(plot_receiver_data, aes(x = factor(PCW_num), y = Interaction_Label, fill = PI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Comm Score", option = "viridis") + # Different color for distinction
  labs(title = "Purkinje Cell Input (Receiver)", x = "PCW", y = "Source") +
  theme_minimal()

# Display or Save
print(p_sender)
print(p_receiver)

library(purrr)

# 1. Calculate Correlation with Time for each interaction
linear_dynamics <- purkinje_lr_clean %>%
  group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
  # We need at least 3-4 timepoints to calculate a meaningful correlation
  filter(n() >= 4) %>% 
  summarise(
    time_correlation = cor(PCW_num, PI, method = "pearson"),
    avg_PI = mean(PI),
    .groups = "drop"
  ) %>%
  # Filter for strong linear trends (|r| > 0.7 is a common threshold)
  filter(abs(time_correlation) > 0.7) %>%
  arrange(desc(time_correlation))

# 2. Separate into 'Increasing' and 'Decreasing'
increasing_pathways <- linear_dynamics %>% filter(time_correlation > 0)
decreasing_pathways <- linear_dynamics %>% filter(time_correlation < 0)

print("Top Pathways that Increase Linearly with Age:")
head(increasing_pathways, 10)

# 1. Get the specific top pairs (including their cell types)
# This ensures we only keep the "Main" conversation for each interaction
top_linear_full_info <- linear_dynamics %>%
  arrange(desc(abs(time_correlation))) %>%
  head(12) %>%
  select(ligand, receptor, lig_cell_type, rec_cell_type)

# 2. Use semi_join to filter the plotting data to ONLY these specific pairs
plot_linear_data <- purkinje_lr_clean %>%
  semi_join(top_linear_full_info, 
            by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
  mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", 
                                    receptor, "(", rec_cell_type, ")"))

# 3. Plot - Now there will only be one dot per timepoint
p1 <- ggplot(plot_linear_data, aes(x = PCW_num, y = PI, color = Interaction_Label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Interaction_Label, scales = "free_y") +
  labs(title = paste("Purkinje Linear Signaling Trajectories -", toupper(target_network)),
       x = "Post-Conceptional Weeks (PCW)", 
       y = "Comm Score (PI)") +
  theme_minimal() +
  theme(legend.position = "none") # Labels are already in facet titles

file_name_purkinje <- paste0("Purkinje_Linear_Signaling_", target_network, "_LinePlot.png")

ggsave(file.path(plot_dir, file_name_purkinje), 
       plot = p1, width = 16, height = 12, dpi = 300)

# ==============================================================================
# ANALYSIS E: GLOBAL LINEAR SIGNALING (ALL CELL TYPES)
# ==============================================================================

# 1. Calculate Correlation for ALL significant interactions
global_linear_dynamics <- sig_lr_data %>%
  group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
  # Filter for interactions present in at least 4 timepoints
  filter(n() >= 4) %>% 
  summarise(
    time_correlation = cor(PCW_num, PI, method = "pearson"),
    avg_PI = mean(PI),
    .groups = "drop"
  ) %>%
  # Focus on very strong linear trends
  filter(abs(time_correlation) > 0.8) %>%
  arrange(desc(abs(time_correlation)))

print("Top 10 Global Linear Developmental Trajectories:")
print(head(global_linear_dynamics, 10))

# 2. Select the top 12 most 'Linear' global pairs for plotting
top_global_info <- global_linear_dynamics %>%
  head(12) %>%
  select(ligand, receptor, lig_cell_type, rec_cell_type)

# 3. Filter the main table and create labels
plot_global_linear_data <- sig_lr_data %>%
  semi_join(top_global_info, 
            by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
  mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", 
                                    receptor, "(", rec_cell_type, ")"))

# 4. Generate Global Line Plot
p_global_linear <- ggplot(plot_global_linear_data, aes(x = PCW_num, y = PI, color = Interaction_Label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Interaction_Label, scales = "free_y", ncol = 4) +
  labs(
    title = paste("Global Linear Signaling Trajectories -", toupper(target_network)),
    x = "Post-Conceptional Weeks (PCW)", 
    y = "Comm Score (PI)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold")
  )

# 5. Save the plot
# Change the filename to include the network
file_name_global <- paste0("Global_Linear_Signaling_", target_network, "_LinePlot.png")

ggsave(file.path(plot_dir, file_name_global), 
       plot = p_global_linear, width = 16, height = 12, dpi = 300)

# ==============================================================================
# ANALYSIS F: GRANULE-CENTRIC LINEAR SIGNALING
# ==============================================================================

# 1. Filter master LR table for Granule involvement (Sender or Receiver)
# Ensure "Granule" matches your exact metadata label
granule_lr_clean <- temporal_lr %>%
  filter(lig_cell_type == "Granule" | rec_cell_type == "Granule") %>%
  filter(p.adj < 0.05) %>%
  filter(PI > 0) %>%         # Active signaling only
  filter(lig_expr > 0.1) %>% # Bio-sanity ligand filter
  filter(rec_expr > 0.1)     # Bio-sanity receptor filter

# 2. Calculate Correlation with Time
granule_linear_dynamics <- granule_lr_clean %>%
  group_by(ligand, receptor, lig_cell_type, rec_cell_type) %>%
  filter(n() >= 4) %>% 
  summarise(
    time_correlation = cor(PCW_num, PI, method = "pearson"),
    avg_PI = mean(PI),
    .groups = "drop"
  ) %>%
  filter(abs(time_correlation) > 0.7) %>%
  arrange(desc(abs(time_correlation)))

print("Top Granule-related Linear Signaling Interactions:")
print(head(granule_linear_dynamics, 12))

# 3. Prepare Plotting Data (Top 12 pairs)
top_granule_info <- granule_linear_dynamics %>%
  head(12) %>%
  select(ligand, receptor, lig_cell_type, rec_cell_type)

plot_granule_linear_data <- granule_lr_clean %>%
  semi_join(top_granule_info, 
            by = c("ligand", "receptor", "lig_cell_type", "rec_cell_type")) %>%
  mutate(Interaction_Label = paste0(ligand, "(", lig_cell_type, ") -> ", 
                                    receptor, "(", rec_cell_type, ")"))

# 4. Generate Line Plot
p_granule_linear <- ggplot(plot_granule_linear_data, aes(x = PCW_num, y = PI, color = Interaction_Label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Interaction_Label, scales = "free_y", ncol = 4) +
  labs(
    title = paste("Granule Linear Signaling Trajectories -", toupper(target_network)),
    x = "Post-Conceptional Weeks (PCW)", 
    y = "Comm Score (PI)"
  ) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(size = 8, face = "bold"))

# 5. Save the plot
file_name_granule <- paste0("Granule_Linear_Signaling_", target_network, "_LinePlot.png")

ggsave(file.path(plot_dir, file_name_granule), 
       plot = p_granule_linear, width = 16, height = 12, dpi = 300)