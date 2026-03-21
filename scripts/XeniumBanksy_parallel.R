Sys.setenv(
  OMP_NUM_THREADS = 12, # 1 for workers >1 for parallelism
  MKL_NUM_THREADS = 12,
  OPENBLAS_NUM_THREADS = 12
)

# Clear the environment
rm(list = ls())

library(here)
library(Seurat)
library(future)
library(future.apply)

# ----------------------------
# 1. Source functions once
# ----------------------------
# Replaced QC/Crop scripts with the new Banksy clustering script
source(here("scripts", "XeniumBanksyCluster.R"))

# ----------------------------
# 2. Sample names
# ----------------------------
sample_names <- c(
  "GZFB_12_X_G_1",
  "FB330_1_X_G",
  "FB78_X_G"
)

# ----------------------------
# 3. Parallel setup
# ----------------------------
workers <- 3  # Safe for 754GB node
plan(multisession, workers = workers) 
options(future.globals.maxSize = 140 * 1024^3)  # 120GB cap
options(future.fork.enable = FALSE)

# ----------------------------
# 4. Ensure logs folder exists
# ----------------------------
log_dir <- here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir)

# ----------------------------
# 5. Pipeline function with logging
# ----------------------------
run_xenium_pipeline <- function(sample) {
  
  log_file <- file.path(log_dir, paste0(sample, "_Banksy.log"))
  
  # Redirect stdout (cat) to log file
  sink(log_file, type = "output", split = FALSE)
  
  start_time <- Sys.time()
  message("========================================")
  message("Starting BANKSY clustering for sample: ", sample)
  message("Time: ", start_time)
  message("========================================")
  
  result <- tryCatch({
    
    # Define the path to the pre-processed RDS file
    rds_file <- file.path(here("outputs", "XeniumRDS"), paste0(sample, "_CB_QC.rds"))
    
    # Failsafe: Ensure the file actually exists before proceeding
    if (!file.exists(rds_file)) {
      stop("Input RDS file not found: ", rds_file)
    }
    
    message("Loading pre-processed data: ", rds_file)
    xenium_obj <- readRDS(rds_file)
    
    message("Initiating BANKSY processing...")
    # Assuming process_xenium_clusters is the function name inside XeniumBanksyCluster.R
    # You can adjust lambda and k_geom here if you want to experiment with spatial weights
    process_xenium_clusters(xenium_obj, 
                            sample_name = sample)
    
    # Clean up memory
    rm(xenium_obj)
    gc()
    
    end_time <- Sys.time()
    message("Finished sample: ", sample)
    message("Time elapsed: ",
            round(difftime(end_time, start_time, units = "mins"), 2),
            " minutes")
    
    TRUE
    
  }, error = function(e) {
    message("ERROR processing sample: ", sample)
    message(conditionMessage(e))
    FALSE
  })
  
  # Restore stdout
  sink(type = "output")
  
  return(result)
}

# ----------------------------
# 6. Run all samples in parallel
# ----------------------------
results <- future_lapply(
  sample_names,
  run_xenium_pipeline,
  future.seed = TRUE
)

plan(sequential)

message("All samples complete. Check the logs directory for individual sample outputs.")