Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1
)

library(here)
library(Seurat)
library(future)
library(future.apply)

# Clear the environment
rm(list = ls())

# ----------------------------
# 1. Source functions once
# ----------------------------
source(here("scripts", "XeniumCropCerebellum.R"))
source(here("scripts", "XeniumQC.R"))
source(here("scripts", "XeniumNormCluster.R"))

# ----------------------------
# 2. Sample names
# ----------------------------
sample_names <- c(
  "GZFB4_X_G",
  "GZFB_1_X_G",
  "GZFB_12_X_G_1"
  # add more later
)

# ----------------------------
# 3. Parallel setup
# ----------------------------
workers <- 5  # Safe for 754GB node
plan(multisession, workers = workers)
options(future.globals.maxSize = 200 * 1024^3)  # 200GB cap

# ----------------------------
# 4. Ensure logs folder exists
# ----------------------------
log_dir <- here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir)

# ----------------------------
# 5. Pipeline function with logging
# ----------------------------
run_xenium_pipeline <- function(sample) {
  
  log_file <- file.path(log_dir, paste0(sample, ".log"))
  
  # Open log connection
  con <- file(log_file, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  
  start_time <- Sys.time()
  message("========================================")
  message("Starting sample: ", sample)
  message("Time: ", start_time)
  message("========================================")
  
  # ---- Main workflow ----
  tryCatch({
    xenium_cereb      <- XeniumCropCerebellum(sample)
    xenium_cereb_QC   <- qc_xenium(xenium_cereb, sample_name = sample)
    process_xenium_clusters(xenium_cereb_QC,
                            sample_name = sample)
    
    rm(xenium_cereb, xenium_cereb_QC)
    gc()
    
    end_time <- Sys.time()
    message("Finished sample: ", sample)
    message("Time elapsed: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
    
    TRUE  # Return success
  }, error = function(e) {
    message("ERROR processing sample: ", sample)
    message(e)
    FALSE  # Return failure
  }, finally = {
    sink(type = "output")
    sink(type = "message")
    close(con)
  })
}

# ----------------------------
# 6. Run all samples in parallel
# ----------------------------
results <- future_lapply(
  sample_names,
  run_xenium_pipeline,
  future.stdout = TRUE
)

plan(sequential)

message("All samples complete.")