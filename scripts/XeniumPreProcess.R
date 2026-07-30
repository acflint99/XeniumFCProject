Sys.setenv(
  OMP_NUM_THREADS = 8, #1 for workers >1 for parallelism
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1
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
source(here("scripts", "XeniumCropCerebellum.R"))
source(here("scripts", "XeniumQC.R"))
source(here("scripts", "XeniumNormCluster_res1.5.R"))

# ----------------------------
# 2. Sample names
# ----------------------------
sample_names <- c(
  "GZFB4_X_G"
)

# ----------------------------
# 3. Parallel setup
# ----------------------------
plan(sequential)
# workers <- 3  # Safe for 754GB node
# plan(multisession, workers = workers) #use for parallelism
options(future.globals.maxSize = 200 * 1024^3)  # 200GB cap
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
  
  log_file <- file.path(log_dir, paste0(sample, ".log"))
  
  # Redirect stdout (cat) to log file
  sink(log_file, type = "output", split = FALSE)
  
  start_time <- Sys.time()
  message("========================================")
  message("Starting sample: ", sample)
  message("Time: ", start_time)
  message("========================================")
  
  result <- tryCatch({
    
    xenium_cereb      <- XeniumCropCerebellum(sample)
    xenium_cereb_QC   <- qc_xenium(xenium_cereb, sample_name = sample)
    process_xenium_clusters(xenium_cereb_QC,
                            sample_name = sample)
    
    rm(xenium_cereb, xenium_cereb_QC)
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

message("All samples complete.")