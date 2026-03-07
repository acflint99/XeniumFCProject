Sys.setenv(
  OMP_NUM_THREADS = 1, #1 for workers >1 for parallelism
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
source(here("scripts", "AnchorBasedTransfer.R"))

# ----------------------------
# 2. Sample names
# ----------------------------
sample_names <- c(
  "GZFB_12_X_G_2",
  "GZFB_12_X_G_3",
  "GZFB_12_X_G_4"
  
)

reference <- "Aldinger"

# ----------------------------
# 3. Parallel setup
# ----------------------------
workers <- 3  # Safe for 754GB node
plan(multisession, workers = workers) #use for parallelism
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
  con <- file(log_file, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  
  start_time <- Sys.time()
  message("========================================")
  message("Starting sample: ", sample)
  message("Time: ", start_time)
  message("========================================")
  
  result <- tryCatch({
    
    xenium_ABT <- run_label_transfer(sample, reference)
    
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
  future.seed = TRUE,
  future.packages = c("Seurat", "ggplot2", "dplyr", "here")
)

plan(sequential)

message("All samples complete.")