# 0. Environment & Libraries
# Set threads to 1 to prevent workers from fighting over the 44 CPUs
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
rm(list = ls())

library(here)
library(Seurat)
library(future)
library(future.apply)

# 1. Source functions and config
source(here("scripts", "AnchorBasedTransfer_RPCA_SCT.R"))
source(here("scripts", "color_palette.R"))

# 2. Define Samples
samples_to_run <- c(
  "GZFB5_X_G" #, "GZFB_12_X_G_1", "GZFB_12_X_G_2"
  #"GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G"
  #"FB330_1_X_G", "FB78_X_G", "GZFB5_X_G", "GZFB_12_X_G_1"
  #"GZFB_12_X_G_2", "GZFB_12_X_G_3", "GZFB_12_X_G_4", "GZFB_12_X_G_5"
  #"GZFB_1_X_G", "GZFB_9_X_G_1", "GZFB_9_X_G_2", "GZFB_9_X_G_3"
)

# 3. Config
reference_name <- "Aldinger" 
output_root    <- "outputs"
log_dir        <- here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir)

# 4. Parallel Setup (Tuned for 440GB RAM)
# Using 10 workers means each worker gets ~44GB of RAM
workers <- 6 
plan(multisession, workers = workers)

# Set global size limit to 100GB to allow large SCT objects to pass between workers
options(future.globals.maxSize = 150 * 1024^3) 

# 5. Wrapper Function
run_parallel_abt <- function(sample) {
  # Force a null device inside the worker to prevent Rplots.pdf
  pdf(NULL)
  null_dev <- dev.cur()
  
  log_file <- file.path(log_dir, paste0(sample, "_ABT_parallel.log"))
  con <- file(log_file, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  
  on.exit({
    while (!is.null(dev.list())) dev.off()
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  
  message(">>> Starting Parallel Task: ", sample)
  
  result <- tryCatch({
    # Calling your run_label_transfer_sct function
    run_label_transfer_sct(
      sample_name = sample, 
      reference_name = reference_name,
      output_root = output_root
    )
    TRUE
  }, error = function(e) {
    message("!!! Error in ", sample, ": ", e$message)
    FALSE
  })
  
  return(result)
}

# 6. Execution
message(paste("Launching", length(samples_to_run), "samples across", workers, "workers..."))

results <- future_lapply(
  samples_to_run, 
  run_parallel_abt, 
  future.seed = TRUE,
  future.packages = c("Seurat", "ggplot2", "dplyr", "here")
)

plan(sequential)
message("All parallel tasks complete.")