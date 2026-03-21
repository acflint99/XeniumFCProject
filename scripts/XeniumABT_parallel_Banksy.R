# 0. Environment & Libraries
# Set threads to 1 to prevent workers from fighting over the 44 CPUs
Sys.setenv(OMP_NUM_THREADS = 6, MKL_NUM_THREADS = 6, OPENBLAS_NUM_THREADS = 6)
rm(list = ls())

library(here)
library(Seurat)
library(future)
library(future.apply)

# 1. Source functions and config
source(here("scripts", "AnchorBasedTransfer_RPCA_Banksy.R"))
source(here("scripts", "color_palette.R"))

# 2. Define Samples
samples_to_run <- c(
  "GZFB_12_X_G_1",
  "FB330_1_X_G",
  "FB78_X_G"
  #"GZFB5_X_G" , "GZFB_12_X_G_1", "GZFB_12_X_G_2"
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

# 4. Parallel Setup (Optimized for 440GB RAM / 44 CPUs)
workers <- 6 # 6 workers @ ~73GB each is safer for ABT than 8 workers
plan(multisession, workers = workers)

# Global size limit for passing anchors/objects
options(future.globals.maxSize = 160 * 1024^3) 

# 5. Updated Wrapper Function
run_parallel_abt <- function(sample) {
  # 1. Immediate Graphics Silence
  options(bitmapType='cairo')
  pdf(NULL)
  
  log_file <- file.path(log_dir, paste0(sample, "_ABT_parallel.log"))
  con <- file(log_file, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  
  # Ensure cleanup happens even if the script crashes
  on.exit({
    while (!is.null(dev.list())) dev.off()
    sink(type = "message")
    sink(type = "output")
    close(con)
    gc() # Force GC after worker finishes
  }, add = TRUE)
  
  message(">>> Starting ABT for: ", sample)
  
  result <- tryCatch({
    run_label_transfer_lognorm(
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