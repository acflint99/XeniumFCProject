Sys.setenv(
  OMP_NUM_THREADS = 1,  
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1
)

library(here)
library(Seurat)
library(future)
library(future.apply)

# source your functions
source(here("scripts", "XeniumCropCerebellum_duo.R")) #changed for duo
source(here("scripts", "XeniumQC.R"))
source(here("scripts", "XeniumNormCluster.R"))

# single sample
sample_names <- c("GZFB_9_X_G_1")

# sequential for single node run
plan(sequential)
options(future.globals.maxSize = 200 * 1024^3)
options(future.fork.enable = FALSE)

# logs folder
log_dir <- here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir)

# run pipeline
run_xenium_pipeline <- function(sample) {
  log_file <- file.path(log_dir, paste0(sample, ".log"))
  sink(log_file, type = "output", split = FALSE)
  on.exit(sink(type="output"), add = TRUE)
  
  start_time <- Sys.time()
  message("Starting sample: ", sample, " at ", start_time)
  
  result <- tryCatch({
    xenium_cereb <- XeniumCropCerebellum(sample)
    xenium_cereb_QC <- qc_xenium(xenium_cereb, sample_name = sample)
    process_xenium_clusters(xenium_cereb_QC, sample_name = sample)
    
    rm(xenium_cereb, xenium_cereb_QC)
    gc()
    
    message("Finished sample: ", sample, " at ", Sys.time())
    TRUE
  }, error = function(e) {
    message("ERROR processing sample: ", sample)
    message(conditionMessage(e))
    FALSE
  })
  
  return(result)
}

results <- lapply(sample_names, run_xenium_pipeline)
message("All samples complete.")