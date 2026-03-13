# Set threads to 1 for stability in Seurat/MKL
Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1
)

# Clear the environment
rm(list = ls())

library(here)
library(Seurat)
library(dplyr)

# ----------------------------
# 1. Source functions
# ----------------------------
source(here("scripts", "AnchorBasedTransfer_RPCA_seq.R"))

# ----------------------------
# 2. Sample Handling
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No sample name provided! Use: Rscript script.R SampleName")
}

# We only take the first argument as our single sample
sample_name <- args[1] 
reference_name <- "Aldinger"

message("Processing single sample: ", sample_name)

# ----------------------------
# 3. Execution (Sequential)
# ----------------------------
start_time <- Sys.time()
message("========================================")
message("Starting Analysis: ", sample_name)
message("Time: ", start_time)
message("========================================")

# Run the transfer directly (No parallel plan needed)
# Wrapping in tryCatch to ensure we see the error if it fails
result <- tryCatch({
  
  # This calls the function in AnchorBasedTransfer_RPCA_seq.R
  run_label_transfer(
    sample_name = sample_name, 
    reference_name = reference_name
  )
  
  TRUE
}, error = function(e) {
  message("CRITICAL ERROR in pipeline: ", conditionMessage(e))
  return(FALSE)
})

# ----------------------------
# 4. Final Wrap-up
# ----------------------------
end_time <- Sys.time()
if (result) {
  message("========================================")
  message("Success! Sample: ", sample_name)
  message("Time elapsed: ", round(difftime(end_time, start_time, units = "mins"), 2), " mins")
  message("========================================")
} else {
  message("Pipeline failed for sample: ", sample_name)
}