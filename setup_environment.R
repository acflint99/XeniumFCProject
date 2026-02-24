# ===============================
# Xenium Project Environment Setup (HPC-safe, home-directory aware)
# ===============================

library(renv)  # load renv first

# -------------------------------
# 0. Check working directory
# -------------------------------
project_root <- getwd()
home_dir <- Sys.getenv("HOME")
setwd("/home/acflint/R/Projects/XeniumFCProject")

if (normalizePath(project_root) == normalizePath(home_dir)) {
  stop("ERROR: You are in your home directory. Please set your working directory to your project folder before initializing renv.")
} else {
  message("Current working directory is safe: ", project_root)
}

# -------------------------------
# 1. Initialize renv for this project
# -------------------------------
if (!file.exists("renv.lock")) {
  message("Initializing renv in project folder...")
  renv::init(bare = TRUE)
} else {
  message("renv.lock already exists, skipping init")
}

# -------------------------------
# 2. HPC-safe installation of BLAS-dependent packages
# -------------------------------
hpc_packages <- c("Matrix", "KernSmooth")

for (pkg in hpc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s from source (HPC-safe)...", pkg))
    install.packages(pkg, type = "source")
  }
}

# -------------------------------
# 3. List of core packages for Xenium analysis
# -------------------------------
packages <- c(
  # Core single-cell + Xenium
  "Seurat",
  
  # Data manipulation
  "here",
  "magrittr",
  "dplyr",
  
  # Visualization
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "ggplot2",
  
  # Spatial / geometry (skip s2 to avoid HPC build issues)
  "sf",   # note: sf will still work without s2 for most Xenium analyses
  
  # File formats
  "arrow",
  "jsonlite",
  
  # Reproducibility / workflow
  "rmarkdown",
  "knitr",
  "targets"  # optional, pipeline management
)

# -------------------------------
# 4. Install remaining packages via renv
# -------------------------------
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s via renv...", pkg))
    renv::install(pkg)
  }
}

# -------------------------------
# 5. Load all libraries
# -------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(sf)  # note: works without s2 for typical Xenium workflows

# -------------------------------
# 6. Snapshot the environment
# -------------------------------
renv::snapshot()

# -------------------------------
# 7. Create outputs folder and save session info
# -------------------------------
if(!dir.exists("outputs")) dir.create("outputs")
writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")

message("HPC-safe environment setup complete. renv.lock created, packages installed and libraries loaded. s2 skipped.")