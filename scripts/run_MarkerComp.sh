#!/bin/bash
#SBATCH --job-name=MarkerComp_R_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/MarkerComp_R_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/MarkerComp_R_job_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=short

# --------------------------
# Load R module
# --------------------------
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# --------------------------
# Move into project directory
# --------------------------
cd ~/R/Projects/XeniumFCProject

# --------------------------
# Disable renv sandbox to avoid long delays
# --------------------------
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

# --------------------------
# Restore project environment (only installs missing packages)
# --------------------------
Rscript -e '
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", repos = "https://cloud.r-project.org/")
}
library(renv)
renv::restore(prompt = FALSE)
'

# --------------------------
# Run your R script using the project library
# --------------------------
Rscript scripts/CrossStudyClusterMarkerComp.R
