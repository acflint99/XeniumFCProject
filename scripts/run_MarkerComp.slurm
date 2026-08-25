#!/bin/bash
#SBATCH --job-name=MarkerComp_R_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/MarkerComp_R_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/MarkerComp_R_job_%j.err
#SBATCH --time=01:00:00
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
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

# --------------------------
# Run your R script using the project library
# --------------------------
Rscript scripts/CrossStudyClusterMarkerComp.R
