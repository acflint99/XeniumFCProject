#!/bin/bash
#SBATCH --job-name=arrow_source
#SBATCH --output=arrow_source.log
#SBATCH --error=arrow_source.log
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# 1. Clean environment and load the correct, modern modules
module purge
module load rc-base
module load GCCcore/13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load OpenSSL/3
module load R/4.4.0-foss-2022b

# 2. Tell Arrow to use the system dependencies we just loaded
# LIBARROW_MINIMAL=false ensures it builds with features like S3 and GCS support
export LIBARROW_BINARY=false
export LIBARROW_MINIMAL=false
export ARROW_R_DEV=TRUE

# IMPORTANT: Remove the manual PATH/LD_LIBRARY_PATH exports you had here.
# The 'module load' commands already handle this correctly for the 13.3.0 versions.

# 3. Navigate to project and prep renv
cd ~/R/Projects/XeniumFCProject

# Remove any stalled installation locks
rm -rf renv/library/linux-rhel-7.9/R-4.4/x86_64-pc-linux-gnu/00LOCK-arrow

# 4. Install Arrow
# We use Ncpus=16 to match your SBATCH request for faster compilation
Rscript -e 'install.packages("arrow", repos = "https://packagemanager.posit.co/cran/__linux__/centos7/latest", type = "source", Ncpus = 16)'

# 5. Verify and Sync renv
Rscript -e 'library(arrow); print(packageVersion("arrow")); renv::snapshot(prompt = FALSE)'
