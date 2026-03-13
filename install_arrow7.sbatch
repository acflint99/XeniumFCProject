#!/bin/bash
#SBATCH --job-name=install_arrow
#SBATCH --output=arrow_source.log
#SBATCH --error=arrow_source.log
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

module purge
module load rc-base
module load R/4.4.0-foss-2022b

# compression libraries
module load zstd/1.5.2-GCCcore-12.2.0
module load zlib/1.2.12-GCCcore-12.2.0
module load libdeflate/1.15-GCCcore-12.2.0

# pkg-config helps Arrow detect libraries
module load pkg-config

cd ~/R/Projects/XeniumFCProject

# remove failed installs
rm -rf renv/library/*/arrow
rm -rf renv/library/*/00LOCK*

# enable full Arrow features
export LIBARROW_MINIMAL=false
export ARROW_WITH_ZSTD=ON

export R_REPOS_OVERRIDE=https://cloud.r-project.org

Rscript -e 'renv::activate(); install.packages("arrow", repos="https://cloud.r-project.org", type="source", Ncpus=16)'