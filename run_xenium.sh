#!/bin/bash
#SBATCH --job-name=xenium_pipeline
#SBATCH --output=logs/xenium_%j.out   # stdout/stderr merged
#SBATCH --error=logs/xenium_%j.err    # optional separate error log
#SBATCH --nodes=1
#SBATCH --ntasks=1                     # 1 R process
#SBATCH --cpus-per-task=16             # number of cores for BLAS / Seurat
#SBATCH --mem=200G                      # memory for 150k cells
#SBATCH --time=06:00:00                # walltime
#SBATCH --partition=medium           # adjust for your cluster

module load R/4.4.0-foss-2022b                     # or your R module
cd /home/acflint/R/Projects/XeniumFCProject

# Run your script
Rscript scripts/XeniumPreProcess_seq.R
