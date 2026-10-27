#!/bin/bash
#SBATCH --job-name=subset_purk_gran
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/scripts/subset_purk_gran_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/scripts/subset_purk_gran_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Stop if any command fails
set -e

# Load R (adjust for your cluster)
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# Run the R script
Rscript --max-ppsize=500000 "/home/acflint/R/Projects/XeniumFCProject/scripts/subset_purk+gran.R"

echo "Subsetting job completed!"
