#!/bin/bash
#SBATCH --job-name=RL_UMAP_Compare
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_UMAP_Compare_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_UMAP_Compare_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R
