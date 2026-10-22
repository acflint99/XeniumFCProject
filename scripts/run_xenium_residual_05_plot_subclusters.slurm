#!/bin/bash
#SBATCH --job-name=Xen_Residual_DotPlot
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
set --
if [[ "${RESIDUAL_DOTPLOT_OVERWRITE:-false}" == true ]]; then set -- --overwrite; fi
Rscript --max-ppsize=500000 scripts/xenium_residual_05_plot_subclusters.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_residual_05_plot_subclusters.R "$@"
