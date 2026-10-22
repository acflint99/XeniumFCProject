#!/bin/bash
#SBATCH --job-name=VZ_PC_DotPlot
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_PC_DotPlot_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_PC_DotPlot_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
# The project .Rprofile activates the existing renv environment.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
if [[ "${PURKINJE_DOTPLOT_OVERWRITE:-false}" == true ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_diagnostic_purkinje_dotplot.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_vz_diagnostic_purkinje_dotplot.R "${driver_options[@]}"
