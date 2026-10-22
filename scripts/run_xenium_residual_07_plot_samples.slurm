#!/bin/bash
#SBATCH --job-name=Xen_Residual_Plots
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Residual_Plots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Residual_Plots_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Empty positional arguments are safe under nounset in older Bash releases.
set --
if [[ "${RESIDUAL_PLOTS_OVERWRITE:-false}" == true ]]; then
  set -- --overwrite
fi
case "${RESIDUAL_PLOTS_MODE:-samples}" in
  samples)
    Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R \
      --dry-run "${SLURM_ARRAY_TASK_ID}"
    Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R \
      "$@" "${SLURM_ARRAY_TASK_ID}"
    ;;
  combine)
    if [[ "${SLURM_ARRAY_TASK_COUNT:-1}" -ne 1 ]]; then
      echo "Residual combine mode requires --array=1." >&2
      exit 2
    fi
    Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R --combine --dry-run
    Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R --combine "$@"
    ;;
  *) echo "RESIDUAL_PLOTS_MODE must be samples or combine." >&2; exit 2 ;;
esac
