#!/bin/bash
#SBATCH --job-name=RLVZ_MarkerPlots
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RLVZ_MarkerPlots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RLVZ_MarkerPlots_%A_%a.err
#SBATCH --array=1-2%2
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Task 1 plots RL; task 2 plots VZ. The project .Rprofile activates renv.
branches=(RL VZ)
branch="${branches[$((SLURM_ARRAY_TASK_ID - 1))]}"

driver_options=("--branch=${branch}")
if [[ "${RL_VZ_MARKER_PLOTS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_vz_rl_plot_rl_vz_marker_dotplot.R "${driver_options[@]}"
