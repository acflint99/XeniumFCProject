#!/bin/bash
#SBATCH --job-name=RL_kRes_Sensitivity
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_kRes_Sensitivity_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_kRes_Sensitivity_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --partition=short
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

# The project .Rprofile activates renv. This diagnostic never writes a
# production RDS and does not install or update packages.
driver_options=()
if [[ "${RL_K_RES_SENSITIVITY_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_rl_05b_k_resolution_sensitivity.R "${driver_options[@]}"
