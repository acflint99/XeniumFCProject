#!/bin/bash
#SBATCH --job-name=RL_PostQC_Diag
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_PostQC_Diag_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_PostQC_Diag_%j.err
#SBATCH --time=12:00:00
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
if [[ "${RL_POST_QC_DIAGNOSTIC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_rl_05a_diagnose_post_qc.R "${driver_options[@]}"
