#!/bin/bash
#SBATCH --job-name=Xen_VZ_Integrate
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Integrate_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Integrate_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=600G
#SBATCH --partition=medium
#SBATCH --constraint=intel
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export XENIUM_STAGE3_CLUSTER_WORKERS="${VZ_INTEGRATE_CLUSTER_WORKERS:-1}"

driver_options=()
if [[ "${VZ_INTEGRATE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
if [[ "${VZ_INTEGRATE_PLOTS_ONLY:-false}" == "true" ]]; then
  driver_options+=(--plots-only)
fi
if [[ "${VZ_INTEGRATE_INCLUDE_UNCORRECTED:-false}" == "true" ]]; then
  driver_options+=(--include-uncorrected)
fi

restart_count="${SLURM_RESTART_COUNT:-0}"
if ! [[ "${restart_count}" =~ ^[0-9]+$ ]]; then
  echo "SLURM_RESTART_COUNT must be a non-negative integer when set." >&2
  exit 2
fi

checkpoint_dir="outputs/xenium/vz/03_integrated/rds"
precluster_checkpoint="${checkpoint_dir}/Xenium_VZ_precluster_checkpoint.rds"
postcluster_checkpoint="${checkpoint_dir}/Xenium_VZ_postcluster_checkpoint.rds"
resume_requested="${VZ_INTEGRATE_RESUME_CHECKPOINT:-false}"

if (( restart_count > 0 )) && \
   [[ "${VZ_INTEGRATE_PLOTS_ONLY:-false}" != "true" ]] && \
   [[ -s "${postcluster_checkpoint}" || -s "${precluster_checkpoint}" ]]; then
  resume_requested="true"
  echo "Slurm restart ${restart_count}: resuming the most advanced available VZ checkpoint."
fi

if [[ "${resume_requested}" == "true" ]]; then
  driver_options+=(--resume-checkpoint)
fi

Rscript --max-ppsize=500000 scripts/xenium_vz_03_integrate.R "${driver_options[@]}"
