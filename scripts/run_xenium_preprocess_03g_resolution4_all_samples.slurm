#!/bin/bash
# Resolution-4 whole-tissue clustering for every sample in the configured manifest.
# Reuses the saved resolution-1.5 SNN graph and writes only beneath
# outputs/xenium/preprocess/03g_resolution4_all_samples/.

#SBATCH --job-name=Xen_res4_all
#SBATCH --array=1-34%4
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_all_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_all_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

manifest_path="config/samples.csv"
if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  exit 2
fi
manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 rows in ${manifest_path}; found ${manifest_count}. Update the Slurm array only after reviewing the manifest." >&2
  exit 2
fi

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=(--all-samples-res4)
if [[ "${RES4_ALL_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
