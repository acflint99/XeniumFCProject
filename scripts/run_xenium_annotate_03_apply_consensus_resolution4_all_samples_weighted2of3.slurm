#!/bin/bash
# Apply weighted-reference 2-of-3 consensus labels to all 34 resolution-4
# all-sample objects and generate standardized DotPlots and spatial plots.

#SBATCH --job-name=Xen_res4_all_w2of3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_all_w2of3_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_all_w2of3_%A_%a.err
#SBATCH --array=1-34%6
#SBATCH --time=00:30:00
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

driver_options=(--all-samples-res4 --weighted-2of3)
if [[ "${RES4_ALL_W2OF3_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
