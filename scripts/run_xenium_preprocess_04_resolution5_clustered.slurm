#!/bin/bash
# Production resolution-5 whole-tissue clustering for every manifest sample.

#SBATCH --job-name=Xen_res5_all
#SBATCH --array=1-41%4
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_all_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_all_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
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
if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]] || (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > manifest_count )); then
  echo "Task ${SLURM_ARRAY_TASK_ID:-unset} is outside the ${manifest_count}-sample manifest. Review --list and array bounds." >&2
  exit 2
fi

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript --max-ppsize=500000 scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --all-samples-res5 --dry-run "${SLURM_ARRAY_TASK_ID}" || exit 2

driver_options=(--all-samples-res5)
if [[ "${RES5_ALL_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_preprocess_03b_resolution2_pilot.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
