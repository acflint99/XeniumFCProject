#!/bin/bash
# Label transfer for all resolution-4 all-sample objects. Submit once for each
# REFERENCE: Aldinger, Sepp, and Science.

#SBATCH --job-name=Xen_res4_all_ABT
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --array=1-34%3
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load GCC/13.2.0
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

case "${REFERENCE:-}" in
  Aldinger|Sepp|Science) ;;
  *)
    echo "REFERENCE must be Aldinger, Sepp, or Science." >&2
    exit 2
    ;;
esac

driver_options=(--all-samples-res4)
if [[ "${RES4_ALL_ABT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_annotate_01_label_transfer_rpca.R \
  "${driver_options[@]}" "${REFERENCE}" "${SLURM_ARRAY_TASK_ID}"
