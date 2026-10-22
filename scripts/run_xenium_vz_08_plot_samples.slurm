#!/bin/bash
#SBATCH --job-name=Xen_VZ_Plots
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Plots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Plots_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

plot_mode="${VZ_PLOTS_MODE:-samples}"
driver_options=()
if [[ "${VZ_PLOTS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

case "${plot_mode}" in
  samples)
    Rscript --max-ppsize=500000 scripts/xenium_vz_08_plot_samples.R \
      "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
    ;;
  combine)
    if [[ "${SLURM_ARRAY_TASK_COUNT:-1}" -ne 1 ]]; then
      echo "VZ combine mode must be submitted with --array=1." >&2
      exit 2
    fi
    Rscript --max-ppsize=500000 scripts/xenium_vz_08_plot_samples.R \
      --combine "${driver_options[@]}"
    ;;
  *)
    echo "VZ_PLOTS_MODE must be 'samples' or 'combine'; received '${plot_mode}'." >&2
    exit 2
    ;;
esac
