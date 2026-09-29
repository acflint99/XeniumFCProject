#!/bin/bash
# Validate or explicitly replace the existing resolution-5 production
# annotation with full-reference transfer and 0.05 margin abstention. Existing
# resolution-5 clustering objects are reused; preprocessing is not rerun.

set -euo pipefail

if [[ "$#" -ne 1 || ( "$1" != "--dry-run" && "$1" != "--overwrite" ) ]]; then
  echo "Usage: bash scripts/submit_xenium_resolution5_all_samples_full_reference_margin005.sh [--dry-run|--overwrite]" >&2
  exit 2
fi

mode="$1"
project_root="/home/acflint/R/Projects/XeniumFCProject"
cd "${project_root}" || exit 1

module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

manifest_path="config/samples.csv"
if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  exit 2
fi
manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 manifest rows; found ${manifest_count}." >&2
  exit 2
fi

Rscript scripts/validate_config.R --check-files
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
  --all-samples-res5 --list

for reference in Aldinger Sepp Science; do
  for task_id in $(seq 1 "${manifest_count}"); do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --all-samples-res5 --dry-run "${reference}" "${task_id}"
  done
done

Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res5 --weighted-2of3 --list
Rscript scripts/xenium_annotate_03_apply_consensus.R \
  --all-samples-res5 --weighted-2of3 --list
Rscript scripts/xenium_annotate_03d_plot_report.R \
  --all-samples-res5 --weighted-2of3 --list

if [[ "${mode}" == "--dry-run" ]]; then
  echo "DRY-RUN PASS: 34 samples, 102 full-reference transfers, and downstream mappings are valid; no jobs submitted."
  exit 0
fi

transfer_jobs=()
for reference in Aldinger Sepp Science; do
  transfer_submit=$(sbatch --parsable \
    --array="1-${manifest_count}%3" \
    --job-name="Xen_res5_full_${reference}" \
    --export=ALL,REFERENCE="${reference}",RES5_ALL_ABT_OVERWRITE=true \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
  transfer_jobs+=("${transfer_submit%%;*}")
done
transfer_dependency=$(IFS=:; echo "${transfer_jobs[*]}")

build_submit=$(sbatch --parsable \
  --dependency="afterok:${transfer_dependency}" \
  --export=ALL,RES5_ALL_CONSENSUS_OVERWRITE=true \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
build_job=${build_submit%%;*}

apply_submit=$(sbatch --parsable \
  --dependency="afterok:${build_job}" \
  --export=ALL,RES5_ALL_W2OF3_OVERWRITE=true \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
apply_job=${apply_submit%%;*}

pages_submit=$(sbatch --parsable \
  --dependency="afterok:${apply_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
pages_job=${pages_submit%%;*}

merge_submit=$(sbatch --parsable \
  --dependency="afterok:${pages_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
merge_job=${merge_submit%%;*}

printf '\nSubmitted full-reference resolution-5 production rerun:\n'
printf 'Label transfers:  %s\n' "${transfer_dependency}"
printf 'Consensus build:  %s\n' "${build_job}"
printf 'Consensus/plots:  %s\n' "${apply_job}"
printf 'Report pages:     %s\n' "${pages_job}"
printf 'Report merge:     %s\n\n' "${merge_job}"

squeue -j "${transfer_dependency//:/,},${build_job},${apply_job},${pages_job},${merge_job}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"
