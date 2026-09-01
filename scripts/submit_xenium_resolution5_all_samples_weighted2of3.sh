#!/bin/bash
# Validate and submit the complete 34-sample resolution-5 weighted 2-of-3
# Xenium workflow. Run from the project root on an allocated compute node.

set -euo pipefail

project_root="/home/acflint/R/Projects/XeniumFCProject"
cd "${project_root}" || exit 1

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

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --all-samples-res5 --list

for task_id in $(seq 1 "${manifest_count}"); do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --all-samples-res5 --dry-run "${task_id}"
done

cluster_submit=$(sbatch --parsable \
  scripts/run_xenium_preprocess_03i_resolution5_all_samples.slurm)
cluster_job=${cluster_submit%%;*}

transfer_jobs=()
for reference in Aldinger Sepp Science; do
  transfer_submit=$(sbatch --parsable \
    --dependency="afterok:${cluster_job}" \
    --array="1-${manifest_count}%6" \
    --job-name="Xen_res5_all_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
  transfer_jobs+=("${transfer_submit%%;*}")
done
transfer_dependency=$(IFS=:; echo "${transfer_jobs[*]}")

build_submit=$(sbatch --parsable \
  --dependency="afterok:${transfer_dependency}" \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
build_job=${build_submit%%;*}

apply_submit=$(sbatch --parsable \
  --dependency="afterok:${build_job}" \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
apply_job=${apply_submit%%;*}

pages_submit=$(sbatch --parsable \
  --dependency="afterok:${apply_job}" \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
pages_job=${pages_submit%%;*}

merge_submit=$(sbatch --parsable \
  --dependency="afterok:${pages_job}" \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
merge_job=${merge_submit%%;*}

printf '\nSubmitted resolution-5 workflow:\n'
printf 'Clustering:       %s\n' "${cluster_job}"
printf 'Label transfers:  %s\n' "${transfer_dependency}"
printf 'Consensus build:  %s\n' "${build_job}"
printf 'Consensus/plots:  %s\n' "${apply_job}"
printf 'Report pages:     %s\n' "${pages_job}"
printf 'Report merge:     %s\n\n' "${merge_job}"

squeue -j "${cluster_job},${transfer_dependency//:/,},${build_job},${apply_job},${pages_job},${merge_job}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"
