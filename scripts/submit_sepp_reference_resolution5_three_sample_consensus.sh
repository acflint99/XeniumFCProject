#!/bin/bash
# After sepp_02 has been rerun with the reviewed reference filtering, propagate
# that reference through only three resolution-5 production samples. Reuse the
# existing Aldinger and Science transfers; replace only the selected Sepp
# transfers, consensus tables, consensus objects, and per-sample plots.

set -Eeuo pipefail

dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
  shift
fi
if [[ "$#" -ne 0 ]]; then
  echo "Usage: bash scripts/submit_sepp_reference_resolution5_three_sample_consensus.sh [--dry-run]" >&2
  exit 2
fi

project_root="/home/acflint/R/Projects/XeniumFCProject"
manifest_path="config/samples.csv"
email="acflint@uab.edu"
target_samples=(
  "GZFB_9_X_G_1"
  "GZFB_20_X_G_1"
  "GZFB_22_X_G_3"
)

cd "${project_root}" || exit 1

if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  exit 2
fi

manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 rows in ${manifest_path}; found ${manifest_count}." >&2
  exit 2
fi

task_ids=()
for sample_id in "${target_samples[@]}"; do
  mapfile -t matches < <(
    awk -F, -v sample_id="${sample_id}" \
      'NR > 1 && $1 == sample_id {print NR - 1}' "${manifest_path}"
  )
  if [[ "${#matches[@]}" -ne 1 ]]; then
    echo "Expected exactly one manifest match for ${sample_id}; found ${#matches[@]}." >&2
    exit 2
  fi
  task_ids+=("${matches[0]}")
done

task_ids_csv="$(IFS=,; echo "${task_ids[*]}")"
export SELECTED_SAMPLE_IDS
SELECTED_SAMPLE_IDS="$(IFS=,; echo "${target_samples[*]}")"

printf 'Selected samples: %s\n' "${SELECTED_SAMPLE_IDS}"
printf 'Manifest task IDs: %s\n' "${task_ids_csv}"

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

Rscript scripts/validate_config.R

for task_id in "${task_ids[@]}"; do
  Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
    --all-samples-res5 --dry-run Sepp "${task_id}"
done

Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res5 --selected-samples --weighted-2of3 --list

if [[ "${dry_run}" == "true" ]]; then
  printf '\nDRY-RUN SUCCESS: inputs and the three-sample task mapping passed; Slurm jobs submitted: 0.\n'
  exit 0
fi

sepp_submit="$(sbatch --parsable \
  --array="${task_ids_csv}%3" \
  --job-name="Xen_res5_3sample_Sepp" \
  --mail-user="${email}" \
  --mail-type=FAIL \
  --export=ALL,REFERENCE=Sepp,RES5_ALL_ABT_OVERWRITE=true \
  scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)"
sepp_job="${sepp_submit%%;*}"

consensus_submit="$(sbatch --parsable \
  --dependency="afterok:${sepp_job}" \
  --mail-user="${email}" \
  --mail-type=FAIL \
  --export=ALL,RES5_SELECTED_CONSENSUS_OVERWRITE=true \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_selected_samples_weighted2of3.slurm)"
consensus_job="${consensus_submit%%;*}"

apply_submit="$(sbatch --parsable \
  --array="${task_ids_csv}%3" \
  --dependency="afterok:${consensus_job}" \
  --job-name="Xen_res5_3sample_w2of3" \
  --mail-user="${email}" \
  --mail-type=END,FAIL \
  --export=ALL,RES5_ALL_W2OF3_OVERWRITE=true \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)"
apply_job="${apply_submit%%;*}"

printf '\nSubmitted targeted Sepp downstream workflow:\n'
printf 'Sepp transfer array: %s (tasks %s)\n' "${sepp_job}" "${task_ids_csv}"
printf 'Consensus build:     %s (three tables only)\n' "${consensus_job}"
printf 'Consensus apply:     %s (tasks %s)\n' "${apply_job}" "${task_ids_csv}"
printf 'Completion email:    %s\n\n' "${email}"

squeue -j "${sepp_job},${consensus_job},${apply_job}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"
