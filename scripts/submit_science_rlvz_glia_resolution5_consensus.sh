#!/bin/bash
# Rebuild the Science reference after RLVZ -> Glia, then rerun only the
# resolution-5 all-sample weighted 2-of-3 workflow through the final consensus
# report. Existing outputs are replaced explicitly; clustering and regional
# analyses are intentionally not submitted.

set -euo pipefail

if [[ "$#" -gt 1 ]] || { [[ "$#" -eq 1 ]] && [[ "$1" != "--dry-run" ]]; }; then
  echo "Usage: bash scripts/submit_science_rlvz_glia_resolution5_consensus.sh [--dry-run]" >&2
  exit 2
fi
dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
fi

project_root="/home/acflint/R/Projects/XeniumFCProject"
cd "${project_root}" || exit 1

manifest_path="config/samples.csv"
if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  exit 2
fi

manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 manifest rows; found ${manifest_count}. Review the manifest and Slurm array bounds." >&2
  exit 2
fi

# Lightweight structural and path-only checks. These do not load large RDS
# objects and do not create or replace analysis outputs.
Rscript scripts/validate_config.R
Rscript scripts/science_01_standardize_cell_types.R --dry-run
Rscript scripts/science_02_subset_gene_panel.R --dry-run
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
  --all-samples-res5 --list

for task_id in $(seq 1 "${manifest_count}"); do
  Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
    --all-samples-res5 --dry-run Science "${task_id}"
done

Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res5 --weighted-2of3 --dry-run

for task_id in $(seq 1 "${manifest_count}"); do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --all-samples-res5 --weighted-2of3 --dry-run "${task_id}"
done

Rscript scripts/xenium_annotate_03d_plot_report.R \
  --all-samples-res5 --weighted-2of3 --dry-run

if [[ "${dry_run}" == "true" ]]; then
  echo "Dry-run complete; no Slurm jobs were submitted."
  exit 0
fi

reference_standardize_submit=$(sbatch --parsable \
  --export=ALL,SCIENCE_REFERENCE_OVERWRITE=true \
  scripts/run_science_01_standardize_cell_types.slurm)
reference_standardize_job=${reference_standardize_submit%%;*}

reference_panel_submit=$(sbatch --parsable \
  --dependency="afterok:${reference_standardize_job}" \
  --export=ALL,SCIENCE_REFERENCE_OVERWRITE=true \
  scripts/run_science_02_subset_gene_panel.slurm)
reference_panel_job=${reference_panel_submit%%;*}

science_transfer_submit=$(sbatch --parsable \
  --dependency="afterok:${reference_panel_job}" \
  --job-name="Xen_res5_all_Science" \
  --export=ALL,REFERENCE=Science,RES5_ALL_ABT_OVERWRITE=true \
  scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
science_transfer_job=${science_transfer_submit%%;*}

consensus_build_submit=$(sbatch --parsable \
  --dependency="afterok:${science_transfer_job}" \
  --export=ALL,RES5_ALL_CONSENSUS_OVERWRITE=true \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_build_job=${consensus_build_submit%%;*}

consensus_apply_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_build_job}" \
  --export=ALL,RES5_ALL_W2OF3_OVERWRITE=true \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_apply_job=${consensus_apply_submit%%;*}

report_pages_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_apply_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
report_pages_job=${report_pages_submit%%;*}

report_merge_submit=$(sbatch --parsable \
  --dependency="afterok:${report_pages_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
report_merge_job=${report_merge_submit%%;*}

printf '\nSubmitted Science RLVZ -> Glia rerun through the consensus report:\n'
printf 'Science standardization: %s\n' "${reference_standardize_job}"
printf 'Science panel subset:    %s\n' "${reference_panel_job}"
printf 'Science label transfer:  %s\n' "${science_transfer_job}"
printf 'Consensus table build:   %s\n' "${consensus_build_job}"
printf 'Consensus label apply:   %s\n' "${consensus_apply_job}"
printf 'Consensus report pages:  %s\n' "${report_pages_job}"
printf 'Consensus report merge:  %s\n\n' "${report_merge_job}"

squeue -j "${reference_standardize_job},${reference_panel_job},${science_transfer_job},${consensus_build_job},${consensus_apply_job},${report_pages_job},${report_merge_job}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"
