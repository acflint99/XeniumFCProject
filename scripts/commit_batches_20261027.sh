#!/usr/bin/env bash
# Run from the Cheaha project root. Preview by default; never submits jobs.
# Expected runtime: seconds to 2 minutes, excluding Git hooks and diff reading.
set -euo pipefail
mode=preview
batch=all
show_diff=no
while (($#)); do
  case "$1" in
    --commit) mode=commit; shift ;;
    --batch) batch="${2:?Supply batch 1..7}"; shift 2 ;;
    --diff) show_diff=yes; shift ;;
    --help) echo 'Usage: bash scripts/commit_batches_20261027.sh [--diff] [--commit] [--batch 1..7]'; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
case "$batch" in all|[1-7]) ;; *) echo 'Batch must be 1..7.' >&2; exit 1 ;; esac
cd "$(git rev-parse --show-toplevel)"
[[ -d scripts && -d config ]] || { echo 'Expected project root with scripts/ and config/.' >&2; exit 1; }
if [[ "$mode" == commit ]]; then
  git diff --cached --quiet || { echo 'Index already contains staged changes; nothing changed.' >&2; exit 1; }
  [[ -z "$(git ls-files -u)" ]] || { echo 'Resolve merge conflicts first.' >&2; exit 1; }
fi
run_batch() {
  local id="$1" title="$2"
  shift 2
  [[ "$batch" == all || "$batch" == "$id" ]] || return 0
  local path
  local pending=()
  printf '\nBatch %s: %s\n' "$id" "$title"
  for path in "$@"; do
    # Literal explicit files only: no directories or recursive staging.
    if [[ -d "$path" || -L "$path" ]]; then
      echo "Expected a regular source file, found directory/symlink: $path" >&2
      exit 1
    fi
    if [[ -n "$(git status --porcelain --untracked-files=all -- "$path")" ]]; then
      pending+=("$path")
    elif [[ ! -f "$path" ]] && ! git ls-files --error-unmatch -- "$path" >/dev/null 2>&1; then
      echo "Required batch file is absent: $path" >&2
      exit 1
    fi
  done
  if ((${#pending[@]} == 0)); then
    echo 'No remaining changes.'
    return 0
  fi
  git status --short -- "${pending[@]}"
  if [[ "$show_diff" == yes ]]; then
    git --no-pager diff HEAD -- "${pending[@]}"
    for path in "${pending[@]}"; do
      if ! git ls-files --error-unmatch -- "$path" >/dev/null 2>&1; then
        git --no-pager diff --no-index -- /dev/null "$path" || [[ $? == 1 ]]
      fi
    done
  fi
  if [[ "$mode" == commit ]]; then
    git diff --cached --quiet || { echo 'Unexpected staged changes; stopping.' >&2; exit 1; }
    git diff --check -- "${pending[@]}"
    git add -A -- "${pending[@]}"
    git diff --cached --check
    git diff --cached --stat
    git commit -m "$title"
  fi
}

run_batch 1 'Harmonize Vascular broad labels across references and downstream reports' \
  config/config.yml \
  config/reference_broad_taxonomy.csv \
  scripts/R/corrected_broad_labels.R \
  scripts/R/science_reviewed_mapping.R \
  scripts/R/subcluster_review_markers.R \
  scripts/aldinger_02_standardize_cell_types.R \
  scripts/science_01_standardize_cell_types.R \
  scripts/science_03_map_reviewed_res2.R \
  scripts/sepp_01_standardize_cell_types.R \
  scripts/reviews/Science_res2_reviewed_labels_20261022.csv \
  scripts/color_palette.R \
  scripts/xenium_rl_plot_sample.R \
  scripts/xenium_vz_rl_optional_cluster_counts.R \
  scripts/test_vascular_labels.R \
  scripts/VASCULAR_LABEL_UPDATE_20261022.md

run_batch 2 'Process retained residual cells after QC and bind downstream provenance' \
  scripts/R/xenium_residual_qc_removal.R \
  scripts/R/residual_post_qc_processing.R \
  scripts/xenium_residual_04_review_subclusters.R \
  scripts/xenium_residual_05_reintegrate_post_qc.R \
  scripts/xenium_residual_05_apply_reviewed_labels.R \
  scripts/xenium_residual_05_plot_subclusters.R \
  scripts/xenium_residual_06_map_to_samples.R \
  scripts/xenium_vz_rl_01_combine_labels.R \
  scripts/run_xenium_residual_05_reintegrate_post_qc.slurm \
  scripts/test_residual_post_qc_processing.R \
  scripts/test_xenium_residual_qc_removal.R \
  scripts/RESIDUAL_POST_QC_PROCESSING.md

run_batch 3 'Extend regional resolution diagnostics to residual post-QC inputs' \
  scripts/R/regional_resolution_inputs.R \
  scripts/xenium_regional_05c_resolution_stability.R \
  scripts/xenium_residual_05c_resolution_stability.R \
  scripts/run_xenium_residual_05c_resolution_stability.slurm \
  scripts/REGIONAL_RESOLUTION_REVIEW.md

run_batch 4 'Require source-bound RL label approval and update selected UMAP handling' \
  scripts/R/rl_label_approval.R \
  scripts/R/rl_umap_selected.R \
  scripts/xenium_rl_06_annotate_subclusters.R \
  scripts/test_rl_label_approval.R \
  scripts/REVIEW_CORRECTIONS_20261022.md

run_batch 5 'Update original-reference similarity and MetaNeighbor comparisons' \
  scripts/R/metaneighbor_original_io.R \
  scripts/R/metaneighbor_results.R \
  scripts/R/reference_similarity_plots.R \
  scripts/xenium_compare_metaneighbor.R \
  scripts/xenium_compare_reference_similarity.R \
  scripts/xenium_redraw_reference_similarity.R \
  scripts/run_xenium_compare_metaneighbor.slurm \
  scripts/run_xenium_compare_reference_similarity.slurm \
  scripts/XENIUM_METANEIGHBOR.md \
  scripts/XENIUM_REFERENCE_SIMILARITY.md

run_batch 6 'Wire post-QC review and reference validation into the refresh workflow' \
  scripts/submit_updated_references_to_merges.sh \
  scripts/test_workflow_dependencies.py \
  scripts/UPDATED_REFERENCES_TO_MERGES.md \
  scripts/WORKFLOW_ALIGNMENT_20261022.md \
  scripts/README.md

run_batch 7 'Record reference preparation launchers and dated recovery utilities' \
  scripts/run_aldinger_03_subset_gene_panel.slurm \
  scripts/run_science_01_standardize_cell_types.slurm \
  scripts/run_science_02_subset_gene_panel.slurm \
  scripts/run_sepp_02_subset_gene_panel.slurm \
  scripts/complete_aldinger_vascular_recovery_20261026.sh \
  scripts/complete_sepp_vascular_recovery_20261026.sh \
  scripts/recover_aldinger_vascular_relabel_20261026.sh \
  scripts/recover_science_annotation_oom_20261026.py \
  scripts/recover_science_panel_cpu_20261026.py \
  scripts/recover_sepp_vascular_relabel_20261026.sh

printf '\nRemaining changes (includes deliberately deferred files):\n'
git status --short
