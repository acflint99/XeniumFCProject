#!/bin/bash
# One-time migration from the former Xenium pipeline filenames to the approved
# naming convention. Run from anywhere inside the project with:
#   bash scripts/migrate_script_names.sh

set -euo pipefail

dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
elif [[ $# -gt 0 ]]; then
  echo "Usage: bash scripts/migrate_script_names.sh [--dry-run]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
cd "${project_root}"

if [[ ! -d scripts || ! -f scripts/README.md ]]; then
  echo "Could not identify the Xenium project root." >&2
  exit 1
fi

replace_reference() {
  local old_name=$1
  local new_name=$2

  OLD_NAME="${old_name}" NEW_NAME="${new_name}" \
    find scripts -maxdepth 1 -type f \
      \( -name '*.R' -o -name '*.r' -o -name '*.sh' -o -name '*.slurm' -o -name 'README.md' \) \
      -exec perl -pi -e 's/\Q$ENV{OLD_NAME}\E/$ENV{NEW_NAME}/g' {} +
}

rename_active() {
  local old_path=$1
  local new_path=$2

  if [[ -e "${old_path}" && -e "${new_path}" ]]; then
    echo "Conflict: both ${old_path} and ${new_path} exist." >&2
    exit 1
  fi

  if [[ -e "${old_path}" ]]; then
    if [[ "${dry_run}" == true ]]; then
      echo "Would rename: ${old_path} -> ${new_path}"
    else
      if git ls-files --error-unmatch "${old_path}" >/dev/null 2>&1; then
        git mv -- "${old_path}" "${new_path}"
      else
        mv -- "${old_path}" "${new_path}"
      fi
      echo "Renamed: ${old_path} -> ${new_path}"
    fi
  elif [[ -e "${new_path}" ]]; then
    echo "Already renamed: ${new_path}"
  else
    echo "Missing both ${old_path} and ${new_path}." >&2
    exit 1
  fi

  if [[ "${dry_run}" == false ]]; then
    replace_reference "$(basename "${old_path}")" "$(basename "${new_path}")"
  fi
}

retire_file() {
  local old_path=$1
  local retired_path="scripts/OLD/$(basename "${old_path}")"

  if [[ -e "${old_path}" && -e "${retired_path}" ]]; then
    echo "Conflict: both ${old_path} and ${retired_path} exist." >&2
    exit 1
  fi

  if [[ -e "${old_path}" ]]; then
    if [[ "${dry_run}" == true ]]; then
      echo "Would retire: ${old_path} -> ${retired_path}"
    else
      mv -- "${old_path}" "${retired_path}"
      echo "Retired: ${old_path} -> ${retired_path}"
    fi
  elif [[ -e "${retired_path}" ]]; then
    echo "Already retired: ${retired_path}"
  else
    echo "Not present; skipping retired file: ${old_path}"
  fi
}

mkdir -p scripts/OLD

# Preprocessing
rename_active scripts/XeniumCropCerebellum.R scripts/xenium_preprocess_01_crop_cerebellum.R
rename_active scripts/XeniumQC.R scripts/xenium_preprocess_02_qc_cells.R
rename_active scripts/XeniumNormCluster_res1.5.R scripts/xenium_preprocess_03_normalize_cluster.R
rename_active scripts/XeniumPreProcess_split_slides.R scripts/xenium_preprocess_split_slides.R
rename_active scripts/run_XeniumPreProcess_split_slides.sh scripts/run_xenium_preprocess_split_slides.slurm
if [[ "${dry_run}" == false ]]; then
  replace_reference run_XeniumPreProcess_split_slides.slurm run_xenium_preprocess_split_slides.slurm
fi

# Reference preparation
rename_active scripts/XeniumGenePanel.R scripts/xenium_extract_gene_panel.R
rename_active scripts/AldingerPreProcess.R scripts/aldinger_01_update_seurat_object.R
rename_active scripts/AldingerPreProcess2.R scripts/aldinger_02_standardize_cell_types.R
rename_active scripts/AldingerGenePanelSubset.R scripts/aldinger_03_subset_gene_panel.R
rename_active scripts/AldingerGenePanelSubset_Plots.R scripts/aldinger_03b_plot_summary.R
rename_active scripts/Aldinger_VZ_Subset.R scripts/aldinger_04_subset_vz.R
rename_active scripts/SeppPreProcess_v2.R scripts/sepp_01_standardize_cell_types.R
rename_active scripts/SeppGenePanelSubset.R scripts/sepp_02_subset_gene_panel.R
rename_active scripts/Sepp_VZ_Subset.R scripts/sepp_03_subset_vz.R
rename_active scripts/SciencePreProcess.R scripts/science_01_standardize_cell_types.R
rename_active scripts/ScienceGenePanelSubset.R scripts/science_02_subset_gene_panel.R
rename_active scripts/run_AldingerGenePanelSubset.slurm scripts/run_aldinger_03_subset_gene_panel.slurm
rename_active scripts/run_SeppGenePanelSubset.slurm scripts/run_sepp_02_subset_gene_panel.slurm
rename_active scripts/run_ScienceGenePanelSubset.slurm scripts/run_science_02_subset_gene_panel.slurm

# Reference annotation and consensus
rename_active scripts/AnchorBasedTransfer_RPCA_res1.5.R scripts/xenium_annotate_01_label_transfer_rpca.R
rename_active scripts/Xenium_ABT_res1.5_comp.R scripts/xenium_annotate_02_build_consensus.R
rename_active scripts/Xenium_ABT_ConsensusLabels.R scripts/xenium_annotate_03_apply_consensus.R
rename_active scripts/Xenium_ABT_GlobalPlot.R scripts/xenium_annotate_03b_plot_spatial.R
rename_active scripts/Xenium_ABT_PropGraph.R scripts/xenium_annotate_03c_plot_proportions.R
rename_active scripts/run_ABT_res1.5.slurm scripts/run_xenium_annotate_01_transfer.slurm
rename_active scripts/run_XeniumConsensusLabels.slurm scripts/run_xenium_annotate_03_apply_consensus.slurm

# VZ branch
rename_active scripts/Xenium_VZ_Subset_res1.5.R scripts/xenium_vz_01_subset_samples.R
rename_active scripts/Xenium_VZ_Merge_res1.5.R scripts/xenium_vz_02_merge_samples.R
rename_active scripts/Xenium_VZ_Merge_Processing_res1.5.R scripts/xenium_vz_03_integrate.R
rename_active scripts/Xenium_VZ_Merge_QC_res1.5.R scripts/xenium_vz_04_review_qc.R
rename_active scripts/Xenium_VZ_Merge_postQC_Processing_res1.5.R scripts/xenium_vz_05_reintegrate_post_qc.R
rename_active scripts/Xenium_VZ_Analysis_res1.5.R scripts/xenium_vz_06_annotate_subclusters.R
rename_active scripts/Xenium_VZ_Mapping_res1.5.R scripts/xenium_vz_07_map_to_samples.R
rename_active scripts/XeniumVZSamplePlots_Loop_res1.5.R scripts/xenium_vz_08_plot_samples.R
rename_active scripts/Xenium_VZ_CountPlot_res1.5.R scripts/xenium_vz_09_plot_counts.R
rename_active scripts/Xenium_VZ_SamplePlots_res1.5.R scripts/xenium_vz_plot_sample.R
rename_active scripts/run_XeniumVZSubset_res1.5.slurm scripts/run_xenium_vz_01_subset.slurm
rename_active scripts/run_XeniumVZMergePostQC.slurm scripts/run_xenium_vz_05_reintegrate_post_qc.slurm
rename_active scripts/run_XeniumVZMapping.slurm scripts/run_xenium_vz_07_map_to_samples.slurm
rename_active scripts/run_XeniumVZSamplePlots.slurm scripts/run_xenium_vz_08_plot_samples.slurm

# RL branch
rename_active scripts/Xenium_RL_Subset_res1.5.R scripts/xenium_rl_01_subset_samples.R
rename_active scripts/Xenium_RL_Merge_res1.5.R scripts/xenium_rl_02_merge_samples.R
rename_active scripts/Xenium_RL_Merge_Processing_res1.5.R scripts/xenium_rl_03_integrate.R
rename_active scripts/Xenium_RL_Merge_QC_res1.5.R scripts/xenium_rl_04_review_qc.R
rename_active scripts/Xenium_RL_Merge_postQC_Processing_res1.5.R scripts/xenium_rl_05_reintegrate_post_qc.R
rename_active scripts/Xenium_RL_Analysis_res1.5.R scripts/xenium_rl_06_annotate_subclusters.R
rename_active scripts/Xenium_RL_Mapping_res1.5.R scripts/xenium_rl_07_map_to_samples.R
rename_active scripts/Xenium_RL_SamplePlots_Loop_res1.5.R scripts/xenium_rl_08_plot_samples.R
rename_active scripts/Xenium_RL_CountPlot_res1.5.R scripts/xenium_rl_09_plot_counts.R
rename_active scripts/Xenium_RL_SamplePlots_res1.5.R scripts/xenium_rl_plot_sample.R
rename_active scripts/run_XeniumRLSubset_res1.5.slurm scripts/run_xenium_rl_01_subset.slurm
rename_active scripts/run_XeniumRLMergePostQC.slurm scripts/run_xenium_rl_05_reintegrate_post_qc.slurm
rename_active scripts/run_XeniumRLMapping.slurm scripts/run_xenium_rl_07_map_to_samples.slurm
rename_active scripts/run_XeniumRLSamplePlots.slurm scripts/run_xenium_rl_08_plot_samples.slurm

# Combined VZ/RL branches
rename_active scripts/Xenium_Combine_Subclusters.R scripts/xenium_vz_rl_01_combine_labels.R
rename_active scripts/Xenium_ABT_GlobalPlot_postQC.R scripts/xenium_vz_rl_optional_consensus_spatial.R
rename_active 'scripts/Xenium_CombSubclusters_Clean&Merge.R' scripts/xenium_vz_rl_03_nonspatial_merge.R
rename_active 'scripts/Xenium_CombSubclusters_Process&Plot.R' scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R
rename_active 'scripts/Xenium_VZ&RL_CountPlot_res1.5.R' scripts/xenium_vz_rl_optional_cluster_counts.R
rename_active scripts/Xenium_VZRL_Subclusters_Spatial_Merge.R scripts/xenium_vz_rl_03_spatial_merge.R
rename_active scripts/Xenium_VZRL_Subclusters_Spatial_Merge_Processing.R scripts/xenium_vz_rl_04_spatial_integrate.R
rename_active scripts/run_XeniumCombineSubclusters.slurm scripts/run_xenium_vz_rl_01_combine_labels.slurm
rename_active scripts/run_XeniumCombinedCleanMerge.slurm scripts/run_xenium_vz_rl_03_nonspatial_merge.slurm
rename_active scripts/run_XeniumCombinedProcess.slurm scripts/run_xenium_vz_rl_04_nonspatial_process.slurm
rename_active scripts/run_XeniumCombinedPlots.slurm scripts/run_xenium_vz_rl_05_nonspatial_plots.slurm
rename_active scripts/run_XeniumVZRLSpatialMerge.slurm scripts/run_xenium_vz_rl_03_spatial_merge.slurm
rename_active scripts/run_XeniumVZRLSpatialProcess.slurm scripts/run_xenium_vz_rl_04_spatial_integrate.slurm

# Cross-study comparisons
rename_active scripts/CrossStudyClusterCorr.R scripts/cross_study_cluster_correlations.R
rename_active scripts/CrossStudyClusterMarkerComp.R scripts/cross_study_cluster_marker_comparison.R
rename_active scripts/CrossStudyModuleScoreTransfer.R scripts/cross_study_module_scores.R
rename_active scripts/XenAld_VZ_Compare.R scripts/xenium_vz_compare_aldinger.R
rename_active scripts/XenSepp_VZ_Compare.R scripts/xenium_vz_compare_sepp.R
rename_active scripts/run_MarkerComp.slurm scripts/run_cross_study_cluster_marker_comparison.slurm

# Retired one-time or superseded drivers. OLD/ is ignored by the active pipeline.
retire_file scripts/XeniumPreProcess.R
retire_file scripts/XeniumPreProcess_parallel.R
retire_file scripts/XeniumPreProcess_seq.R
retire_file scripts/XeniumABT_seq.R
retire_file scripts/SaveRawClusterPlots.R
retire_file scripts/XeniumRLSamplePlots_Loop.R
retire_file scripts/run_ABT_res1.5.sh
retire_file scripts/run_MarkerComp.sh
retire_file scripts/run_XeniumRLAnalysis.slurm
retire_file scripts/run_XeniumRLSubset_res1.5.sh
retire_file scripts/run_XeniumVZSubset_res1.5.sh
retire_file scripts/run_VZRL_Subcluster_Merge.sh

if [[ "${dry_run}" == true ]]; then
  echo "Dry-run completed; no files were changed."
else
  echo "Filename migration completed. Review git status before staging anything."
fi
