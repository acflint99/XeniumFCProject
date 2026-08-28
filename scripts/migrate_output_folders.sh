#!/bin/bash
# Safely migrate existing pipeline outputs to the approved directory layout.
# Run from the project root or scripts directory. Use --dry-run first.

set -euo pipefail

dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
elif [[ $# -gt 0 ]]; then
  echo "Usage: bash scripts/migrate_output_folders.sh [--dry-run]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
output_root="${project_root}/outputs"

manifest="${project_root}/config/samples.csv"
if [[ ! -f "${manifest}" ]]; then
  manifest="${script_dir}/config/samples.csv"
fi
if [[ ! -f "${manifest}" ]]; then
  echo "Could not find config/samples.csv." >&2
  exit 1
fi
if [[ ! -d "${output_root}" ]]; then
  echo "Output root does not exist: ${output_root}" >&2
  exit 1
fi

first_column=$(head -n 1 "${manifest}" | cut -d, -f1 | tr -d '"\r')
if [[ "${first_column}" != "sample_id" ]]; then
  echo "Expected sample_id to be the first column in ${manifest}." >&2
  exit 1
fi

mapfile -t sample_ids < <(
  awk -F, 'NR > 1 {
    gsub(/\r/, "", $1)
    gsub(/^"|"$/, "", $1)
    if (length($1)) print $1
  }' "${manifest}"
)
if [[ ${#sample_ids[@]} -eq 0 ]]; then
  echo "No sample IDs found in ${manifest}." >&2
  exit 1
fi

old_crop_plots="${output_root}/XeniumCropPlots"
old_rds="${output_root}/XeniumRDS"
old_qc_reports="${output_root}/XeniumQCPlots"
old_clustered_rds="${output_root}/Xenium_Res1.5_RDS"
old_clustered_plots="${output_root}/Xenium_Res1.5_Plots"

new_crop_plots="${output_root}/xenium/preprocess/01_cropped/plots"
new_crop_rds="${output_root}/xenium/preprocess/01_cropped/rds"
new_qc_reports="${output_root}/xenium/preprocess/02_qc/reports"
new_qc_rds="${output_root}/xenium/preprocess/02_qc/rds"
new_clustered_rds="${output_root}/xenium/preprocess/03_clustered/rds"
new_clustered_plots="${output_root}/xenium/preprocess/03_clustered/plots"

old_aldinger_rds="${output_root}/AldingerRDS"
old_aldinger_plots="${output_root}/AldingerPlots"
old_aldinger_vz_rds="${output_root}/Aldinger_VZ_RDS"
old_aldinger_vz_plots="${output_root}/Aldinger_VZ_Plots"
old_sepp_rds="${output_root}/SeppRDS"
old_sepp_plots="${output_root}/SeppPlots"
old_sepp_vz_rds="${output_root}/Sepp_VZ_RDS"
old_sepp_vz_plots="${output_root}/Sepp_VZ_Plots"
old_science_rds="${output_root}/ScienceRDS"
old_science_plots="${output_root}/SciencePlots"

new_aldinger_rds="${output_root}/references/aldinger/rds"
new_aldinger_plots="${output_root}/references/aldinger/plots"
new_aldinger_vz_rds="${output_root}/references/aldinger/vz/rds"
new_aldinger_vz_plots="${output_root}/references/aldinger/vz/plots"
new_sepp_rds="${output_root}/references/sepp/rds"
new_sepp_plots="${output_root}/references/sepp/plots"
new_sepp_vz_rds="${output_root}/references/sepp/vz/rds"
new_sepp_vz_plots="${output_root}/references/sepp/vz/plots"
new_science_rds="${output_root}/references/science/rds"
new_science_plots="${output_root}/references/science/plots"

external_aldinger_updated="/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat_updated.rds"
project_aldinger_updated="${new_aldinger_rds}/Aldinger_seurat_updated.rds"

annotation_old_dirs=(
  "${output_root}/Xenium_AldingerABT_Res1.5_RDS"
  "${output_root}/Xenium_AldingerABT_Res1.5_Plots"
  "${output_root}/Xenium_AldingerABT_Res1.5_Tables"
  "${output_root}/Xenium_SeppABT_Res1.5_RDS"
  "${output_root}/Xenium_SeppABT_Res1.5_Plots"
  "${output_root}/Xenium_SeppABT_Res1.5_Tables"
  "${output_root}/Xenium_ScienceABT_Res1.5_RDS"
  "${output_root}/Xenium_ScienceABT_Res1.5_Plots"
  "${output_root}/Xenium_ScienceABT_Res1.5_Tables"
  "${output_root}/Xenium_Comp_ABT_Res1.5_Tables"
  "${output_root}/Xenium_ConsensusABT_Res1.5_RDS"
  "${output_root}/Xenium_ConsensusABT_Res1.5_Plots"
  "${output_root}/Xenium_ConsensusABT_Res1.5_GlobalPlots"
  "${output_root}/Xenium_ConsensusABT_Res1.5_PCW_Plots"
  "${output_root}/validation_through_consensus"
)
annotation_new_dirs=(
  "${output_root}/xenium/annotation/01_label_transfer/aldinger/rds"
  "${output_root}/xenium/annotation/01_label_transfer/aldinger/plots"
  "${output_root}/xenium/annotation/01_label_transfer/aldinger/tables"
  "${output_root}/xenium/annotation/01_label_transfer/sepp/rds"
  "${output_root}/xenium/annotation/01_label_transfer/sepp/plots"
  "${output_root}/xenium/annotation/01_label_transfer/sepp/tables"
  "${output_root}/xenium/annotation/01_label_transfer/science/rds"
  "${output_root}/xenium/annotation/01_label_transfer/science/plots"
  "${output_root}/xenium/annotation/01_label_transfer/science/tables"
  "${output_root}/xenium/annotation/02_consensus/tables"
  "${output_root}/xenium/annotation/03_consensus_labels/rds"
  "${output_root}/xenium/annotation/03_consensus_labels/plots/samples"
  "${output_root}/xenium/annotation/03_consensus_labels/plots/spatial"
  "${output_root}/xenium/annotation/03_consensus_labels/plots/proportions"
  "${output_root}/validation/through_consensus"
)

consensus_plot_root="${output_root}/xenium/annotation/03_consensus_labels/plots"
consensus_plot_sources=(
  "${output_root}/Xenium_ConsensusABT_Res1.5_Plots"
  "${output_root}/Xenium_ConsensusABT_Res1.5_GlobalPlots"
  "${consensus_plot_root}/samples"
  "${consensus_plot_root}/spatial"
)

branch_old_dirs=(
  "${output_root}/XenAld_VZ_Res1.5_RDS"
  "${output_root}/XenAld_VZ_Res1.5_Plots"
  "${output_root}/Xenium_AldingerABT_VZ_QC_Res1.5_RDS"
  "${output_root}/XenAld_VZ_QC_Res1.5_Plots"
  "${output_root}/XenAld_VZ_QC_Res1.5_Tables"
  "${output_root}/XenAld_VZ_postQC_Res1.5_RDS"
  "${output_root}/XenAld_VZ_postQC_Res1.5_Plots"
  "${output_root}/XenAld_VZ_Subclusters_Res1.5_RDS"
  "${output_root}/XenAld_VZ_Subclusters_Res1.5_Plots"
  "${output_root}/XenAld_VZ_Subclusters_Res1.5_Tables"
  "${output_root}/Xenium_AldingerABT_VZsubclusters_Res1.5_RDS"
  "${output_root}/Xenium_AldingerABT_VZsubclusters_Res1.5_Results"
  "${output_root}/XenAld_RL_Res1.5_RDS"
  "${output_root}/XenAld_RL_Res1.5_Plots"
  "${output_root}/Xenium_AldingerABT_VZsub_RL_QC_Res1.5_RDS"
  "${output_root}/XenAld_RL_QC_Res1.5_Plots"
  "${output_root}/XenAld_RL_QC_Res1.5_Tables"
  "${output_root}/XenAld_RL_postQC_Res1.5_RDS"
  "${output_root}/XenAld_RL_postQC_Res1.5_Plots"
  "${output_root}/XenAld_RL_postQC_Res1.5_Tables"
  "${output_root}/XenAld_RL_Subclusters_Res1.5_RDS"
  "${output_root}/XenAld_RL_Subclusters_Res1.5_Plots"
  "${output_root}/XenAld_RL_Subclusters_Res1.5_Tables"
  "${output_root}/Xenium_AldingerABT_VZ&RLsubclusters_QC_Res1.5_RDS"
  "${output_root}/Xenium_AldingerABT_RLsubcluster_Res1.5_Results"
)
branch_new_dirs=(
  "${output_root}/xenium/vz/03_integrated/rds"
  "${output_root}/xenium/vz/03_integrated/plots"
  "${output_root}/xenium/vz/04_qc/rds"
  "${output_root}/xenium/vz/04_qc/plots"
  "${output_root}/xenium/vz/04_qc/tables"
  "${output_root}/xenium/vz/05_post_qc/rds"
  "${output_root}/xenium/vz/05_post_qc/plots"
  "${output_root}/xenium/vz/06_subclusters/rds"
  "${output_root}/xenium/vz/06_subclusters/plots"
  "${output_root}/xenium/vz/06_subclusters/tables"
  "${output_root}/xenium/vz/07_mapped/rds"
  "${output_root}/xenium/vz/08_sample_reports"
  "${output_root}/xenium/rl/03_integrated/rds"
  "${output_root}/xenium/rl/03_integrated/plots"
  "${output_root}/xenium/rl/04_qc/rds"
  "${output_root}/xenium/rl/04_qc/plots"
  "${output_root}/xenium/rl/04_qc/tables"
  "${output_root}/xenium/rl/05_post_qc/rds"
  "${output_root}/xenium/rl/05_post_qc/plots"
  "${output_root}/xenium/rl/05_post_qc/tables"
  "${output_root}/xenium/rl/06_subclusters/rds"
  "${output_root}/xenium/rl/06_subclusters/plots"
  "${output_root}/xenium/rl/06_subclusters/tables"
  "${output_root}/xenium/rl/07_mapped/rds"
  "${output_root}/xenium/rl/08_sample_reports"
)

old_vz_subsets="${output_root}/XenAld_VZ_Subsets_Res1.5_RDS"
new_vz_subsets="${output_root}/xenium/vz/01_subsets/rds"
new_vz_merged_rds="${output_root}/xenium/vz/02_merged/rds"
new_vz_merged_tables="${output_root}/xenium/vz/02_merged/tables"
old_rl_subsets="${output_root}/XenAld_RL_Subsets_Res1.5_RDS"
new_rl_subsets="${output_root}/xenium/rl/01_subsets/rds"
new_rl_merged_rds="${output_root}/xenium/rl/02_merged/rds"
new_rl_merged_tables="${output_root}/xenium/rl/02_merged/tables"

combined_old_dirs=(
  "${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS"
  "${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Plots"
  "${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Tables"
  "${output_root}/Xenium_ConsensusABT_Res1.5_postQC_GlobalPlots"
  "${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_Plots"
  "${output_root}/XenAld_VZ&RL_Subclusters_Res1.5_Plots"
  "${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_Plots"
)
combined_new_dirs=(
  "${output_root}/xenium/vz_rl/01_combined_labels/rds"
  "${output_root}/xenium/vz_rl/01_combined_labels/plots/combined"
  "${output_root}/xenium/vz_rl/01_combined_labels/tables"
  "${output_root}/xenium/vz_rl/01_combined_labels/plots/consensus"
  "${output_root}/xenium/vz_rl/03_processed/plots"
  "${output_root}/xenium/vz_rl/03_processed/plots/cluster_counts"
  "${output_root}/xenium/vz_rl/spatial/02_integrated/plots"
)

old_combined_mixed="${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS"
new_combined_merged_rds="${output_root}/xenium/vz_rl/02_merged/rds"
new_combined_merged_tables="${output_root}/xenium/vz_rl/02_merged/tables"
new_combined_processed_rds="${output_root}/xenium/vz_rl/03_processed/rds"
new_combined_processed_tables="${output_root}/xenium/vz_rl/03_processed/tables"

old_spatial_mixed="${output_root}/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS"
new_spatial_merged_rds="${output_root}/xenium/vz_rl/spatial/01_merged/rds"
new_spatial_merged_tables="${output_root}/xenium/vz_rl/spatial/01_merged/tables"
new_spatial_integrated_rds="${output_root}/xenium/vz_rl/spatial/02_integrated/rds"
new_spatial_integrated_tables="${output_root}/xenium/vz_rl/spatial/02_integrated/tables"

require_one_location() {
  local old_path=$1
  local new_path=$2

  if [[ -e "${old_path}" && -e "${new_path}" ]]; then
    echo "Collision: both old and new paths exist:" >&2
    echo "- ${old_path}" >&2
    echo "- ${new_path}" >&2
    exit 1
  fi
  if [[ ! -e "${old_path}" && ! -e "${new_path}" ]]; then
    echo "Missing expected output from both locations:" >&2
    echo "- ${old_path}" >&2
    echo "- ${new_path}" >&2
    exit 1
  fi
}

require_no_collision() {
  local old_path=$1
  local new_path=$2

  if [[ -e "${old_path}" && -e "${new_path}" ]]; then
    echo "Collision: both old and new paths exist:" >&2
    echo "- ${old_path}" >&2
    echo "- ${new_path}" >&2
    exit 1
  fi
}

move_directory() {
  local old_path=$1
  local new_path=$2

  if [[ -e "${old_path}" && -e "${new_path}" ]]; then
    echo "Collision: both ${old_path} and ${new_path} exist." >&2
    exit 1
  elif [[ -d "${old_path}" ]]; then
    if [[ "${dry_run}" == true ]]; then
      echo "Would move directory: ${old_path} -> ${new_path}"
    else
      mkdir -p "$(dirname "${new_path}")"
      mv -- "${old_path}" "${new_path}"
      echo "Moved directory: ${old_path} -> ${new_path}"
    fi
  elif [[ -d "${new_path}" ]]; then
    echo "Already migrated: ${new_path}"
  else
    echo "Missing both ${old_path} and ${new_path}." >&2
    exit 1
  fi
}

move_optional_directory() {
  local old_path=$1
  local new_path=$2

  if [[ -d "${old_path}" ]]; then
    if [[ "${dry_run}" == true ]]; then
      echo "Would move directory: ${old_path} -> ${new_path}"
    else
      mkdir -p "$(dirname "${new_path}")"
      mv -- "${old_path}" "${new_path}"
      echo "Moved directory: ${old_path} -> ${new_path}"
    fi
  elif [[ -d "${new_path}" ]]; then
    echo "Already migrated: ${new_path}"
  else
    echo "Optional output directory not present; skipping: ${old_path}"
  fi
}

require_flattened_directories_no_collision() {
  local destination_dir=$1
  shift
  local source_dir source_file destination_file
  declare -A planned_destinations=()

  for source_dir in "$@"; do
    [[ -d "${source_dir}" ]] || continue
    while IFS= read -r -d '' source_file; do
      destination_file="${destination_dir}/$(basename "${source_file}")"
      if [[ -e "${destination_file}" || -n "${planned_destinations[${destination_file}]:-}" ]]; then
        echo "Collision while flattening consensus plot directories:" >&2
        echo "- source: ${source_file}" >&2
        echo "- destination: ${destination_file}" >&2
        exit 1
      fi
      planned_destinations["${destination_file}"]="${source_file}"
    done < <(find "${source_dir}" -mindepth 1 -maxdepth 1 -type f -print0)
  done
}

move_optional_directory_contents() {
  local source_dir=$1
  local destination_dir=$2
  local source_file

  [[ -d "${source_dir}" ]] || return 0
  while IFS= read -r -d '' source_file; do
    move_file "${source_file}" "${destination_dir}/$(basename "${source_file}")"
  done < <(find "${source_dir}" -mindepth 1 -maxdepth 1 -type f -print0)

  if [[ "${dry_run}" == false ]]; then
    if find "${source_dir}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
      echo "Left unexpected entries for manual review: ${source_dir}"
    else
      rmdir -- "${source_dir}"
      echo "Removed empty consensus plot subfolder: ${source_dir}"
    fi
  fi
}

move_file() {
  local old_path=$1
  local new_path=$2

  if [[ -e "${old_path}" && -e "${new_path}" ]]; then
    echo "Collision: both old and new files exist:" >&2
    echo "- ${old_path}" >&2
    echo "- ${new_path}" >&2
    exit 1
  elif [[ -e "${old_path}" ]]; then
    if [[ "${dry_run}" == true ]]; then
      echo "Would move file: ${old_path} -> ${new_path}"
    else
      mkdir -p "$(dirname "${new_path}")"
      mv -- "${old_path}" "${new_path}"
    fi
  fi
}

move_subset_tree() {
  local old_path=$1
  local new_subset_path=$2
  local merged_rds_path=$3
  local merged_table_path=$4
  local branch_label=$5

  move_optional_directory "${old_path}" "${new_subset_path}"

  local legacy_merged="${new_subset_path}/Merged"
  if [[ "${dry_run}" == true && -d "${old_path}/Merged" ]]; then
    legacy_merged="${old_path}/Merged"
  fi
  move_file \
    "${legacy_merged}/Xenium_Merged_${branch_label}Subsets.rds" \
    "${merged_rds_path}/Xenium_Merged_${branch_label}Subsets.rds"
  move_file \
    "${legacy_merged}/Xenium_Merged_${branch_label}Subsets_manifest.csv" \
    "${merged_table_path}/Xenium_Merged_${branch_label}Subsets_manifest.csv"

  if [[ "${dry_run}" == false && -d "${legacy_merged}" ]]; then
    if find "${legacy_merged}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
      echo "Left unexpected merged files in place for manual review: ${legacy_merged}"
    else
      rmdir -- "${legacy_merged}"
      echo "Removed empty legacy directory: ${legacy_merged}"
    fi
  fi
}

move_count_plots() {
  local source_dir=$1
  local destination_dir=$2
  local source_file

  [[ -d "${source_dir}" ]] || return 0
  while IFS= read -r -d '' source_file; do
    move_file "${source_file}" "${destination_dir}/$(basename "${source_file}")"
  done < <(find "${source_dir}" -maxdepth 1 -type f \
    \( -name '*ClusterCountPlot.tif' -o -name '*ClusterCountPlot.pdf' \) -print0)
}

move_mixed_directory() {
  local old_path=$1
  local new_base_path=$2
  shift 2

  move_optional_directory "${old_path}" "${new_base_path}"

  local source_dir="${new_base_path}"
  if [[ "${dry_run}" == true && -d "${old_path}" ]]; then
    source_dir="${old_path}"
  fi

  while [[ $# -gt 0 ]]; do
    local file_name=$1
    local destination_dir=$2
    move_file "${source_dir}/${file_name}" "${destination_dir}/${file_name}"
    shift 2
  done
}

copy_external_file() {
  local source_path=$1
  local destination_path=$2

  if [[ -e "${destination_path}" ]]; then
    echo "Project copy already exists: ${destination_path}"
  elif [[ "${dry_run}" == true ]]; then
    echo "Would copy external file: ${source_path} -> ${destination_path}"
  else
    mkdir -p "$(dirname "${destination_path}")"
    cp -p -- "${source_path}" "${destination_path}"
    echo "Copied external file: ${source_path} -> ${destination_path}"
  fi
}

# Preflight every expected output before changing the filesystem.
reference_old_dirs=(
  "${old_aldinger_rds}" "${old_aldinger_plots}"
  "${old_aldinger_vz_rds}" "${old_aldinger_vz_plots}"
  "${old_sepp_rds}" "${old_sepp_plots}"
  "${old_sepp_vz_rds}" "${old_sepp_vz_plots}"
  "${old_science_rds}" "${old_science_plots}"
)
reference_new_dirs=(
  "${new_aldinger_rds}" "${new_aldinger_plots}"
  "${new_aldinger_vz_rds}" "${new_aldinger_vz_plots}"
  "${new_sepp_rds}" "${new_sepp_plots}"
  "${new_sepp_vz_rds}" "${new_sepp_vz_plots}"
  "${new_science_rds}" "${new_science_plots}"
)
for index in "${!reference_old_dirs[@]}"; do
  require_no_collision "${reference_old_dirs[$index]}" "${reference_new_dirs[$index]}"
done
for index in "${!annotation_old_dirs[@]}"; do
  require_no_collision "${annotation_old_dirs[$index]}" "${annotation_new_dirs[$index]}"
done
require_flattened_directories_no_collision \
  "${consensus_plot_root}" "${consensus_plot_sources[@]}"
for index in "${!branch_old_dirs[@]}"; do
  require_no_collision "${branch_old_dirs[$index]}" "${branch_new_dirs[$index]}"
done
for index in "${!combined_old_dirs[@]}"; do
  require_no_collision "${combined_old_dirs[$index]}" "${combined_new_dirs[$index]}"
done
require_no_collision "${old_vz_subsets}" "${new_vz_subsets}"
require_no_collision "${old_rl_subsets}" "${new_rl_subsets}"
require_no_collision "${old_combined_mixed}" "${new_combined_merged_rds}"
require_no_collision "${old_spatial_mixed}" "${new_spatial_merged_rds}"

combined_mixed_source="${old_combined_mixed}"
[[ -d "${combined_mixed_source}" ]] || combined_mixed_source="${new_combined_merged_rds}"
require_no_collision \
  "${combined_mixed_source}/XenAld_VZRL_clean_merge_manifest.csv" \
  "${new_combined_merged_tables}/XenAld_VZRL_clean_merge_manifest.csv"
require_no_collision \
  "${combined_mixed_source}/XenAld_VZRL_clean_merge_processed.rds" \
  "${new_combined_processed_rds}/XenAld_VZRL_clean_merge_processed.rds"
require_no_collision \
  "${combined_mixed_source}/XenAld_VZRL_clean_merge_processed_manifest.csv" \
  "${new_combined_processed_tables}/XenAld_VZRL_clean_merge_processed_manifest.csv"

spatial_mixed_source="${old_spatial_mixed}"
[[ -d "${spatial_mixed_source}" ]] || spatial_mixed_source="${new_spatial_merged_rds}"
require_no_collision \
  "${spatial_mixed_source}/XenAld_VZRL_spatial_merged_manifest.csv" \
  "${new_spatial_merged_tables}/XenAld_VZRL_spatial_merged_manifest.csv"
require_no_collision \
  "${spatial_mixed_source}/XenAld_VZRL_spatial_integrated.rds" \
  "${new_spatial_integrated_rds}/XenAld_VZRL_spatial_integrated.rds"
require_no_collision \
  "${spatial_mixed_source}/XenAld_VZRL_spatial_integrated_manifest.csv" \
  "${new_spatial_integrated_tables}/XenAld_VZRL_spatial_integrated_manifest.csv"

vz_legacy_merged="${old_vz_subsets}/Merged"
[[ -d "${vz_legacy_merged}" ]] || vz_legacy_merged="${new_vz_subsets}/Merged"
rl_legacy_merged="${old_rl_subsets}/Merged"
[[ -d "${rl_legacy_merged}" ]] || rl_legacy_merged="${new_rl_subsets}/Merged"
require_no_collision \
  "${vz_legacy_merged}/Xenium_Merged_VZSubsets.rds" \
  "${new_vz_merged_rds}/Xenium_Merged_VZSubsets.rds"
require_no_collision \
  "${vz_legacy_merged}/Xenium_Merged_VZSubsets_manifest.csv" \
  "${new_vz_merged_tables}/Xenium_Merged_VZSubsets_manifest.csv"
require_no_collision \
  "${rl_legacy_merged}/Xenium_Merged_RLSubsets.rds" \
  "${new_rl_merged_rds}/Xenium_Merged_RLSubsets.rds"
require_no_collision \
  "${rl_legacy_merged}/Xenium_Merged_RLSubsets_manifest.csv" \
  "${new_rl_merged_tables}/Xenium_Merged_RLSubsets_manifest.csv"

if [[ ! -f "${external_aldinger_updated}" && ! -f "${project_aldinger_updated}" ]]; then
  echo "Missing the updated Aldinger object from both locations:" >&2
  echo "- ${external_aldinger_updated}" >&2
  echo "- ${project_aldinger_updated}" >&2
  exit 1
fi

for sample_id in "${sample_ids[@]}"; do
  require_one_location \
    "${old_crop_plots}/${sample_id}_CB_nCount_FeatPlot.tif" \
    "${new_crop_plots}/${sample_id}_CB_nCount_FeatPlot.tif"
  require_one_location \
    "${old_rds}/${sample_id}_CB.rds" \
    "${new_crop_rds}/${sample_id}_CB.rds"
  require_one_location \
    "${old_qc_reports}/${sample_id}_QCplots.pdf" \
    "${new_qc_reports}/${sample_id}_QCplots.pdf"
  require_one_location \
    "${old_qc_reports}/${sample_id}_QC_thresholds.txt" \
    "${new_qc_reports}/${sample_id}_QC_thresholds.txt"
  require_one_location \
    "${old_rds}/${sample_id}_CB_QC.rds" \
    "${new_qc_rds}/${sample_id}_CB_QC.rds"
  require_one_location \
    "${old_clustered_rds}/${sample_id}_CB_QC_cluster.rds" \
    "${new_clustered_rds}/${sample_id}_CB_QC_cluster.rds"

  for suffix in \
    _UMAP.tif \
    _RawCluster_UMAP.tif \
    _GlobalRawClustersSpatialPlot.tif \
    _FacetRawClustersSpatialPlot.tif
  do
    require_one_location \
      "${old_clustered_plots}/${sample_id}${suffix}" \
      "${new_clustered_plots}/${sample_id}${suffix}"
  done
done

# Reference folders contain one study and artifact class, so move them intact.
for index in "${!reference_old_dirs[@]}"; do
  move_optional_directory "${reference_old_dirs[$index]}" "${reference_new_dirs[$index]}"
done
copy_external_file "${external_aldinger_updated}" "${project_aldinger_updated}"

# Label-transfer, consensus, and validation outputs have one unambiguous new
# destination each. Missing optional plot/report directories are acceptable.
for index in "${!annotation_old_dirs[@]}"; do
  move_optional_directory "${annotation_old_dirs[$index]}" "${annotation_new_dirs[$index]}"
done
for consensus_plot_source in "${consensus_plot_sources[@]}"; do
  move_optional_directory_contents "${consensus_plot_source}" "${consensus_plot_root}"
done

# The VZ and RL branches now use identical numbered stages. Their legacy
# subset folders also held merged objects, so split those two artifact classes.
move_subset_tree \
  "${old_vz_subsets}" "${new_vz_subsets}" \
  "${new_vz_merged_rds}" "${new_vz_merged_tables}" "VZ"
move_subset_tree \
  "${old_rl_subsets}" "${new_rl_subsets}" \
  "${new_rl_merged_rds}" "${new_rl_merged_tables}" "RL"

for index in "${!branch_old_dirs[@]}"; do
  move_optional_directory "${branch_old_dirs[$index]}" "${branch_new_dirs[$index]}"
done

# Cluster-count plots were historically mixed into stage 06 plot folders.
move_count_plots \
  "${output_root}/xenium/vz/06_subclusters/plots" \
  "${output_root}/xenium/vz/09_cluster_counts/plots"
move_count_plots \
  "${output_root}/xenium/rl/06_subclusters/plots" \
  "${output_root}/xenium/rl/09_cluster_counts/plots"

# Combined outputs use the same stage-first structure. Two legacy RDS folders
# mixed manifests and processed objects, so relocate only those known files.
for index in "${!combined_old_dirs[@]}"; do
  move_optional_directory "${combined_old_dirs[$index]}" "${combined_new_dirs[$index]}"
done
move_mixed_directory \
  "${old_combined_mixed}" "${new_combined_merged_rds}" \
  "XenAld_VZRL_clean_merge_manifest.csv" "${new_combined_merged_tables}" \
  "XenAld_VZRL_clean_merge_processed.rds" "${new_combined_processed_rds}" \
  "XenAld_VZRL_clean_merge_processed_manifest.csv" "${new_combined_processed_tables}"
move_mixed_directory \
  "${old_spatial_mixed}" "${new_spatial_merged_rds}" \
  "XenAld_VZRL_spatial_merged_manifest.csv" "${new_spatial_merged_tables}" \
  "XenAld_VZRL_spatial_integrated.rds" "${new_spatial_integrated_rds}" \
  "XenAld_VZRL_spatial_integrated_manifest.csv" "${new_spatial_integrated_tables}"

# Folders containing one artifact class can be moved atomically.
move_directory "${old_crop_plots}" "${new_crop_plots}"
move_directory "${old_qc_reports}" "${new_qc_reports}"
move_directory "${old_clustered_rds}" "${new_clustered_rds}"
move_directory "${old_clustered_plots}" "${new_clustered_plots}"

# XeniumRDS historically mixed cropped and QC-filtered objects, so split only
# the manifest-defined files and leave any unexpected files untouched.
for sample_id in "${sample_ids[@]}"; do
  move_file "${old_rds}/${sample_id}_CB.rds" "${new_crop_rds}/${sample_id}_CB.rds"
  move_file "${old_rds}/${sample_id}_CB_QC.rds" "${new_qc_rds}/${sample_id}_CB_QC.rds"
done

if [[ "${dry_run}" == false && -d "${old_rds}" ]]; then
  if find "${old_rds}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    echo "Left unexpected files in place for manual review: ${old_rds}"
  else
    rmdir -- "${old_rds}"
    echo "Removed empty legacy directory: ${old_rds}"
  fi
fi

if [[ "${dry_run}" == true ]]; then
  echo "Dry-run completed for ${#sample_ids[@]} samples; no files were changed."
else
  echo "Approved output-folder migration completed for ${#sample_ids[@]} samples."
fi
