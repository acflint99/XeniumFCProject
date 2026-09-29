#!/bin/bash
# Quarantine retired whole-tissue resolutions and superseded trial outputs only
# after the complete 34-sample resolution-5 production chain is present.

set -euo pipefail

mode=${1:---dry-run}
confirmation=${2:-}
case "${mode}" in
  --dry-run|--quarantine|--delete-quarantine) ;;
  *)
    echo "Usage: bash scripts/cleanup_non_resolution5_outputs.sh [--dry-run|--quarantine|--delete-quarantine DELETE_NON_RES5_QUARANTINE]" >&2
    exit 2
    ;;
esac

if [[ "${mode}" == "--delete-quarantine" &&
      "${confirmation}" != "DELETE_NON_RES5_QUARANTINE" ]]; then
  echo "Permanent deletion requires the exact confirmation token DELETE_NON_RES5_QUARANTINE." >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
output_root="${project_root}/outputs"
sample_manifest="${project_root}/config/samples.csv"
quarantine="${project_root}/quarantine_non_resolution5_2026-09-29"
retained_manifest="${output_root}/validation/resolution5_cleanup_quarantine_2026-09-29.tsv"

[[ -d "${output_root}" ]] || {
  echo "Output root not found: ${output_root}" >&2
  exit 1
}
[[ -f "${sample_manifest}" ]] || {
  echo "Sample manifest not found: ${sample_manifest}" >&2
  exit 1
}

expected_samples=$(awk -F, '
  NR > 1 {
    gsub(/\r/, "", $1)
    gsub(/^"|"$/, "", $1)
    if (length($1)) count++
  }
  END { print count + 0 }
' "${sample_manifest}")
[[ "${expected_samples}" -eq 34 ]] || {
  echo "Expected 34 biological samples; found ${expected_samples}." >&2
  exit 1
}

count_matching_files() {
  local directory=$1
  local pattern=$2
  if [[ ! -d "${directory}" ]]; then
    echo 0
    return
  fi
  find "${directory}" -maxdepth 1 -type f -name "${pattern}" -print | wc -l
}

require_count() {
  local description=$1
  local directory=$2
  local pattern=$3
  local observed
  observed=$(count_matching_files "${directory}" "${pattern}")
  printf '%-42s %s/%s\n' "${description}" "${observed}" "${expected_samples}"
  [[ "${observed}" -eq "${expected_samples}" ]] || {
    echo "STOP: incomplete resolution-5 production output set." >&2
    exit 1
  }
}

resolution5_preprocess="${output_root}/xenium/preprocess/04_resolution5_clustered"
resolution5_annotation="${output_root}/xenium/annotation"
resolution5_transfer="${resolution5_annotation}/01_label_transfer"
resolution5_consensus_tables="${resolution5_annotation}/02_consensus_weighted_2of3/tables"
resolution5_consensus="${resolution5_annotation}/03_consensus_labels_weighted_2of3"

if [[ ! -d "${resolution5_preprocess}" && \
      -d "${output_root}/xenium/preprocess/03i_resolution5_all_samples" ]]; then
  echo "STOP: resolution-5 preprocessing still uses its historical directory name." >&2
  echo "Run scripts/rename_resolution5_preprocess_output.sh first." >&2
  exit 1
fi

echo "Checking the authoritative resolution-5 replacement set..."
require_count \
  "Resolution-5 clustering RDS" \
  "${resolution5_preprocess}/rds" \
  '*_whole_tissue_Res5.0.rds'

for reference in aldinger sepp science; do
  reference_title=${reference^}
  require_count \
    "${reference_title} resolution-5 transfer RDS" \
    "${resolution5_transfer}/${reference}/rds" \
    "*_${reference_title}_annotated.rds"
  require_count \
    "${reference_title} transfer provenance" \
    "${resolution5_transfer}/${reference}/tables" \
    "*_${reference_title}_transfer_provenance.csv"
done

require_count \
  "Resolution-5 consensus tables" \
  "${resolution5_consensus_tables}" \
  '*_comparison_merged.csv'
require_count \
  "Resolution-5 consensus RDS" \
  "${resolution5_consensus}/rds" \
  '*_Consensus_annotated.rds'

for reference_path in \
  "${output_root}/references/aldinger/rds/Aldinger_newClusters_newUMAPv2_5k.rds" \
  "${output_root}/references/sepp/rds/Sepp_newClusters_newUMAPv2_5k.rds" \
  "${output_root}/references/science/rds/Science_newClusters_newUMAPv2_5k.rds"
do
  [[ -s "${reference_path}" ]] || {
    echo "STOP: current reference object is missing or empty: ${reference_path}" >&2
    exit 1
  }
done

# Exact, non-overlapping retired paths. Resolution-neutral prerequisites
# (01_cropped, 02_qc, 03_clustered) and current reference directories are not
# included because resolution-5 production depends on them.
retired_paths=(
  "outputs/export_smoke_test"
  "outputs/OLD"
  "outputs/validation_through_consensus_before_39824127"
  "outputs/validation/through_consensus"
  "outputs/validation/archive/validation_through_consensus_before_39824127"
  "outputs/references/science_old"
  "outputs/xenium/annotation/legacy_resolution1_5"
  "outputs/xenium/annotation/resolution4_all_samples"
  "outputs/xenium/annotation/resolution5_all_samples/01_label_transfer/science_old"
  "outputs/xenium/annotation/resolution5_all_samples/02_consensus_weighted_2of3_old"
  "outputs/xenium/annotation/resolution5_all_samples/03_consensus_labels_weighted_2of3_old"
  "outputs/xenium/preprocess/03b_resolution2_pilot"
  "outputs/xenium/preprocess/03c_resolution3_pilot"
  "outputs/xenium/preprocess/03d_resolution4_pilot"
  "outputs/xenium/preprocess/03e_resolution5_pilot"
  "outputs/xenium/preprocess/03g_resolution4_all_samples"
  "outputs/xenium/preprocess/03h_resolution5_selected_sample"
)

if [[ "${mode}" == "--delete-quarantine" ]]; then
  [[ -d "${quarantine}" ]] || {
    echo "Quarantine does not exist: ${quarantine}" >&2
    exit 1
  }
  [[ "${quarantine}" == "${project_root}/quarantine_non_resolution5_2026-09-29" ]] || {
    echo "Refusing unexpected deletion target: ${quarantine}" >&2
    exit 1
  }
  echo "Permanently deleting the reviewed quarantine: ${quarantine}"
  find "${quarantine}" -depth -delete
  echo "Deleted quarantine. This cannot be recovered from the project filesystem."
  exit 0
fi

echo
echo "Retired paths selected for ${mode}:"
found=0
for relative_path in "${retired_paths[@]}"; do
  source_path="${project_root}/${relative_path}"
  if [[ -e "${source_path}" ]]; then
    found=$((found + 1))
    du -sh "${source_path}"
  else
    echo "MISSING/SKIP  ${relative_path}"
  fi
done

[[ "${found}" -gt 0 ]] || {
  echo "No listed retired paths remain under the active project tree."
  exit 0
}

if [[ "${mode}" == "--dry-run" ]]; then
  echo
  echo "Dry run only. Re-run with --quarantine to move these exact paths."
  exit 0
fi

[[ ! -e "${quarantine}" ]] || {
  echo "STOP: quarantine already exists: ${quarantine}" >&2
  exit 1
}
mkdir -p "${quarantine}"

for relative_path in "${retired_paths[@]}"; do
  source_path="${project_root}/${relative_path}"
  [[ -e "${source_path}" ]] || continue
  destination_path="${quarantine}/${relative_path}"
  mkdir -p "$(dirname "${destination_path}")"
  mv -- "${source_path}" "${destination_path}"
  echo "Moved: ${relative_path}"
done

mkdir -p "$(dirname "${retained_manifest}")"
find "${quarantine}" -type f \
  -printf '%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%P\n' \
  | LC_ALL=C sort > "${retained_manifest}"
cp -- "${retained_manifest}" "${quarantine}/quarantine_manifest.tsv"

echo
du -sh "${quarantine}"
echo "Quarantine complete: ${quarantine}"
echo "Retained manifest: ${retained_manifest}"
echo "Review the manifest and run production dry-runs before permanent deletion."
