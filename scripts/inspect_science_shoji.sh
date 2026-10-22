#!/usr/bin/env bash
# Read-only targeted previews for the user-supplied Science /shoji HDF5 file.
# Expected runtime: 10 sec–2 min, low confidence (HPC filesystem dependent).
# --dry-run: path/tool checks only, under 5 sec. No R or new packages required.
set -euo pipefail

usage() {
  echo 'Usage: bash inspect_science_shoji.sh ORIGINAL_H5 [--dry-run]'
}
if [[ $# -lt 1 || $# -gt 2 ]]; then usage >&2; exit 2; fi
input_path=$1
if [[ ${2:-} != '' && ${2:-} != '--dry-run' ]]; then usage >&2; exit 2; fi
if [[ ! -r "$input_path" ]]; then echo "Input is not readable: $input_path" >&2; exit 1; fi
if ! command -v h5dump >/dev/null 2>&1; then echo 'h5dump is not available in this environment.' >&2; exit 1; fi
if [[ ${2:-} == '--dry-run' ]]; then
  echo "DRY-RUN PASS: readable input and h5dump available; no HDF5 contents loaded: $input_path"
  exit 0
fi

# Metadata attributes describe tensor dimensions and storage semantics.
# Dollar signs below are literal attribute names, protected by single quotes.
for attribute in \
  '/shoji/Dimension$cells' '/shoji/Dimension$genes' '/shoji/Dimension$clusters' \
  '/shoji/Tensor$Expression' '/shoji/Tensor$MeanExpression' \
  '/shoji/Tensor$Clusters' '/shoji/Tensor$ClusterID' \
  '/shoji/Tensor$Class' '/shoji/Tensor$CellClass' \
  '/shoji/Tensor$Gene' '/shoji/Tensor$Accession' \
  '/shoji/Tensor$AnnotationPosterior' '/shoji/Tensor$AnnotationName'; do
  echo "ATTRIBUTE: $attribute"
  h5dump -a "$attribute" "$input_path"
done

# These label/cluster tables are small: 51 annotations or 617 clusters,
# according to the supplied header. They are printed in full for mapping review.
for dataset in AnnotationName AnnotationDescription AnnotationDefinition ClusterID Class NCells; do
  echo "COMPLETE SMALL REFERENCE TABLE: /shoji/$dataset"
  h5dump -d "/shoji/$dataset" "$input_path"
done

# Larger cell/gene metadata arrays: inspect only the first eight entries.
for dataset in CellID Clusters CellClass Gene Accession Donor SampleID Region Subregion Subdivision Tissue Age TotalUMIs; do
  echo "FIRST EIGHT ENTRIES: /shoji/$dataset"
  h5dump -d "/shoji/$dataset" -s 0 -c 8 "$input_path"
done

# Counts and mean-expression previews are deliberately tiny (2 rows x 8 columns).
# Header order is cells x genes for Expression and clusters x genes for MeanExpression.
for dataset in Expression MeanExpression; do
  echo "TWO ROWS BY EIGHT GENES ONLY: /shoji/$dataset"
  h5dump -d "/shoji/$dataset" -s 0,0 -c 2,8 "$input_path"
done

echo 'Inspection complete. The original file is unchanged; no full expression matrix was read.'
echo 'ClusterID/Class, Clusters/CellClass, and AnnotationPosterior/AnnotationName mappings still require validation.'
echo 'Integer-looking Expression previews alone do not establish raw-count provenance.'
