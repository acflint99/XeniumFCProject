# Give each sample's FOVs stable names before any batching. Seurat's automatic
# duplicate-image suffixes can collide with existing suffixes at the next merge.
namespace_vz_rl_images <- function(object, sample_id) {
  if (!inherits(object, "Seurat") || length(sample_id) != 1L ||
      is.na(sample_id) || !nzchar(sample_id)) {
    stop("Expected a Seurat object and one nonblank sample ID.")
  }
  original_names <- SeuratObject::Images(object)
  if (!length(original_names) || anyNA(original_names) ||
      any(!nzchar(original_names)) || anyDuplicated(original_names)) {
    stop(sample_id, " requires nonempty, unique image/FOV names.")
  }
  new_names <- paste0(sample_id, "__", original_names)
  if (any(new_names %in% setdiff(names(object), original_names))) {
    stop(sample_id, " sample-qualified FOV names conflict with other object names.")
  }
  # Rename only the named image-list entries. The FOVs themselves, their keys,
  # coordinates, and cell memberships remain unchanged.
  names(object@images) <- new_names
  object
}

# Merge a small ordered collection with the same Seurat options used by the
# former pairwise spatial merge. Keeping this in one function makes it possible
# to compare batched and pairwise results on a small fixture before a full run.
merge_vz_rl_objects <- function(objects, project = "Xenium_VZRL_Spatial") {
  if (!is.list(objects) || !length(objects) ||
      any(!vapply(objects, inherits, logical(1), what = "Seurat"))) {
    stop("Expected a nonempty list of Seurat objects.")
  }
  if (length(objects) == 1L) return(objects[[1L]])
  expected_images <- unlist(lapply(objects, SeuratObject::Images), use.names = FALSE)
  if (anyDuplicated(expected_images)) {
    stop("Image/FOV names overlap across merge inputs; qualify them by sample before batching.")
  }
  merged <- merge(
    x = objects[[1L]], y = objects[-1L], project = project,
    collapse = FALSE, merge.data = TRUE, merge.dr = FALSE
  )
  if (!identical(SeuratObject::Images(merged), expected_images)) {
    stop("Spatial merge lost, renamed, or reordered image/FOV entries.")
  }
  merged
}
