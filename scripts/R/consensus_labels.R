# Shared helpers for detailed consensus labels.

is_unknown_consensus_label <- function(labels) {
  labels <- trimws(as.character(labels))
  !is.na(labels) & (tolower(labels) == "unknown" | grepl("^Unknown-[0-9]{2,}$", labels))
}

sample_qualify_unknown_consensus_labels <- function(labels, sample_ids) {
  labels <- trimws(as.character(labels))
  sample_ids <- trimws(as.character(sample_ids))
  if (length(sample_ids) == 1L && length(labels) != 1L) {
    sample_ids <- rep(sample_ids, length(labels))
  }
  if (length(labels) != length(sample_ids)) {
    stop("Consensus labels and sample IDs must have the same length.")
  }
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Consensus labels must not be blank or missing.")
  }
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("Sample IDs must not be blank or missing.")
  }

  qualified <- labels
  unknown <- is_unknown_consensus_label(labels)
  qualified[unknown] <- paste0(labels[unknown], "_", sample_ids[unknown])
  qualified
}

pool_unknown_consensus_labels <- function(labels) {
  labels <- trimws(as.character(labels))
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Consensus labels must not be blank or missing.")
  }
  labels[is_unknown_consensus_label(labels)] <- "Unknown"
  labels
}

add_regional_unknown_labels <- function(object, sample_ids) {
  if (!"consensus_label" %in% colnames(object[[]])) {
    stop("Object lacks consensus_label metadata.")
  }
  original <- as.character(object$consensus_label)
  qualified <- sample_qualify_unknown_consensus_labels(original, sample_ids)
  pooled <- pool_unknown_consensus_labels(original)
  unknown <- is_unknown_consensus_label(original)

  object$regional_consensus_label <- qualified
  object$regional_consensus_label_pooled <- pooled
  object$regional_unknown_label <- ifelse(unknown, qualified, "Non-Unknown")
  object
}

validate_regional_unknown_labels <- function(object, sample_ids) {
  required <- c(
    "consensus_label", "regional_consensus_label",
    "regional_consensus_label_pooled", "regional_unknown_label"
  )
  missing <- setdiff(required, colnames(object[[]]))
  if (length(missing)) return(FALSE)

  original <- as.character(object$consensus_label)
  qualified <- sample_qualify_unknown_consensus_labels(original, sample_ids)
  pooled <- pool_unknown_consensus_labels(original)
  unknown_only <- ifelse(
    is_unknown_consensus_label(original), qualified, "Non-Unknown"
  )

  identical(unname(as.character(object$regional_consensus_label)), unname(qualified)) &&
    identical(unname(as.character(object$regional_consensus_label_pooled)), unname(pooled)) &&
    identical(unname(as.character(object$regional_unknown_label)), unname(unknown_only))
}

order_consensus_levels <- function(labels, known_order) {
  labels <- unique(trimws(as.character(labels)))
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Consensus labels must not be blank or missing.")
  }

  unknown_labels <- labels[is_unknown_consensus_label(labels)]
  unknown_number <- suppressWarnings(as.integer(sub("^Unknown-", "", unknown_labels)))
  unknown_labels <- unknown_labels[
    order(is.na(unknown_number), unknown_number, unknown_labels)
  ]
  known_labels <- labels[!is_unknown_consensus_label(labels)]

  c(
    intersect(known_order, known_labels),
    setdiff(known_labels, known_order),
    unknown_labels
  )
}
