# Shared helpers for detailed consensus labels.

is_unknown_consensus_label <- function(labels) {
  labels <- trimws(as.character(labels))
  !is.na(labels) & (tolower(labels) == "unknown" | grepl("^Unknown-[0-9]{2,}$", labels))
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
