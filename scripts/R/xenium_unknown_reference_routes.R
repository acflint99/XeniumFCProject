# Helpers for deterministic broad routing of sample-qualified Unknown
# populations from the three existing reference-transfer votes.

unknown_reference_names <- function() {
  c("aldinger", "sepp", "science")
}

unknown_reference_vote_columns <- function() {
  paste0(unknown_reference_names(), "_reference_vote")
}

unknown_reference_route_values <- function() {
  c("route_vz", "route_rl", "route_other", "route_unresolved")
}

collapse_reference_vote <- function(labels, config = load_pipeline_config()) {
  labels <- trimws(as.character(labels))
  routes <- rep(NA_character_, length(labels))
  abstain <- is.na(labels) | !nzchar(labels) | tolower(labels) == "unknown"
  routes[abstain] <- "abstain"

  vz_labels <- unlist(config$regional_subsets$vz_broad_labels, use.names = FALSE)
  rl_labels <- unlist(config$regional_subsets$rl_broad_labels, use.names = FALSE)
  nonregional_labels <- unlist(
    config$regional_subsets$nonregional_broad_labels, use.names = FALSE
  )
  routes[!abstain & labels %in% vz_labels] <- "route_vz"
  routes[!abstain & labels %in% rl_labels] <- "route_rl"
  routes[!abstain & labels %in% nonregional_labels] <- "route_other"

  unrecognized <- unique(labels[!abstain & is.na(routes)])
  if (length(unrecognized)) {
    stop(
      "Reference vote label(s) are absent from the configured VZ, RL, and ",
      "nonregional taxonomies: ", paste(sort(unrecognized), collapse = ", ")
    )
  }
  routes
}

derive_reference_routes <- function(table, config = load_pipeline_config()) {
  vote_columns <- unknown_reference_vote_columns()
  required <- c("seurat_clusters", "consensus_label", vote_columns)
  missing <- setdiff(required, names(table))
  if (length(missing)) {
    stop("Consensus table lacks route column(s): ", paste(missing, collapse = ", "))
  }
  if (!nrow(table)) stop("Consensus table has no cluster rows.")
  if (anyDuplicated(as.character(table$seurat_clusters))) {
    stop("Consensus table contains duplicate seurat_clusters values.")
  }

  collapsed_columns <- paste0(unknown_reference_names(), "_broad_route")
  for (i in seq_along(vote_columns)) {
    table[[collapsed_columns[[i]]]] <- collapse_reference_vote(
      table[[vote_columns[[i]]]], config
    )
  }

  min_support <- as.integer(
    config$regional_subsets$unknown_route_min_reference_support
  )
  route_result <- lapply(seq_len(nrow(table)), function(i) {
    votes <- as.character(unlist(table[i, collapsed_columns, drop = FALSE]))
    valid <- votes[votes != "abstain"]
    if (!length(valid)) {
      return(list(
        route = "route_unresolved", support = 0L, nonabstain = 0L,
        status = "all_references_abstained"
      ))
    }
    counts <- base::table(valid)
    support <- max(counts)
    winners <- names(counts)[counts == support]
    if (support < min_support || length(winners) != 1L) {
      return(list(
        route = "route_unresolved", support = as.integer(support),
        nonabstain = length(valid), status = "insufficient_broad_support"
      ))
    }
    list(
      route = winners[[1]], support = as.integer(support),
      nonabstain = length(valid),
      status = if (support == 3L) "unanimous_3of3" else "majority_2of3"
    )
  })

  table$reference_route <- vapply(route_result, `[[`, character(1), "route")
  table$reference_route_support_n <- vapply(
    route_result, `[[`, integer(1), "support"
  )
  table$reference_route_nonabstain_n <- vapply(
    route_result, `[[`, integer(1), "nonabstain"
  )
  table$reference_route_status <- vapply(
    route_result, `[[`, character(1), "status"
  )
  table
}

unknown_reference_route_paths <- function(config = load_pipeline_config()) {
  production_paths <- get_production_output_paths(config)
  stage_dir <- file.path(
    production_paths$annotation_root, "04_unknown_reference_routing"
  )
  table_dir <- file.path(stage_dir, "tables")
  list(
    stage_dir = stage_dir,
    table_dir = table_dir,
    route_manifest = file.path(
      table_dir, "Xenium_Unknown_Reference_Routes.csv"
    ),
    input_manifest = file.path(
      table_dir, "Xenium_Unknown_Reviewed_Input_manifest.csv"
    ),
    validation = file.path(
      table_dir, "Xenium_Reference_Route_Validation.csv"
    ),
    summary = file.path(
      table_dir, "Xenium_Unknown_Reference_Route_Summary.csv"
    )
  )
}

validate_unknown_reference_route_manifest <- function(routes, sample_ids) {
  required <- c(
    "sample_id", "seurat_clusters", "consensus_label",
    unknown_reference_vote_columns(),
    paste0(unknown_reference_names(), "_broad_route"),
    "reference_route", "reference_route_support_n",
    "reference_route_nonabstain_n", "reference_route_status", "qc_decision",
    "consensus_table_path", "consensus_table_size", "consensus_table_mtime",
    "broad_consensus_table_path", "broad_consensus_table_size",
    "broad_consensus_table_mtime"
  )
  missing <- setdiff(required, names(routes))
  if (length(missing)) {
    stop("Unknown route manifest lacks column(s): ", paste(missing, collapse = ", "))
  }
  if (!nrow(routes)) stop("Unknown route manifest has no rows.")
  if (any(!routes$sample_id %in% sample_ids)) {
    stop("Unknown route manifest contains sample IDs absent from config/samples.csv.")
  }
  keys <- paste(
    trimws(as.character(routes$sample_id)),
    trimws(as.character(routes$consensus_label)), sep = "\u001f"
  )
  if (anyDuplicated(keys)) stop("Unknown route manifest contains duplicate keys.")
  if (any(!is_unknown_consensus_label(routes$consensus_label))) {
    stop("Unknown route manifest contains a non-Unknown consensus label.")
  }
  if (any(!routes$reference_route %in% unknown_reference_route_values())) {
    stop("Unknown route manifest contains an unsupported reference_route value.")
  }
  support <- suppressWarnings(as.integer(routes$reference_route_support_n))
  nonabstain <- suppressWarnings(as.integer(routes$reference_route_nonabstain_n))
  if (anyNA(support) || anyNA(nonabstain) || any(support < 0L | support > 3L) ||
      any(nonabstain < 0L | nonabstain > 3L) || any(support > nonabstain)) {
    stop("Unknown route manifest contains invalid reference-support counts.")
  }
  resolved <- routes$reference_route != "route_unresolved"
  if (any(resolved & support < 2L)) {
    stop("A resolved Unknown route has fewer than two supporting references.")
  }
  table_size <- suppressWarnings(as.numeric(routes$consensus_table_size))
  table_mtime <- suppressWarnings(as.numeric(routes$consensus_table_mtime))
  if (any(!is.finite(table_size)) || any(table_size <= 0) ||
      any(!is.finite(table_mtime))) {
    stop("Unknown route manifest contains invalid consensus-table provenance.")
  }
  broad_table_size <- suppressWarnings(as.numeric(routes$broad_consensus_table_size))
  broad_table_mtime <- suppressWarnings(as.numeric(routes$broad_consensus_table_mtime))
  if (any(!is.finite(broad_table_size)) || any(broad_table_size <= 0) ||
      any(!is.finite(broad_table_mtime))) {
    stop("Unknown route manifest contains invalid broad-consensus provenance.")
  }
  decisions <- tolower(trimws(as.character(routes$qc_decision)))
  if (any(!decisions %in% c("keep", "remove"))) {
    stop("Unknown route manifest contains a QC decision other than keep/remove.")
  }
  invisible(routes)
}

load_unknown_reference_routes <- function(config = load_pipeline_config()) {
  path <- unknown_reference_route_paths(config)$route_manifest
  if (!file.exists(path)) stop("Unknown reference-route manifest not found: ", path)
  routes <- read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  sample_ids <- load_sample_manifest(config)$sample_id
  validate_unknown_reference_route_manifest(routes, sample_ids)
  route_samples <- unique(as.character(routes$sample_id))
  provenance <- routes[match(route_samples, routes$sample_id), c(
    "sample_id", "consensus_table_path", "consensus_table_size",
    "consensus_table_mtime", "broad_consensus_table_path",
    "broad_consensus_table_size", "broad_consensus_table_mtime"
  ), drop = FALSE]
  for (sample_id in route_samples) {
    sample_rows <- routes[routes$sample_id == sample_id, , drop = FALSE]
    signatures <- unique(paste(
      sample_rows$consensus_table_path,
      sample_rows$consensus_table_size,
      sample_rows$consensus_table_mtime, sep = "\u001f"
    ))
    if (length(signatures) != 1L) {
      stop("Reference-route provenance is inconsistent within sample ", sample_id, ".")
    }
  }
  table_paths <- as.character(provenance$consensus_table_path)
  missing_tables <- table_paths[!file.exists(table_paths)]
  if (length(missing_tables)) {
    stop("Reference-route consensus table(s) are missing:\n- ",
         paste(missing_tables, collapse = "\n- "))
  }
  table_info <- file.info(table_paths)
  current <- vapply(seq_along(table_paths), function(i) {
    isTRUE(all.equal(
      as.numeric(provenance$consensus_table_size[[i]]),
      as.numeric(table_info$size[[i]])
    )) && isTRUE(all.equal(
      as.numeric(provenance$consensus_table_mtime[[i]]),
      as.numeric(table_info$mtime[[i]])
    ))
  }, logical(1))
  if (any(!current)) {
    stop(
      "Reference-route consensus-table provenance changed for sample(s): ",
      paste(provenance$sample_id[!current], collapse = ", ")
    )
  }
  broad_table_paths <- as.character(provenance$broad_consensus_table_path)
  missing_broad_tables <- broad_table_paths[!file.exists(broad_table_paths)]
  if (length(missing_broad_tables)) {
    stop("Reference-route broad-consensus table(s) are missing:\n- ",
         paste(missing_broad_tables, collapse = "\n- "))
  }
  broad_info <- file.info(broad_table_paths)
  broad_current <- vapply(seq_along(broad_table_paths), function(i) {
    isTRUE(all.equal(
      as.numeric(provenance$broad_consensus_table_size[[i]]),
      as.numeric(broad_info$size[[i]])
    )) && isTRUE(all.equal(
      as.numeric(provenance$broad_consensus_table_mtime[[i]]),
      as.numeric(broad_info$mtime[[i]])
    ))
  }, logical(1))
  if (any(!broad_current)) {
    stop(
      "Reference-route broad-consensus provenance changed for sample(s): ",
      paste(provenance$sample_id[!broad_current], collapse = ", ")
    )
  }
  routes
}

load_unknown_reviewed_input_manifest <- function(config = load_pipeline_config()) {
  path <- unknown_reference_route_paths(config)$input_manifest
  if (!file.exists(path)) stop("Reviewed Unknown input manifest not found: ", path)
  manifest <- read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  sample_ids <- load_sample_manifest(config)$sample_id
  inspection <- inspect_unknown_qc_routing_manifest(manifest, sample_ids)
  if (!inspection$valid) {
    stop("Reviewed Unknown input manifest is invalid:\n- ",
         paste(inspection$errors, collapse = "\n- "))
  }
  paths <- as.character(manifest$mapping_input_path)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Reviewed Unknown input file(s) are missing:\n- ", paste(missing, collapse = "\n- "))
  }
  info <- file.info(paths)
  current <- vapply(seq_along(paths), function(i) {
    isTRUE(all.equal(
      as.numeric(manifest$mapping_input_size[[i]]), as.numeric(info$size[[i]])
    )) && isTRUE(all.equal(
      as.numeric(manifest$mapping_input_mtime[[i]]), as.numeric(info$mtime[[i]])
    ))
  }, logical(1))
  if (any(!current)) {
    stop(
      "Reviewed Unknown input provenance changed for sample(s): ",
      paste(sample_ids[!current], collapse = ", ")
    )
  }
  review_paths <- unique(as.character(manifest$reviewed_qc_path))
  if (length(review_paths) != 1L || !file.exists(review_paths[[1]])) {
    stop("The reviewed Unknown-QC manifest recorded by Stage 03 is unavailable.")
  }
  review_info <- file.info(review_paths[[1]])
  review_current <- vapply(seq_len(nrow(manifest)), function(i) {
    isTRUE(all.equal(
      as.numeric(manifest$reviewed_qc_size[[i]]), as.numeric(review_info$size)
    )) && isTRUE(all.equal(
      as.numeric(manifest$reviewed_qc_mtime[[i]]), as.numeric(review_info$mtime)
    ))
  }, logical(1))
  if (any(!review_current)) {
    stop("The reviewed Unknown-QC manifest changed after Stage 03 application.")
  }
  manifest
}
