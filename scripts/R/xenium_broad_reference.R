# Shared helpers for panel-constrained, broad-level reference scoring.

broad_reference_names <- function() c("Aldinger", "Sepp", "Science")

broad_reference_categories <- function() c("VZ", "RL", "Other")

broad_category_to_route <- function(categories) {
  categories <- trimws(as.character(categories))
  routes <- c(VZ = "route_vz", RL = "route_rl", Other = "route_other")
  result <- unname(routes[categories])
  if (anyNA(result)) {
    stop(
      "Unsupported broad category or missing value: ",
      paste(sort(unique(categories[is.na(result)])), collapse = ", ")
    )
  }
  result
}

load_reference_broad_taxonomy <- function(config = load_pipeline_config()) {
  path <- resolve_config_path(config$manifests$reference_broad_taxonomy, config)
  if (!file.exists(path)) stop("Reference broad-taxonomy manifest not found: ", path)
  taxonomy <- read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  required <- c("fine_label", "broad_label", "review_notes")
  missing <- setdiff(required, names(taxonomy))
  if (length(missing)) {
    stop("Reference broad taxonomy lacks column(s): ", paste(missing, collapse = ", "))
  }
  taxonomy <- taxonomy[required]
  taxonomy$fine_label <- trimws(as.character(taxonomy$fine_label))
  taxonomy$broad_label <- trimws(as.character(taxonomy$broad_label))
  taxonomy$review_notes <- trimws(as.character(taxonomy$review_notes))
  if (!nrow(taxonomy) || anyNA(taxonomy$fine_label) ||
      any(!nzchar(taxonomy$fine_label)) || anyDuplicated(taxonomy$fine_label)) {
    stop("Reference broad taxonomy must contain unique, nonblank fine labels.")
  }
  if (anyNA(taxonomy$broad_label) ||
      any(!taxonomy$broad_label %in% broad_reference_categories())) {
    stop("Reference broad taxonomy may use only VZ, RL, and Other.")
  }
  attr(taxonomy, "path") <- normalizePath(path, mustWork = TRUE)
  taxonomy
}

validate_taxonomy_against_config <- function(taxonomy,
                                             config = load_pipeline_config()) {
  configured <- data.frame(
    fine_label = c(
      unlist(config$regional_subsets$vz_broad_labels, use.names = FALSE),
      unlist(config$regional_subsets$rl_broad_labels, use.names = FALSE),
      unlist(config$regional_subsets$nonregional_broad_labels, use.names = FALSE)
    ),
    broad_label = c(
      rep("VZ", length(config$regional_subsets$vz_broad_labels)),
      rep("RL", length(config$regional_subsets$rl_broad_labels)),
      rep("Other", length(config$regional_subsets$nonregional_broad_labels))
    ),
    stringsAsFactors = FALSE
  )
  configured <- configured[order(configured$fine_label), , drop = FALSE]
  observed <- taxonomy[c("fine_label", "broad_label")]
  observed <- observed[order(observed$fine_label), , drop = FALSE]
  if (!identical(configured, observed)) {
    stop(
      "config/reference_broad_taxonomy.csv and config.yml regional label ",
      "definitions disagree. Review both before continuing."
    )
  }
  invisible(TRUE)
}

prediction_score_matrix <- function(object) {
  metadata <- object[[]]
  score_columns <- grep(
    "^prediction\\.score\\.", names(metadata), value = TRUE
  )
  score_columns <- setdiff(
    score_columns, c("prediction.score.max", "prediction.score.id")
  )
  if (length(score_columns) < 2L) {
    stop("Annotated object contains fewer than two class prediction-score columns.")
  }
  matrix <- as.matrix(metadata[, score_columns, drop = FALSE])
  storage.mode(matrix) <- "double"
  if (any(!is.finite(matrix))) stop("Prediction-score matrix contains non-finite values.")
  labels <- sub("^prediction\\.score\\.", "", score_columns)
  if (anyDuplicated(labels) || any(!nzchar(labels))) {
    stop("Prediction-score columns do not encode unique, nonblank class labels.")
  }
  colnames(matrix) <- labels
  matrix
}

collapse_prediction_scores <- function(score_matrix, taxonomy) {
  missing_taxonomy <- setdiff(colnames(score_matrix), taxonomy$fine_label)
  unused_taxonomy <- setdiff(taxonomy$fine_label, colnames(score_matrix))
  if (length(missing_taxonomy)) {
    stop(
      "Reference score label(s) are absent from the broad taxonomy: ",
      paste(sort(missing_taxonomy), collapse = ", ")
    )
  }
  broad <- matrix(
    0, nrow = nrow(score_matrix), ncol = length(broad_reference_categories()),
    dimnames = list(rownames(score_matrix), broad_reference_categories())
  )
  label_index <- match(colnames(score_matrix), taxonomy$fine_label)
  for (category in broad_reference_categories()) {
    columns <- which(taxonomy$broad_label[label_index] == category)
    if (length(columns)) {
      broad[, category] <- rowSums(score_matrix[, columns, drop = FALSE])
    }
  }
  if (any(!is.finite(broad)) || any(broad < -1e-8)) {
    stop("Collapsed broad prediction scores are invalid.")
  }
  row_totals <- rowSums(broad)
  if (any(!is.finite(row_totals)) || any(row_totals <= 0)) {
    stop("At least one cell has no finite positive broad prediction-score mass.")
  }
  attr(broad, "unused_taxonomy_labels") <- unused_taxonomy
  broad
}

summarize_cluster_broad_scores <- function(object, taxonomy) {
  if (!"seurat_clusters" %in% names(object[[]])) {
    stop("Annotated object lacks seurat_clusters metadata.")
  }
  clusters <- trimws(as.character(object$seurat_clusters))
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    stop("Annotated object contains blank or missing cluster IDs.")
  }
  fine_scores <- prediction_score_matrix(object)
  broad_scores <- collapse_prediction_scores(fine_scores, taxonomy)
  cell_winner_index <- max.col(broad_scores, ties.method = "first")
  cell_winner <- colnames(broad_scores)[cell_winner_index]
  cluster_levels <- unique(clusters)
  numeric_clusters <- suppressWarnings(as.numeric(cluster_levels))
  cluster_levels <- cluster_levels[order(
    is.na(numeric_clusters), numeric_clusters, cluster_levels
  )]

  rows <- lapply(cluster_levels, function(cluster_id) {
    index <- which(clusters == cluster_id)
    means <- colMeans(broad_scores[index, , drop = FALSE])
    if (!is.finite(sum(means)) || sum(means) <= 0) {
      stop("Cluster ", cluster_id, " has no positive broad prediction-score mass.")
    }
    ordering <- order(means, decreasing = TRUE)
    winner <- names(means)[ordering[[1]]]
    runner_up <- names(means)[ordering[[2]]]
    probabilities <- means / sum(means)
    positive <- probabilities > 0
    entropy <- -sum(probabilities[positive] * log(probabilities[positive])) /
      log(length(probabilities))
    data.frame(
      seurat_clusters = cluster_id,
      n_cells = length(index),
      broad_winner = winner,
      broad_winner_route = broad_category_to_route(winner),
      broad_winner_mean_score = unname(means[[winner]]),
      broad_runner_up = runner_up,
      broad_runner_up_mean_score = unname(means[[runner_up]]),
      broad_top2_mean_score_margin = unname(means[[winner]] - means[[runner_up]]),
      broad_winner_cell_fraction = mean(cell_winner[index] == winner),
      broad_score_entropy = entropy,
      broad_mean_score_vz = unname(means[["VZ"]]),
      broad_mean_score_rl = unname(means[["RL"]]),
      broad_mean_score_other = unname(means[["Other"]]),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (anyDuplicated(result$seurat_clusters)) {
    stop("Broad-score summary produced duplicate cluster IDs.")
  }
  result
}

load_reviewed_broad_thresholds <- function(config = load_pipeline_config(),
                                           require_approved = TRUE) {
  path <- resolve_config_path(
    config$manifests$reviewed_broad_reference_thresholds, config
  )
  if (!file.exists(path)) stop("Reviewed broad-threshold manifest not found: ", path)
  table <- read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA")
  )
  required <- c(
    "reference", "min_broad_winner_score", "min_broad_top2_margin",
    "min_broad_winner_cell_fraction", "approved", "review_notes"
  )
  missing <- setdiff(required, names(table))
  if (length(missing)) stop("Broad-threshold manifest lacks: ", paste(missing, collapse = ", "))
  table <- table[required]
  table$reference <- trimws(as.character(table$reference))
  if (anyDuplicated(table$reference) ||
      !setequal(table$reference, broad_reference_names())) {
    stop("Broad-threshold manifest must contain exactly one row per reference.")
  }
  approved <- suppressWarnings(as.logical(table$approved))
  if (anyNA(approved)) stop("Broad-threshold approved values must be TRUE or FALSE.")
  table$approved <- approved
  numeric_columns <- c(
    "min_broad_winner_score", "min_broad_top2_margin",
    "min_broad_winner_cell_fraction"
  )
  for (column in numeric_columns) {
    table[[column]] <- suppressWarnings(as.numeric(table[[column]]))
  }
  if (require_approved) {
    if (any(!table$approved)) {
      stop("All broad-reference thresholds must be explicitly approved before application.")
    }
    if (any(!is.finite(as.matrix(table[numeric_columns]))) ||
        any(as.matrix(table[numeric_columns]) < 0) ||
        any(as.matrix(table[numeric_columns]) > 1)) {
      stop("Approved broad-reference thresholds must be finite values in [0, 1].")
    }
    if (any(!nzchar(trimws(as.character(table$review_notes))))) {
      stop("Every approved broad-reference threshold requires review notes.")
    }
  }
  attr(table, "path") <- normalizePath(path, mustWork = TRUE)
  table
}
