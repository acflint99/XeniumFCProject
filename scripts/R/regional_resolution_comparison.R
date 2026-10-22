# Exact cell-aligned transitions between partitions; integer labels are arbitrary.
regional_resolution_contingency <- function(from, to, cell_ids) {
  if (!length(cell_ids) || anyNA(cell_ids) || any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids) || length(from) != length(cell_ids) ||
      length(to) != length(cell_ids) || anyNA(from) || anyNA(to) ||
      any(!nzchar(as.character(from))) || any(!nzchar(as.character(to)))) {
    stop("Contingency requires unique cell IDs and complete, aligned cluster labels.")
  }
  result <- as.data.frame(table(from_cluster = as.character(from),
                                to_cluster = as.character(to)), stringsAsFactors = FALSE)
  names(result)[3] <- "cells"
  result <- result[result$cells > 0L, , drop = FALSE]
  if (sum(result$cells) != length(cell_ids)) stop("Contingency lost cells.")
  result$fraction_of_from <- result$cells / ave(result$cells, result$from_cluster, FUN = sum)
  result$fraction_of_to <- result$cells / ave(result$cells, result$to_cluster, FUN = sum)
  result
}
