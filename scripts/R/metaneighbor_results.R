# MetaNeighborUS result validation/formatting. No inferred renaming of source labels.
mn_pair_scores <- function(aurocs, query_labels, reference_labels, reference) {
  q <- paste("xenium", query_labels, sep = "|")
  r <- paste(reference, reference_labels, sep = "|")
  expected <- c(q, r)
  if (!is.matrix(aurocs) || anyDuplicated(rownames(aurocs)) || anyDuplicated(colnames(aurocs)) ||
      anyNA(match(expected, rownames(aurocs))) || anyNA(match(expected, colnames(aurocs)))) stop("MetaNeighbor result omitted/duplicated original populations.")
  forward <- aurocs[q, r, drop = FALSE]
  reverse <- t(aurocs[r, q, drop = FALSE])
  if (any(!is.finite(forward)) || any(!is.finite(reverse)) || any(forward < 0 | forward > 1) || any(reverse < 0 | reverse > 1)) stop("Non-finite/out-of-range cross-dataset AUROC.")
  dimnames(forward) <- dimnames(reverse) <- list(query_labels, reference_labels)
  mean_score <- (forward + reverse) / 2
  full <- expand.grid(query_label = query_labels, reference_label = reference_labels, stringsAsFactors = FALSE)
  full$query_recovery_AUROC <- as.vector(forward)
  full$reference_recovery_AUROC <- as.vector(reverse)
  full$mean_bidirectional_AUROC <- as.vector(mean_score)
  # Reciprocal maxima retain ties; ties are reported, never a unique forced match.
  q_best <- forward >= apply(forward, 1L, max) - 1e-12
  r_best <- sweep(reverse, 2L, apply(reverse, 2L, max), "-") >= -1e-12
  full$query_direction_best <- as.vector(q_best)
  full$reference_direction_best <- as.vector(r_best)
  full$reciprocal_best <- full$query_direction_best & full$reference_direction_best
  best <- lapply(seq_along(query_labels), function(i) {
    ix <- order(-forward[i, ], reference_labels)
    top <- ix[1L]; second <- ix[2L]
    data.frame(query_label = query_labels[i], best_reference_label = reference_labels[top],
               runner_up = reference_labels[second], query_recovery_AUROC = forward[i, top],
               reference_recovery_AUROC = reverse[i, top], mean_bidirectional_AUROC = mean_score[i, top],
               top2_query_AUROC_margin = forward[i, top] - forward[i, second],
               tied_best_references = sum(q_best[i, ]), reciprocal_best = r_best[i, top],
               interpretation = "Transcriptional correspondence; AUROC is not annotation accuracy or a cell-type probability")
  })
  list(full = full, best = do.call(rbind, best), forward = forward, reverse = reverse, mean = mean_score)
}
mn_heatmaps <- function(scores, file, title) {
  ql <- rownames(scores$mean); rl <- colnames(scores$mean)
  bottom <- max(1.5, .065 * max(nchar(rl)) + .45)
  left <- max(1.6, .07 * max(nchar(ql)) + .3)
  grDevices::pdf(file, width = max(10, left + length(rl) * .23 + 1.8),
                 height = max(7, bottom + length(ql) * .24 + 1.1), useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  colors <- grDevices::colorRampPalette(c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F"))(255)
  for (metric in c("forward", "reverse", "mean")) {
    m <- scores[[metric]][rev(ql), , drop = FALSE]
    graphics::layout(matrix(1:2, nrow = 1), widths = c(1, .13))
    graphics::par(mai = c(bottom, left, 1.25, .15))
    graphics::image(seq_along(rl), seq_along(ql), t(m), col = colors, zlim = c(0,1), axes = FALSE, xlab = "", ylab = "")
    best <- m >= apply(m, 1L, max) - 1e-12
    mark <- which(best, arr.ind = TRUE)
    graphics::text(mark[, "col"], mark[, "row"], "*", cex = .9, font = 2, col = "#FFFFFF")
    graphics::axis(1, seq_along(rl), rl, las = 2, cex.axis = .8, tick = FALSE)
    graphics::axis(2, seq_along(ql), rev(ql), las = 2, cex.axis = .8, tick = FALSE)
    page_title <- c(forward = "Recovery of Xenium populations from reference", reverse = "Recovery of reference populations from Xenium", mean = "Mean of both directional AUROCs")[[metric]]
    graphics::title(main = title, cex.main = 1)
    graphics::mtext("* = row maximum (ties retained)", side = 3, line = 1.35, adj = 0, cex = .62)
    graphics::mtext(page_title, side = 3, line = .5, cex = .75)
    graphics::box(col = "#AAAAAA")
    graphics::par(mai = c(bottom, .12, .95, .62))
    z <- seq(0, 1, length.out = 256)
    graphics::image(c(0,1), z, matrix(rep(z, each = 2), nrow = 2), col = colors, zlim = c(0,1), axes = FALSE, xlab = "", ylab = "")
    graphics::axis(4, at = seq(0,1,.25), las = 2, cex.axis = .75)
    graphics::mtext("AUROC", side = 3, line = .3, cex = .75)
    graphics::box(col = "#AAAAAA")
  }
}
