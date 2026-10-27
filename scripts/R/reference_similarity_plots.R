# Read the shared literal order vectors without sourcing unrelated plot helpers.
# This keeps color_palette.R authoritative and adds no plotting dependencies.
reference_similarity_order <- function(labels, level, palette_file = here::here("scripts", "color_palette.R"), legacy_query_labels = FALSE) {
  labels <- as.character(labels)
  if (!length(labels) || anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)) stop("Invalid/duplicate heatmap query labels.")
  vector_name <- switch(level, comb_subcluster = "master_subcluster_order",
                        corrected_broad_label = "celltype_order", stop("Unknown query level: ", level))
  expressions <- as.list(parse(palette_file))
  assignments <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
                          identical(x[[2L]], as.name(vector_name)), expressions)
  if (length(assignments) != 1L) stop("Expected one shared order definition: ", vector_name)
  rhs <- assignments[[1L]][[3L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("c"))) stop("Shared order must be a literal character vector: ", vector_name)
  items <- as.list(rhs)[-1L]
  if (!length(items) || any(!vapply(items, function(x) is.character(x) && length(x) == 1L && !is.na(x), FALSE))) stop("Shared order contains nonliteral values: ", vector_name)
  preferred <- unlist(items, use.names = FALSE)
  if (anyDuplicated(preferred) || any(!nzchar(trimws(preferred)))) stop("Invalid shared order: ", vector_name)
  extras <- setdiff(labels, preferred)
  # Numbered Unknown groups are established dynamic identities; keep them last.
  unknown <- extras[grepl("^Unknown-[0-9]{2,}$", extras)]
  legacy <- if (isTRUE(legacy_query_labels)) intersect(extras, c("Endothelial", "Pericytes")) else character()
  unexpected <- setdiff(extras, c(unknown, legacy))
  if (length(unexpected)) stop("Query labels absent from ", vector_name, ": ", paste(unexpected, collapse = ", "), ". Update the shared order explicitly; no label is omitted.")
  ordered <- preferred[preferred %in% labels]
  if (length(legacy)) {
    # Legacy-only display placement: immediately before the harmonized Vascular identity.
    legacy <- c("Endothelial", "Pericytes")[c("Endothelial", "Pericytes") %in% legacy]
    vascular <- match("Vascular", ordered)
    if (is.na(vascular)) ordered <- c(ordered, legacy) else ordered <- append(ordered, legacy, after = vascular - 1L)
  }
  c(ordered, sort(unknown))
}

# Plotting only: no change to correlation values. Typical redraw: <1-2 minutes.
write_reference_heatmaps <- function(raw, specific, file, title, reference_display = colnames(raw), query_level = NULL, legacy_query_labels = FALSE) {
  if (!identical(dim(raw), dim(specific)) || any(!is.finite(raw)) || any(!is.finite(specific)) ||
      is.null(rownames(raw)) || length(reference_display) != ncol(raw)) stop("Invalid heatmap matrix/labels.")
  top_order <- if (is.null(query_level)) rownames(raw) else reference_similarity_order(rownames(raw), query_level, legacy_query_labels = legacy_query_labels)
  if (!identical(dimnames(raw), dimnames(specific))) stop("Heatmap matrix labels differ.")
  # Base image() draws the first row at the bottom; reverse to put shared first at top.
  raw <- raw[rev(top_order), , drop = FALSE]
  specific <- specific[rev(top_order), , drop = FALSE]
  bottom <- max(1.5, .065 * max(nchar(reference_display)) + .45)
  left <- max(1.6, .070 * max(nchar(rownames(raw))) + .3)
  grDevices::pdf(file, width = max(10, left + ncol(raw) * .23 + 1.9),
                 height = max(7, bottom + nrow(raw) * .24 + 1.1), useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  palette <- grDevices::colorRampPalette(c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F"))(255)
  for (method in c("raw", "specificity")) {
    m <- if (method == "raw") raw else specific
    limits <- if (method == "raw") range(m) else c(-1, 1) * max(abs(m))
    if (diff(limits) <= .Machine$double.eps) limits <- limits + c(-.01, .01)
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, .13))
    graphics::par(mai = c(bottom, left, 1.25, .15))
    graphics::image(seq_len(ncol(m)), seq_len(nrow(m)), t(m), zlim = limits,
      col = palette, axes = FALSE, xlab = "", ylab = "")
    best <- m >= apply(m, 1L, max) - 1e-12
    mark <- which(best, arr.ind = TRUE)
    graphics::text(mark[, "col"], mark[, "row"], "*", cex = .9, font = 2, col = "#FFFFFF")
    graphics::axis(1, seq_len(ncol(m)), reference_display, las = 2, cex.axis = .8, tick = FALSE)
    graphics::axis(2, seq_len(nrow(m)), rownames(m), las = 2, cex.axis = .8, tick = FALSE)
    graphics::box(col = "#AAAAAA")
    graphics::title(main = paste(title, if (method == "raw") "Raw expression similarity" else "Population specificity", sep = " | "), cex.main = 1)
    graphics::mtext("* = row maximum (ties retained)", side = 3, line = 1.35, adj = 0, cex = .62)
    graphics::mtext(if (method == "raw") "Blue = lower correlation; red = higher. Colors span the observed range; blue may be positive."
                    else "Gene-centered expression: blue = negative; white = zero; red = positive correlation.", side = 3, line = .4, cex = .68)
    graphics::par(mai = c(bottom, .12, .95, .62))
    z <- seq(limits[1], limits[2], length.out = 256)
    graphics::image(c(0, 1), z, matrix(rep(z, each = 2), nrow = 2), col = palette,
                    zlim = limits, axes = FALSE, xlab = "", ylab = "")
    ticks <- sort(unique(c(limits, pretty(limits, n = 4))))
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    graphics::axis(4, at = ticks, labels = format(round(ticks, 2), nsmall = 2), las = 2, cex.axis = .75)
    graphics::mtext("Spearman rho", side = 3, line = .3, cex = .75)
    graphics::box(col = "#AAAAAA")
  }
  invisible(top_order)
}

# Cache the unchanged reference's ranks once; Pearson of ranks is Spearman.
make_reference_spearman <- function(reference) {
  normalize_ranks <- function(x) {
    ranked <- apply(x, 2L, rank, ties.method = "average")
    ranked <- matrix(ranked, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    ranked <- sweep(ranked, 2L, colMeans(ranked), "-")
    norms <- sqrt(colSums(ranked * ranked))
    if (any(!is.finite(norms)) || any(norms == 0)) stop("Undefined correlation: constant profile.")
    sweep(ranked, 2L, norms, "/")
  }
  rr <- normalize_ranks(reference)
  function(query) {
    if (nrow(query) != nrow(reference) || any(!is.finite(query))) stop("Invalid correlation input.")
    ans <- crossprod(normalize_ranks(query), rr)
    # Bound floating-point roundoff at the correlation endpoints.
    ans[] <- pmax(-1, pmin(1, ans))
    ans
  }
}
