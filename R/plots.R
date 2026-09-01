#' Plot a biosynthetic pathway
#'
#' @param x A `glyco_pathway` object.
#' @param show_rate_labels Show effective rate-class labels on reaction edges.
#' @param ... Additional graphical parameters passed to the graph plot.
#' @export
plot.glyco_pathway <- function(x, show_rate_labels = TRUE, ...) {
  validate_pathway(x)
  if (requireNamespace("igraph", quietly = TRUE)) {
    graph <- igraph::graph_from_data_frame(x$edges[, c("from", "to")],
                                           directed = TRUE,
                                           vertices = x$nodes[, c("id", "label")])
    edge_labels <- if (show_rate_labels) x$edges$rate_class else NA_character_
    layout <- tryCatch(igraph::layout_with_sugiyama(graph)$layout,
                       error = function(e) igraph::layout_in_circle(graph))
    graphics::plot(graph, layout = layout, vertex.shape = "circle",
                   vertex.color = "white", vertex.frame.color = "black",
                   vertex.label.color = "black", vertex.label.cex = 0.75,
                   vertex.size = 22, edge.color = "grey35",
                   edge.arrow.size = 0.35, edge.label = edge_labels,
                   edge.label.color = "black", edge.label.cex = 0.55,
                   margin = 0.08, ...)
  } else {
    n <- nrow(x$nodes)
    theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
    coordinates <- cbind(cos(theta), sin(theta))
    rownames(coordinates) <- x$nodes$id
    graphics::plot(coordinates, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, ...)
    from <- coordinates[x$edges$from, , drop = FALSE]
    to <- coordinates[x$edges$to, , drop = FALSE]
    graphics::arrows(from[, 1], from[, 2], to[, 1], to[, 2], length = 0.07, col = "grey35")
    graphics::points(coordinates, pch = 21, bg = "white", col = "black", cex = 2.4)
    graphics::text(coordinates, labels = x$nodes$label, cex = 0.65)
    if (show_rate_labels) {
      midpoint <- (from + to) / 2
      graphics::text(midpoint, labels = x$edges$rate_class, cex = 0.5, pos = 3)
    }
  }
  invisible(x)
}

#' Plot a panel-design audit
#'
#' @param x Object to plot.
#' @param ... Additional graphical parameters.
#' @export
plot.glyco_panel_audit <- function(x, ...) {
  values <- x$singular_values
  graphics::plot(seq_along(values), values, type = "b", pch = 16, col = "black",
                 xlab = "Singular-value index", ylab = "Local sensitivity",
                 main = "Resolvable rate directions", ...)
  graphics::abline(h = x$rank_tolerance, lty = 2, col = "grey55")
  invisible(x)
}

#' @rdname plot.glyco_panel_audit
#' @export
plot.glyco_robustness_audit <- function(x, ...) {
  values <- x$metrics$condition_number
  finite <- is.finite(values) & values > 0
  if (!any(finite)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No locally full-rank operating points")
    return(invisible(x))
  }
  graphics::boxplot(log10(values[finite]), horizontal = TRUE, col = "grey85",
                    border = "black", xlab = expression(log[10](condition~number)),
                    main = "Operating-point robustness", ...)
  invisible(x)
}

#' @rdname plot.glyco_panel_audit
#' @export
plot.glyco_panel_set <- function(x, ...) {
  values <- x$metrics$amplification
  finite <- is.finite(values) & values > 0
  if (!any(finite)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No locally full-rank panels")
    return(invisible(x))
  }
  graphics::boxplot(log10(values[finite]), horizontal = TRUE, col = "grey85",
                    border = "black", xlab = expression(log[10](noise~amplification~g)),
                    main = "Candidate-panel precision", ...)
  invisible(x)
}

#' @rdname plot.glyco_panel_audit
#' @export
plot.glyco_pathway_fit <- function(x, ...) {
  limits <- range(c(x$observed, x$predicted), finite = TRUE)
  graphics::plot(x$observed, x$predicted, pch = 16, col = "black", asp = 1,
                 xlim = limits, ylim = limits, xlab = "Observed composition",
                 ylab = "Predicted composition", main = "Forward-model fit", ...)
  graphics::abline(0, 1, col = "grey55", lty = 2)
  invisible(x)
}

#' Plot the local sensitivity matrix
#'
#' @param audit A `glyco_panel_audit` object.
#' @param ... Additional graphical parameters passed to `graphics::image`.
#' @return The input invisibly.
#' @export
plot_sensitivity <- function(audit, ...) {
  if (!inherits(audit, "glyco_panel_audit")) .gpd_stop("audit must be a glyco_panel_audit object.")
  z <- audit$jacobian
  palette <- grDevices::colorRampPalette(c("#2F4858", "white", "#A44A3F"))(101)
  bound <- max(abs(z))
  if (bound == 0) bound <- 1
  graphics::image(seq_len(nrow(z)), seq_len(ncol(z)), z,
                  zlim = c(-bound, bound), col = palette, axes = FALSE,
                  xlab = "Measured glycoforms", ylab = "Rate classes", ...)
  graphics::axis(1, at = seq_len(nrow(z)), labels = rownames(z), las = 2, cex.axis = 0.7)
  graphics::axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.75)
  graphics::box()
  invisible(audit)
}
