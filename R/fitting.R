.prepare_observed <- function(observed, pathway, panel = NULL) {
  if (is.data.frame(observed)) {
    .require_columns(observed, c("glycoform", "proportion"), "observed")
    values <- observed$proportion
    names(values) <- as.character(observed$glycoform)
    observed <- values
  }
  original_names <- names(observed)
  observed <- as.numeric(observed)
  names(observed) <- original_names
  if (is.null(names(observed))) .gpd_stop("observed must be named or have glycoform and proportion columns.")
  panel <- as.character(panel %||% names(observed))
  missing <- setdiff(panel, names(observed))
  if (length(missing)) .gpd_stop("observed is missing panel class(es): %s.", paste(missing, collapse = ", "))
  unknown <- setdiff(panel, pathway$observed_classes)
  if (length(unknown)) .gpd_stop("Panel contains unknown class(es): %s.", paste(unknown, collapse = ", "))
  setNames(.normalise(observed[panel], "observed composition"), panel)
}

#' Fit effective pathway rates to an observed composition
#'
#' This function estimates effective log-rate parameters for a specified graph.
#' The estimates remain conditional on that graph, observation mapping, and
#' operating-point assumptions.
#'
#' @param pathway A `glyco_pathway` object.
#' @param observed Named numeric composition or a data frame with columns
#'   `glycoform` and `proportion`.
#' @param panel Classes used in fitting. Defaults to names in `observed`.
#' @param start Initial named log-rate vector.
#' @param lower,upper Optimisation bounds on the log-rate scale.
#' @param n_starts Number of optimisation starts.
#' @param seed Random seed for additional starts.
#' @param loss Squared error on proportions or log proportions.
#' @param pseudocount Pseudocount for log-composition loss.
#' @param control List passed to `stats::optim`.
#' @inheritParams forward_pathway
#'
#' @return A `glyco_pathway_fit` object.
#' @export
fit_pathway <- function(pathway, observed, panel = NULL, start = NULL,
                        lower = -5, upper = 5, n_starts = 8L, seed = 1L,
                        loss = c("composition", "log_composition"),
                        pseudocount = 1e-8, secretion = NULL, entry = NULL,
                        control = list(maxit = 2000)) {
  validate_pathway(pathway)
  loss <- match.arg(loss)
  observed <- .prepare_observed(observed, pathway, panel)
  k <- length(pathway$rate_classes)
  start <- .named_numeric(start, pathway$rate_classes, "start",
                          default = setNames(rep(0, k), pathway$rate_classes))
  lower <- if (length(lower) == 1L) rep(lower, k) else lower
  upper <- if (length(upper) == 1L) rep(upper, k) else upper
  lower <- .named_numeric(lower, pathway$rate_classes, "lower")
  upper <- .named_numeric(upper, pathway$rate_classes, "upper")
  if (any(lower >= upper)) .gpd_stop("Every lower bound must be below its upper bound.")
  objective <- function(par) {
    prediction <- panel_composition(predict_composition(pathway, par, secretion, entry), names(observed))
    if (loss == "composition") {
      sum((prediction - observed)^2)
    } else {
      sum((log(prediction + pseudocount) - log(observed + pseudocount))^2)
    }
  }
  set.seed(seed)
  n_starts <- max(1L, as.integer(n_starts))
  starts <- matrix(NA_real_, nrow = n_starts, ncol = k,
                   dimnames = list(NULL, pathway$rate_classes))
  starts[1L, ] <- pmin(pmax(start, lower), upper)
  if (n_starts > 1L) {
    starts[-1L, ] <- matrix(stats::runif((n_starts - 1L) * k,
                                         rep(lower, each = n_starts - 1L),
                                         rep(upper, each = n_starts - 1L)),
                            nrow = n_starts - 1L, byrow = FALSE)
  }
  fits <- lapply(seq_len(n_starts), function(i) {
    stats::optim(starts[i, ], objective, method = "L-BFGS-B", lower = lower,
                 upper = upper, control = control)
  })
  objective_values <- vapply(fits, `[[`, numeric(1), "value")
  best <- fits[[which.min(objective_values)]]
  names(best$par) <- pathway$rate_classes
  predicted <- panel_composition(predict_composition(pathway, best$par, secretion, entry), names(observed))
  audit <- check_panel(pathway, names(observed), best$par, secretion, entry)
  structure(
    list(pathway = pathway, observed = observed, predicted = predicted,
         log_rates = best$par, rates = exp(best$par), objective = best$value,
         convergence = best$convergence, message = best$message %||% "",
         fit_correlation = .safe_correlation(observed, predicted), loss = loss,
         panel = names(observed), audit = audit, all_objectives = objective_values,
         conditional = TRUE),
    class = "glyco_pathway_fit"
  )
}

#' Compare fitted pathway specifications
#'
#' @param pathways Named list of `glyco_pathway` objects.
#' @param ... Additional arguments passed to [fit_pathway()].
#' @inheritParams fit_pathway
#'
#' @return A `glyco_pathway_comparison` object.
#' @export
compare_pathways <- function(pathways, observed, panel = NULL, n_starts = 8L,
                             seed = 1L, loss = c("composition", "log_composition"), ...) {
  if (!is.list(pathways) || !length(pathways)) .gpd_stop("pathways must be a non-empty list.")
  if (is.null(names(pathways))) names(pathways) <- paste0("model_", seq_along(pathways))
  loss <- match.arg(loss)
  fits <- lapply(seq_along(pathways), function(i) {
    fit_pathway(pathways[[i]], observed, panel, n_starts = n_starts,
                seed = seed + i - 1L, loss = loss, ...)
  })
  names(fits) <- names(pathways)
  metrics <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    data.frame(model = name, rate_classes = length(fit$pathway$rate_classes),
               objective = fit$objective, fit_correlation = fit$fit_correlation,
               rank = fit$audit$rank, full_rank = fit$audit$full_rank,
               condition_number = fit$audit$condition_number,
               amplification = fit$audit$amplification, stringsAsFactors = FALSE)
  }))
  structure(list(fits = fits, metrics = metrics, conditional = TRUE),
            class = "glyco_pathway_comparison")
}

#' @export
print.glyco_pathway_fit <- function(x, ...) {
  cat("<glyco_pathway_fit>\n")
  cat("  panel classes:      ", length(x$panel), "\n", sep = "")
  cat("  fit correlation:    ", .format_finite(x$fit_correlation), "\n", sep = "")
  cat("  objective:          ", .format_finite(x$objective), "\n", sep = "")
  cat("  local rank at fit:  ", x$audit$rank, "/", x$audit$n_rate_classes, "\n", sep = "")
  cat("  condition number:   ", .format_finite(x$audit$condition_number), "\n", sep = "")
  cat("  conditional on the specified forward model and fitted local operating point.\n")
  invisible(x)
}

#' @export
print.glyco_pathway_comparison <- function(x, ...) {
  cat("<glyco_pathway_comparison>\n")
  print(x$metrics, row.names = FALSE)
  cat("  Fit and recoverability are separate model-design criteria.\n")
  invisible(x)
}
