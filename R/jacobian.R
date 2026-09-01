#' Compute the local log-rate Jacobian
#'
#' @inheritParams forward_pathway
#' @param eps Finite-difference step on the log-rate scale.
#' @param method Central or forward finite differences.
#'
#' @return A matrix with observed glycoforms in rows and rate classes in columns.
#' @export
lograte_jacobian <- function(pathway, log_rates = NULL, secretion = NULL, entry = NULL,
                             eps = 1e-6, method = c("central", "forward")) {
  validate_pathway(pathway)
  method <- match.arg(method)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    .gpd_stop("eps must be a positive finite scalar.")
  }
  log_rates <- .named_numeric(
    log_rates, pathway$rate_classes, "log_rates",
    default = setNames(rep(0, length(pathway$rate_classes)), pathway$rate_classes)
  )
  base <- predict_composition(pathway, log_rates, secretion, entry)
  jac <- matrix(NA_real_, nrow = length(base), ncol = length(log_rates),
                dimnames = list(names(base), names(log_rates)))
  for (j in seq_along(log_rates)) {
    plus <- log_rates
    plus[[j]] <- plus[[j]] + eps
    if (method == "central") {
      minus <- log_rates
      minus[[j]] <- minus[[j]] - eps
      jac[, j] <- (predict_composition(pathway, plus, secretion, entry) -
                    predict_composition(pathway, minus, secretion, entry)) / (2 * eps)
    } else {
      jac[, j] <- (predict_composition(pathway, plus, secretion, entry) - base) / eps
    }
  }
  attr(jac, "composition") <- base
  attr(jac, "log_rates") <- log_rates
  jac
}

.panel_jacobian_from_full <- function(full_jacobian, full_composition, panel) {
  panel <- as.character(panel)
  missing <- setdiff(panel, rownames(full_jacobian))
  if (length(missing)) .gpd_stop("Panel contains unknown glycoform(s): %s.", paste(missing, collapse = ", "))
  mass <- sum(full_composition[panel])
  if (!is.finite(mass) || mass <= 0) .gpd_stop("The selected panel has zero predicted mass.")
  selected_j <- full_jacobian[panel, , drop = FALSE]
  selected_mu <- full_composition[panel]
  (selected_j * mass - outer(selected_mu, colSums(selected_j))) / (mass^2)
}

#' Compute the within-panel-normalised log-rate Jacobian
#'
#' @inheritParams lograte_jacobian
#' @param panel Character or integer measured-glycoform panel.
#'
#' @return A panel-by-rate-class Jacobian matrix.
#' @export
panel_jacobian <- function(pathway, panel, log_rates = NULL, secretion = NULL, entry = NULL,
                           eps = 1e-6, method = c("central", "forward")) {
  method <- match.arg(method)
  full <- lograte_jacobian(pathway, log_rates, secretion, entry, eps, method)
  composition <- attr(full, "composition")
  if (is.numeric(panel)) panel <- names(composition)[panel]
  out <- .panel_jacobian_from_full(full, composition, panel)
  attr(out, "composition") <- panel_composition(composition, panel)
  attr(out, "full_composition") <- composition
  attr(out, "log_rates") <- attr(full, "log_rates")
  out
}
