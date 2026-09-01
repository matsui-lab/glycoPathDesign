.noise_amplification <- function(jacobian, composition, rank_tol = 1e-9, rcond = 1e-10) {
  k <- ncol(jacobian)
  if (.matrix_rank(jacobian, rank_tol) < k) return(Inf)
  y <- .normalise(composition, "composition")
  m <- diag(y, nrow = length(y)) %*% (diag(length(y)) - outer(rep(1, length(y)), y))
  decomposition <- svd(m)
  keep <- decomposition$d > rcond * max(decomposition$d)
  if (!any(keep)) return(Inf)
  whitening <- sweep(t(decomposition$u[, keep, drop = FALSE]), 1L,
                     decomposition$d[keep], "/")
  whitened <- whitening %*% jacobian
  values <- svd(whitened, nu = 0, nv = 0)$d
  if (length(values) < k || values[[k]] <= rcond * values[[1L]]) return(Inf)
  1 / values[[k]]
}

.audit_from_context <- function(pathway, panel, full_jacobian, full_composition,
                                target_sd, assay_cv, rank_tol, rcond, log_rates,
                                rate_classes = NULL) {
  if (is.numeric(panel)) panel <- names(full_composition)[panel]
  panel <- as.character(panel)
  jacobian <- .panel_jacobian_from_full(full_jacobian, full_composition, panel)
  rate_classes <- as.character(rate_classes %||% colnames(jacobian))
  unknown_rates <- setdiff(rate_classes, colnames(jacobian))
  if (length(unknown_rates)) {
    .gpd_stop("Unknown rate class(es): %s.", paste(unknown_rates, collapse = ", "))
  }
  if (!length(rate_classes)) .gpd_stop("At least one rate class must be selected.")
  jacobian <- jacobian[, rate_classes, drop = FALSE]
  composition <- panel_composition(full_composition, panel)
  singular <- svd(jacobian, nu = 0, nv = 0)$d
  k <- ncol(jacobian)
  rank <- as.integer(sum(singular > rank_tol))
  full_rank <- rank == k
  sigma_max <- if (length(singular)) singular[[1L]] else NA_real_
  sigma_min <- if (full_rank) singular[[k]] else 0
  condition <- if (full_rank) sigma_max / sigma_min else Inf
  amplification <- .noise_amplification(jacobian, composition, rank_tol, rcond)
  max_cv <- if (is.finite(amplification)) target_sd / amplification else 0
  achieved_sd <- if (is.null(assay_cv)) NA_real_ else amplification * assay_cv
  structure(
    list(
      pathway = pathway,
      panel = panel,
      panel_size = length(panel),
      rate_classes = colnames(jacobian),
      n_rate_classes = k,
      independent_dimensions = max(0L, length(panel) - 1L),
      dimension_sufficient = length(panel) - 1L >= k,
      rank = rank,
      full_rank = full_rank,
      singular_values = singular,
      sigma_max = sigma_max,
      sigma_min = sigma_min,
      condition_number = condition,
      amplification = amplification,
      target_sd = target_sd,
      max_assay_cv = max_cv,
      assay_cv = assay_cv,
      achieved_sd = achieved_sd,
      precision_target_met = if (is.null(assay_cv)) NA else is.finite(achieved_sd) && achieved_sd <= target_sd,
      composition = composition,
      jacobian = jacobian,
      log_rates = log_rates,
      rank_tolerance = rank_tol,
      rcond = rcond,
      conditional = TRUE
    ),
    class = "glyco_panel_audit"
  )
}

#' Evaluate a proposed glycoform measurement panel
#'
#' Checks the dimension bound, local rank, singular spectrum, condition number,
#' and coefficient-of-variation requirement for recovering the specified local
#' log-rate directions.
#'
#' @inheritParams panel_jacobian
#' @param target_sd Target worst-direction log-rate standard deviation.
#' @param assay_cv Optional assumed coefficient of variation on class intensities.
#' @param rank_tol Absolute singular-value threshold for rank.
#' @param rcond Relative tolerance used while whitening compositional noise.
#' @param rate_classes Effective reaction-rate classes to distinguish. Defaults
#'   to every rate class in the supplied pathway. Unselected classes remain in
#'   the forward model but are held fixed in the local inverse problem.
#'
#' @return A `glyco_panel_audit` object.
#' @export
check_panel <- function(pathway, panel, log_rates = NULL, secretion = NULL, entry = NULL,
                        target_sd = 0.25, assay_cv = NULL, rank_tol = 1e-9,
                        rcond = 1e-10, eps = 1e-6,
                        method = c("central", "forward"), rate_classes = NULL) {
  method <- match.arg(method)
  if (!is.numeric(target_sd) || length(target_sd) != 1L || target_sd <= 0) {
    .gpd_stop("target_sd must be positive.")
  }
  if (!is.null(assay_cv) && (!is.numeric(assay_cv) || length(assay_cv) != 1L || assay_cv < 0)) {
    .gpd_stop("assay_cv must be a non-negative scalar.")
  }
  full <- lograte_jacobian(pathway, log_rates, secretion, entry, eps, method)
  composition <- attr(full, "composition")
  log_rates <- attr(full, "log_rates")
  .audit_from_context(pathway, panel, full, composition, target_sd, assay_cv,
                      rank_tol, rcond, log_rates, rate_classes)
}

#' Audit one panel across multiple operating points
#'
#' @param pathway A `glyco_pathway` object.
#' @param panel Character measured-glycoform panel.
#' @param operating_points Matrix or data frame with one row per log-rate vector
#'   and one column per rate class. If omitted, points are sampled from a normal
#'   distribution around zero.
#' @param n Number of sampled points when `operating_points` is omitted.
#' @param log_rate_sd Sampling standard deviation.
#' @param seed Random seed used only for sampling.
#' @inheritParams check_panel
#'
#' @return A `glyco_robustness_audit` object.
#' @export
check_operating_points <- function(pathway, panel, operating_points = NULL, n = 100L,
                                   log_rate_sd = 0.7, seed = 1L, target_sd = 0.25,
                                   assay_cv = NULL, rank_tol = 1e-9, rcond = 1e-10,
                                   eps = 1e-6, rate_classes = NULL) {
  validate_pathway(pathway)
  k <- length(pathway$rate_classes)
  if (is.null(operating_points)) {
    set.seed(seed)
    operating_points <- matrix(stats::rnorm(as.integer(n) * k, sd = log_rate_sd), ncol = k)
    colnames(operating_points) <- pathway$rate_classes
  } else {
    operating_points <- as.matrix(operating_points)
    if (is.null(colnames(operating_points)) && ncol(operating_points) == k) {
      colnames(operating_points) <- pathway$rate_classes
    }
    missing <- setdiff(pathway$rate_classes, colnames(operating_points))
    if (length(missing)) .gpd_stop("operating_points is missing rate class(es): %s.", paste(missing, collapse = ", "))
    operating_points <- operating_points[, pathway$rate_classes, drop = FALSE]
  }
  rows <- vector("list", nrow(operating_points))
  audits <- vector("list", nrow(operating_points))
  for (i in seq_len(nrow(operating_points))) {
    audit <- check_panel(pathway, panel, operating_points[i, ], target_sd = target_sd,
                         assay_cv = assay_cv, rank_tol = rank_tol, rcond = rcond,
                         eps = eps, rate_classes = rate_classes)
    audits[[i]] <- audit
    rows[[i]] <- data.frame(
      operating_point = i,
      rank = audit$rank,
      full_rank = audit$full_rank,
      condition_number = audit$condition_number,
      amplification = audit$amplification,
      max_assay_cv = audit$max_assay_cv,
      stringsAsFactors = FALSE
    )
  }
  metrics <- do.call(rbind, rows)
  finite_condition <- metrics$condition_number[is.finite(metrics$condition_number)]
  finite_g <- metrics$amplification[is.finite(metrics$amplification)]
  structure(
    list(
      pathway = pathway,
      panel = panel,
      operating_points = operating_points,
      audits = audits,
      metrics = metrics,
      full_rank_fraction = mean(metrics$full_rank),
      condition_median = if (length(finite_condition)) stats::median(finite_condition) else Inf,
      condition_max = if (length(finite_condition)) max(finite_condition) else Inf,
      amplification_median = if (length(finite_g)) stats::median(finite_g) else Inf,
      max_cv_tenth_percentile = if (any(metrics$max_assay_cv > 0)) {
        as.numeric(stats::quantile(metrics$max_assay_cv, 0.1, names = FALSE))
      } else 0,
      conditional = TRUE
    ),
    class = "glyco_robustness_audit"
  )
}

#' Determine the locally supported biological resolution
#'
#' @param audit A `glyco_panel_audit` object.
#' @param fit_correlation Optional correlation between fitted and observed
#'   compositions.
#' @param minimum_fit Optional user-declared minimum acceptable fit correlation.
#' @param condition_warning User-declared condition-number warning threshold.
#'
#' @return A one-row data frame with a machine-readable code and interpretation.
#' @export
supported_resolution <- function(audit, fit_correlation = NULL, minimum_fit = NULL,
                                 condition_warning = 100) {
  if (!inherits(audit, "glyco_panel_audit")) .gpd_stop("audit must be a glyco_panel_audit object.")
  code <- "local_rate_directions_supported"
  label <- "Individual rate-class directions are locally supported under the stated assumptions."
  if (!audit$dimension_sufficient) {
    code <- "dimensionally_insufficient"
    label <- "The panel has fewer independent composition directions than rate classes."
  } else if (!audit$full_rank) {
    code <- "rank_deficient"
    label <- "Some rate-class directions are locally indistinguishable in this panel."
  } else if (!is.null(minimum_fit) && !is.null(fit_correlation) &&
             (is.na(fit_correlation) || fit_correlation < minimum_fit)) {
    code <- "forward_model_inadequate"
    label <- "The declared forward-model fit criterion is not met."
  } else if (!is.null(audit$assay_cv) && !isTRUE(audit$precision_target_met)) {
    code <- "precision_target_not_met"
    label <- "The panel is full rank, but the assumed assay precision does not meet the target."
  } else if (is.finite(audit$condition_number) && audit$condition_number > condition_warning) {
    code <- "full_rank_noise_sensitive"
    label <- "The panel is full rank but strongly amplifies noise in weakly observed directions."
  }
  data.frame(
    code = code,
    interpretation = label,
    panel_size = audit$panel_size,
    rate_classes = audit$n_rate_classes,
    rank = audit$rank,
    condition_number = audit$condition_number,
    amplification = audit$amplification,
    max_assay_cv = audit$max_assay_cv,
    fit_correlation = fit_correlation %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @export
print.glyco_panel_audit <- function(x, ...) {
  cat("<glyco_panel_audit>\n")
  cat("  panel:              ", paste(x$panel, collapse = ", "), "\n", sep = "")
  cat("  dimension:          ", x$independent_dimensions, " independent / ",
      x$n_rate_classes, " required\n", sep = "")
  cat("  local rank:         ", x$rank, "/", x$n_rate_classes, "\n", sep = "")
  cat("  condition number:   ", .format_finite(x$condition_number), "\n", sep = "")
  cat("  amplification g:    ", .format_finite(x$amplification), "\n", sep = "")
  cat("  maximum assay CV:   ", formatC(100 * x$max_assay_cv, digits = 3, format = "fg"),
      "% for target SD ", x$target_sd, "\n", sep = "")
  cat("  interpretation:     ", supported_resolution(x)$interpretation, "\n", sep = "")
  cat("  conditional on the supplied graph, operating point, observation mapping, and error model.\n")
  invisible(x)
}

#' @export
print.glyco_robustness_audit <- function(x, ...) {
  cat("<glyco_robustness_audit>\n")
  cat("  operating points:   ", nrow(x$metrics), "\n", sep = "")
  cat("  full-rank fraction: ", formatC(x$full_rank_fraction, digits = 3, format = "f"), "\n", sep = "")
  cat("  median condition:   ", .format_finite(x$condition_median), "\n", sep = "")
  cat("  median g:           ", .format_finite(x$amplification_median), "\n", sep = "")
  cat("  conditional on the sampled operating-point distribution.\n")
  invisible(x)
}
