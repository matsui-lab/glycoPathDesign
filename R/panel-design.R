#' Enumerate or sample candidate glycoform panels
#'
#' Candidate panels are evaluated at one specified operating point. If the
#' number of combinations exceeds `max_panels`, a reproducible random sample is
#' evaluated instead of silently truncating lexicographic enumeration.
#'
#' @param pathway A `glyco_pathway` object.
#' @param panel_size Number of measured glycoform classes.
#' @param candidates Candidate observed classes. Defaults to all classes.
#' @param must_include Classes required in every panel.
#' @param max_panels Maximum number of panels evaluated.
#' @param seed Random seed used when panels are sampled.
#' @inheritParams check_panel
#'
#' @return A `glyco_panel_set` object.
#' @export
enumerate_panels <- function(pathway, panel_size, candidates = NULL,
                             must_include = character(), max_panels = 100000L,
                             seed = 1L, log_rates = NULL, secretion = NULL,
                             entry = NULL, target_sd = 0.25, assay_cv = NULL,
                             rank_tol = 1e-9, rcond = 1e-10, eps = 1e-6,
                             rate_classes = NULL) {
  validate_pathway(pathway)
  candidates <- as.character(candidates %||% pathway$observed_classes)
  must_include <- unique(as.character(must_include))
  unknown <- setdiff(unique(c(candidates, must_include)), pathway$observed_classes)
  if (length(unknown)) .gpd_stop("Unknown observed class(es): %s.", paste(unknown, collapse = ", "))
  candidates <- unique(c(must_include, candidates))
  panel_size <- as.integer(panel_size)
  max_panels <- as.integer(max_panels)
  if (length(panel_size) != 1L || is.na(panel_size) || panel_size < 2L || panel_size > length(candidates)) {
    .gpd_stop("panel_size must be between 2 and the number of candidates.")
  }
  if (length(must_include) > panel_size) .gpd_stop("must_include is larger than panel_size.")
  if (max_panels < 1L) .gpd_stop("max_panels must be positive.")

  optional <- setdiff(candidates, must_include)
  choose_n <- panel_size - length(must_include)
  total <- choose(length(optional), choose_n)
  if (!is.finite(total)) .gpd_stop("The number of candidate combinations is too large.")
  sampled <- total > max_panels
  n_eval <- as.integer(min(total, max_panels))
  if (!sampled) {
    if (choose_n == 0L) {
      panels <- matrix(must_include, nrow = length(must_include), ncol = 1L)
    } else {
      panels <- utils::combn(optional, choose_n)
      if (length(must_include)) {
        panels <- rbind(matrix(rep(must_include, ncol(panels)), nrow = length(must_include)), panels)
      }
    }
  } else {
    set.seed(seed)
    keys <- character()
    panels_list <- vector("list", n_eval)
    i <- 0L
    attempts <- 0L
    max_attempts <- max(1000L, 50L * n_eval)
    while (i < n_eval && attempts < max_attempts) {
      attempts <- attempts + 1L
      draw <- sort(sample(optional, choose_n, replace = FALSE))
      key <- paste(draw, collapse = "\r")
      if (!key %in% keys) {
        i <- i + 1L
        keys[[i]] <- key
        panels_list[[i]] <- c(must_include, draw)
      }
    }
    if (i < n_eval) .gpd_stop("Could not obtain the requested number of unique sampled panels.")
    panels <- do.call(cbind, panels_list)
  }

  full_jacobian <- lograte_jacobian(pathway, log_rates, secretion, entry, eps)
  full_composition <- attr(full_jacobian, "composition")
  log_rates <- attr(full_jacobian, "log_rates")
  rows <- vector("list", ncol(panels))
  audits <- vector("list", ncol(panels))
  for (i in seq_len(ncol(panels))) {
    panel <- panels[, i]
    audit <- .audit_from_context(pathway, panel, full_jacobian, full_composition,
                                 target_sd, assay_cv, rank_tol, rcond, log_rates,
                                 rate_classes)
    audits[[i]] <- audit
    rows[[i]] <- data.frame(
      panel_id = i,
      panel = paste(panel, collapse = " | "),
      panel_size = length(panel),
      rank = audit$rank,
      full_rank = audit$full_rank,
      condition_number = audit$condition_number,
      amplification = audit$amplification,
      max_assay_cv = audit$max_assay_cv,
      precision_target_met = audit$precision_target_met,
      stringsAsFactors = FALSE
    )
  }
  metrics <- do.call(rbind, rows)
  structure(
    list(pathway = pathway, panels = panels, metrics = metrics, audits = audits,
         total_combinations = total, evaluated = ncol(panels), sampled = sampled,
         target_sd = target_sd, assay_cv = assay_cv, conditional = TRUE),
    class = "glyco_panel_set"
  )
}

#' Rank candidate panels by recoverability
#'
#' @param panel_set A `glyco_panel_set` object.
#' @param objective Rank by noise amplification or condition number.
#' @param full_rank_only Retain only locally full-rank panels.
#' @param n Number of panels returned.
#'
#' @return A data frame ordered from best to worst.
#' @export
optimise_panels <- function(panel_set, objective = c("amplification", "condition_number"),
                            full_rank_only = TRUE, n = 20L) {
  if (!inherits(panel_set, "glyco_panel_set")) .gpd_stop("panel_set must be a glyco_panel_set object.")
  objective <- match.arg(objective)
  out <- panel_set$metrics
  if (full_rank_only) out <- out[out$full_rank, , drop = FALSE]
  out <- out[order(out[[objective]], out$condition_number, na.last = TRUE), , drop = FALSE]
  utils::head(out, as.integer(n))
}

#' Find the first searched panel size with a locally full-rank candidate
#'
#' @param sizes Panel sizes to examine.
#' @param stop_at_first Stop after finding the first qualifying size.
#' @inheritParams enumerate_panels
#'
#' @return A list containing the size summary and evaluated panel sets. The
#'   `minimality_proven` field is `TRUE` only when every smaller candidate size
#'   was exhaustively enumerated; a sampled negative result is never treated as
#'   a proof of minimality.
#' @export
minimum_panel_size <- function(pathway, sizes = 2:length(pathway$observed_classes),
                               candidates = NULL, must_include = character(),
                               max_panels = 100000L, seed = 1L, log_rates = NULL,
                               secretion = NULL, entry = NULL, target_sd = 0.25,
                               assay_cv = NULL, rank_tol = 1e-9, rcond = 1e-10,
                               eps = 1e-6, stop_at_first = TRUE,
                               rate_classes = NULL) {
  sets <- list()
  rows <- list()
  for (size in unique(as.integer(sizes))) {
    current <- enumerate_panels(pathway, size, candidates, must_include, max_panels,
                                seed, log_rates, secretion, entry, target_sd,
                                assay_cv, rank_tol, rcond, eps, rate_classes)
    sets[[as.character(size)]] <- current
    rows[[length(rows) + 1L]] <- data.frame(
      panel_size = size,
      total_combinations = current$total_combinations,
      evaluated = current$evaluated,
      sampled = current$sampled,
      full_rank_panels = sum(current$metrics$full_rank),
      best_amplification = suppressWarnings(min(current$metrics$amplification, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    if (stop_at_first && any(current$metrics$full_rank)) break
  }
  summary <- do.call(rbind, rows)
  summary$best_amplification[!is.finite(summary$best_amplification)] <- Inf
  found_size <- if (any(summary$full_rank_panels > 0)) {
    min(summary$panel_size[summary$full_rank_panels > 0])
  } else NA_integer_
  prior <- if (is.na(found_size)) seq_len(nrow(summary)) else which(summary$panel_size <= found_size)
  proven <- !is.na(found_size) && all(!summary$sampled[prior]) &&
    all(summary$full_rank_panels[summary$panel_size < found_size] == 0)
  structure(list(summary = summary, panel_sets = sets,
                 minimum_size = found_size, minimality_proven = proven,
                 search_interpretation = if (proven) {
                   "Minimum size proven within the supplied candidate set and assumptions."
                 } else if (!is.na(found_size)) {
                   "First full-rank size found; smaller sampled searches do not prove minimality."
                 } else {
                   "No full-rank panel was found in the evaluated searches."
                 }),
            class = "glyco_minimum_panel")
}

#' @export
print.glyco_panel_set <- function(x, ...) {
  cat("<glyco_panel_set>\n")
  cat("  combinations:      ", format(x$total_combinations, scientific = FALSE), "\n", sep = "")
  cat("  evaluated:         ", x$evaluated, if (x$sampled) " (random sample)" else " (exhaustive)", "\n", sep = "")
  cat("  full-rank panels:  ", sum(x$metrics$full_rank), "\n", sep = "")
  best <- optimise_panels(x, n = 1L)
  if (nrow(best)) cat("  best amplification:", .format_finite(best$amplification), "\n")
  cat("  conditional on the supplied graph and operating point.\n")
  invisible(x)
}

#' @export
print.glyco_minimum_panel <- function(x, ...) {
  cat("<glyco_minimum_panel>\n")
  print(x$summary, row.names = FALSE)
  cat("  first full-rank size:", ifelse(is.na(x$minimum_size), "not found", x$minimum_size), "\n")
  cat("  interpretation:      ", x$search_interpretation, "\n", sep = "")
  invisible(x)
}
