.rate_matrix <- function(pathway, log_rates) {
  node_ids <- pathway$nodes$id
  rates <- exp(.named_numeric(log_rates, pathway$rate_classes, "log_rates",
                              default = setNames(rep(0, length(pathway$rate_classes)), pathway$rate_classes)))
  matrix_out <- matrix(0, length(node_ids), length(node_ids), dimnames = list(node_ids, node_ids))
  from <- match(pathway$edges$from, node_ids)
  to <- match(pathway$edges$to, node_ids)
  values <- rates[pathway$edges$rate_class] * pathway$edges$weight
  for (i in seq_along(values)) matrix_out[from[[i]], to[[i]]] <- matrix_out[from[[i]], to[[i]]] + values[[i]]
  matrix_out
}

.observation_matrix <- function(pathway) {
  node_ids <- pathway$nodes$id
  classes <- pathway$observed_classes
  out <- matrix(0, nrow = length(classes), ncol = length(node_ids),
                dimnames = list(classes, node_ids))
  node_index <- match(pathway$observations$node, node_ids)
  class_index <- match(pathway$observations$glycoform, classes)
  out[cbind(class_index, node_index)] <- 1
  out
}

#' Evaluate the steady-state pathway forward model
#'
#' @param pathway A `glyco_pathway` object.
#' @param log_rates Named log effective rates. Unnamed vectors are accepted when
#'   their order matches `pathway$rate_classes`.
#' @param secretion Optional positive node-specific secretion vector. Defaults
#'   to the pathway definition.
#' @param entry Optional entry distribution. Defaults to the pathway definition.
#'
#' @return A list containing latent occupancy, state secretion flux, observed
#'   flux, and the normalised observed composition.
#' @export
forward_pathway <- function(pathway, log_rates = NULL, secretion = NULL, entry = NULL) {
  validate_pathway(pathway)
  node_ids <- pathway$nodes$id
  log_rates <- .named_numeric(
    log_rates, pathway$rate_classes, "log_rates",
    default = setNames(rep(0, length(pathway$rate_classes)), pathway$rate_classes)
  )
  secretion <- .named_numeric(secretion, node_ids, "secretion", default = pathway$secretion)
  if (any(secretion <= 0)) .gpd_stop("secretion must be positive.")
  entry <- .named_numeric(entry, node_ids, "entry", default = pathway$entry)
  entry <- setNames(.normalise(entry, "entry"), node_ids)

  reaction <- .rate_matrix(pathway, log_rates)
  generator <- diag(secretion + rowSums(reaction), nrow = length(node_ids)) - t(reaction)
  occupancy <- tryCatch(
    solve(generator, entry),
    error = function(e) .gpd_stop("The forward model could not be solved: %s", conditionMessage(e))
  )
  if (any(!is.finite(occupancy)) || any(occupancy < -1e-10)) {
    .gpd_stop("The forward model produced an inadmissible steady state.")
  }
  occupancy <- pmax(occupancy, 0)
  state_flux <- secretion * occupancy
  observation <- .observation_matrix(pathway)
  observed_flux <- as.vector(observation %*% state_flux)
  names(observed_flux) <- rownames(observation)
  composition <- setNames(.normalise(observed_flux, "observed secretion flux"), names(observed_flux))
  structure(
    list(
      composition = composition,
      observed_flux = observed_flux,
      state_flux = setNames(state_flux, node_ids),
      occupancy = setNames(occupancy, node_ids),
      reaction_matrix = reaction,
      log_rates = log_rates,
      secretion = secretion,
      entry = entry
    ),
    class = "glyco_forward"
  )
}

#' Predict a normalised glycoform composition
#'
#' @inheritParams forward_pathway
#' @return A named numeric composition summing to one.
#' @export
predict_composition <- function(pathway, log_rates = NULL, secretion = NULL, entry = NULL) {
  forward_pathway(pathway, log_rates, secretion, entry)$composition
}

#' Select and renormalise a measured glycoform panel
#'
#' @param composition Named full composition.
#' @param panel Character or integer panel specification.
#' @return Named within-panel-normalised composition.
#' @export
panel_composition <- function(composition, panel) {
  if (is.numeric(panel)) panel <- names(composition)[panel]
  panel <- as.character(panel)
  missing <- setdiff(panel, names(composition))
  if (length(missing)) .gpd_stop("Panel contains unknown glycoform(s): %s.", paste(missing, collapse = ", "))
  if (anyDuplicated(panel)) .gpd_stop("Panel glycoforms must be unique.")
  setNames(.normalise(composition[panel], "panel composition"), panel)
}

#' @export
print.glyco_forward <- function(x, ...) {
  cat("<glyco_forward>\n")
  cat("  observed classes:", length(x$composition), "\n")
  cat("  composition sum:", format(sum(x$composition), digits = 6), "\n")
  cat("  largest class:   ", names(which.max(x$composition)), " (",
      format(max(x$composition), digits = 4), ")\n", sep = "")
  invisible(x)
}
