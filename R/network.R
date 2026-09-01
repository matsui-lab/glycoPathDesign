#' Construct a glycoform biosynthetic pathway model
#'
#' @param nodes Data frame or CSV path with columns `id` and optionally `label`.
#' @param edges Data frame or CSV path with columns `from`, `to`, `rate_class`,
#'   and optionally `weight`.
#' @param observations Optional data frame or CSV path with columns `node` and
#'   `glycoform`. Multiple latent nodes may map to one measured glycoform.
#' @param entry Entry node identifier, or a named non-negative vector over nodes.
#' @param secretion Positive scalar, named vector, or vector of length equal to
#'   the number of nodes.
#' @param metadata Optional named list.
#'
#' @return An object of class `glyco_pathway`.
#' @export
glyco_pathway <- function(nodes, edges, observations = NULL, entry = NULL,
                          secretion = 0.5, metadata = list()) {
  nodes <- .as_data_frame(nodes, "nodes")
  edges <- .as_data_frame(edges, "edges")
  .require_columns(nodes, "id", "nodes")
  .require_columns(edges, c("from", "to", "rate_class"), "edges")
  nodes$id <- as.character(nodes$id)
  if (!"label" %in% names(nodes)) nodes$label <- nodes$id
  nodes$label <- as.character(nodes$label)
  if (!length(nodes$id) || anyNA(nodes$id) || any(!nzchar(nodes$id))) {
    .gpd_stop("nodes$id must contain non-empty identifiers.")
  }
  if (anyDuplicated(nodes$id)) .gpd_stop("nodes$id must be unique.")

  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$rate_class <- as.character(edges$rate_class)
  if (!"weight" %in% names(edges)) edges$weight <- 1
  edges$weight <- as.numeric(edges$weight)
  unknown <- setdiff(unique(c(edges$from, edges$to)), nodes$id)
  if (length(unknown)) .gpd_stop("edges reference unknown node(s): %s.", paste(unknown, collapse = ", "))
  if (anyNA(edges$rate_class) || any(!nzchar(edges$rate_class))) {
    .gpd_stop("Every edge must have a non-empty rate_class.")
  }
  if (any(!is.finite(edges$weight)) || any(edges$weight <= 0)) {
    .gpd_stop("Edge weights must be finite and positive.")
  }

  if (is.null(observations)) {
    observations <- data.frame(node = nodes$id, glycoform = nodes$id, stringsAsFactors = FALSE)
  } else {
    observations <- .as_data_frame(observations, "observations")
    .require_columns(observations, c("node", "glycoform"), "observations")
    observations$node <- as.character(observations$node)
    observations$glycoform <- as.character(observations$glycoform)
    unknown_obs <- setdiff(observations$node, nodes$id)
    if (length(unknown_obs)) {
      .gpd_stop("observations reference unknown node(s): %s.", paste(unknown_obs, collapse = ", "))
    }
    if (anyDuplicated(observations$node)) .gpd_stop("Each latent node may occur only once in observations.")
    if (anyNA(observations$glycoform) || any(!nzchar(observations$glycoform))) {
      .gpd_stop("observations$glycoform must contain non-empty identifiers.")
    }
  }

  node_ids <- nodes$id
  if (is.null(entry)) entry <- node_ids[[1L]]
  if (is.character(entry) && length(entry) == 1L) {
    if (!entry %in% node_ids) .gpd_stop("Entry node '%s' is not present in nodes.", entry)
    entry_vector <- setNames(numeric(length(node_ids)), node_ids)
    entry_vector[[entry]] <- 1
  } else {
    entry_vector <- .named_numeric(entry, node_ids, "entry")
    entry_vector <- setNames(.normalise(entry_vector, "entry"), node_ids)
  }

  if (length(secretion) == 1L && is.null(names(secretion))) secretion <- rep(secretion, length(node_ids))
  secretion <- .named_numeric(secretion, node_ids, "secretion")
  if (any(secretion <= 0)) .gpd_stop("secretion must be positive for every node.")

  rate_classes <- sort(unique(edges$rate_class))
  observed_classes <- unique(observations$glycoform)
  object <- list(
    nodes = nodes,
    edges = edges,
    observations = observations,
    entry = entry_vector,
    secretion = setNames(secretion, node_ids),
    rate_classes = rate_classes,
    observed_classes = observed_classes,
    metadata = metadata
  )
  class(object) <- "glyco_pathway"
  validate_pathway(object)
  object
}

#' Validate a glycoform biosynthetic pathway model
#'
#' @param pathway A `glyco_pathway` object.
#' @return The input invisibly, or an error when invalid.
#' @export
validate_pathway <- function(pathway) {
  if (!inherits(pathway, "glyco_pathway")) .gpd_stop("pathway must inherit from 'glyco_pathway'.")
  required <- c("nodes", "edges", "observations", "entry", "secretion", "rate_classes", "observed_classes")
  missing <- setdiff(required, names(pathway))
  if (length(missing)) .gpd_stop("pathway is missing component(s): %s.", paste(missing, collapse = ", "))
  if (length(pathway$rate_classes) < 1L) .gpd_stop("The pathway has no rate classes.")
  if (length(pathway$observed_classes) < 2L) .gpd_stop("At least two observed glycoform classes are required.")
  invisible(pathway)
}

#' Read a pathway from CSV files
#'
#' @inheritParams glyco_pathway
#' @return A `glyco_pathway` object.
#' @export
read_pathway <- function(nodes, edges, observations = NULL, entry = NULL,
                         secretion = 0.5, metadata = list()) {
  glyco_pathway(nodes, edges, observations, entry, secretion, metadata)
}

#' Read a pathway directory written by `write_pathway()`
#'
#' @param directory Directory containing `nodes.csv`, `edges.csv`,
#'   `observations.csv`, and `boundary.csv`.
#' @param metadata Optional named list attached to the pathway.
#'
#' @return A `glyco_pathway` object.
#' @export
read_pathway_directory <- function(directory, metadata = list()) {
  files <- file.path(directory, c("nodes.csv", "edges.csv", "observations.csv", "boundary.csv"))
  missing <- files[!file.exists(files)]
  if (length(missing)) .gpd_stop("Pathway directory is missing: %s.", paste(basename(missing), collapse = ", "))
  boundary <- utils::read.csv(files[[4L]], stringsAsFactors = FALSE, check.names = FALSE)
  .require_columns(boundary, c("node", "entry", "secretion"), "boundary.csv")
  entry <- setNames(boundary$entry, boundary$node)
  secretion <- setNames(boundary$secretion, boundary$node)
  glyco_pathway(files[[1L]], files[[2L]], files[[3L]], entry, secretion, metadata)
}

#' Write a pathway definition to a directory
#'
#' @param pathway A `glyco_pathway` object.
#' @param directory Output directory.
#' @return The normalised directory path invisibly.
#' @export
write_pathway <- function(pathway, directory) {
  validate_pathway(pathway)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(pathway$nodes, file.path(directory, "nodes.csv"), row.names = FALSE)
  utils::write.csv(pathway$edges, file.path(directory, "edges.csv"), row.names = FALSE)
  utils::write.csv(pathway$observations, file.path(directory, "observations.csv"), row.names = FALSE)
  utils::write.csv(data.frame(node = names(pathway$entry), entry = unname(pathway$entry),
                              secretion = unname(pathway$secretion)),
                   file.path(directory, "boundary.csv"), row.names = FALSE)
  invisible(normalizePath(directory, mustWork = TRUE))
}

#' Canonical eight-rate-class N-glycan pathway
#'
#' Returns the 22-state, 25-edge illustrative pathway used in the accompanying
#' manuscript. It has Man9 entry, uniform secretion of 0.5, and identity
#' observation mapping. It is a design example, not a universal glycosylation
#' model.
#'
#' @return A `glyco_pathway` object.
#' @export
canonical_pathway <- function() {
  node_ids <- c("M9", "M8", "M7", "M6", "M5", "M5Gn", "M4Gn", "M3Gn",
                "M3Gn2", "M3Gn3", "G1", "G2", "G2S1", "G2S2", "G1S1",
                "M3Gn2F", "G1F", "G2F", "G2S1F", "G2S2F", "Gt2", "Gt2S")
  edge_values <- rbind(
    c("M9", "M8", "MAN1"), c("M8", "M7", "MAN1"), c("M7", "M6", "MAN1"),
    c("M6", "M5", "MAN1"), c("M5", "M5Gn", "MGAT1"),
    c("M5Gn", "M4Gn", "MAN2"), c("M4Gn", "M3Gn", "MAN2"),
    c("M3Gn", "M3Gn2", "MGAT2"), c("M3Gn2", "M3Gn3", "MGAT4"),
    c("M3Gn2", "G1", "B4GALT"), c("G1", "G2", "B4GALT"),
    c("M3Gn3", "Gt2", "B4GALT"), c("G2", "G2S1", "ST6GAL"),
    c("G2S1", "G2S2", "ST6GAL"), c("G1", "G1S1", "ST6GAL"),
    c("Gt2", "Gt2S", "ST6GAL"), c("M3Gn2", "M3Gn2F", "FUT8"),
    c("G1", "G1F", "FUT8"), c("G2", "G2F", "FUT8"),
    c("G2S1", "G2S1F", "FUT8"), c("G2S2", "G2S2F", "FUT8"),
    c("M3Gn2F", "G1F", "B4GALT"), c("G1F", "G2F", "B4GALT"),
    c("G2F", "G2S1F", "ST6GAL"), c("G2S1F", "G2S2F", "ST6GAL")
  )
  edges <- data.frame(from = edge_values[, 1], to = edge_values[, 2],
                      rate_class = edge_values[, 3], stringsAsFactors = FALSE)
  glyco_pathway(
    nodes = data.frame(id = node_ids, label = node_ids, stringsAsFactors = FALSE),
    edges = edges,
    entry = "M9",
    secretion = 0.5,
    metadata = list(name = "Canonical N-glycan pathway", conditional = TRUE,
                    note = "Illustrative 22-state, eight-rate-class design model")
  )
}

#' @export
print.glyco_pathway <- function(x, ...) {
  cat("<glyco_pathway>\n")
  cat("  latent states:   ", nrow(x$nodes), "\n", sep = "")
  cat("  observed classes:", length(x$observed_classes), "\n")
  cat("  edges:           ", nrow(x$edges), "\n", sep = "")
  cat("  rate classes:    ", length(x$rate_classes), " (",
      paste(x$rate_classes, collapse = ", "), ")\n", sep = "")
  cat("  entry:           ", paste(names(x$entry)[x$entry > 0], collapse = ", "), "\n", sep = "")
  invisible(x)
}
