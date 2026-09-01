.audit_table <- function(audit) {
  resolution <- supported_resolution(audit)
  data.frame(
    metric = c("Panel size", "Independent composition dimensions", "Rate classes",
               "Local rank", "Full rank", "Condition number", "Noise amplification g",
               "Target log-rate SD", "Maximum assay CV", "Supported resolution"),
    value = c(audit$panel_size, audit$independent_dimensions, audit$n_rate_classes,
              audit$rank, audit$full_rank, audit$condition_number, audit$amplification,
              audit$target_sd, audit$max_assay_cv, resolution$interpretation),
    stringsAsFactors = FALSE
  )
}

#' Write a portable panel-design report
#'
#' @param object A `glyco_panel_audit`, `glyco_panel_set`,
#'   `glyco_robustness_audit`, or `glyco_pathway_fit` object.
#' @param file Output `.csv`, `.json`, or `.html` file.
#'
#' @return The normalised output path invisibly.
#' @export
write_design_report <- function(object, file) {
  extension <- tolower(tools::file_ext(file))
  table <- if (inherits(object, "glyco_panel_audit")) {
    .audit_table(object)
  } else if (inherits(object, "glyco_panel_set") || inherits(object, "glyco_robustness_audit")) {
    object$metrics
  } else if (inherits(object, "glyco_pathway_fit")) {
    data.frame(glycoform = names(object$observed), observed = unname(object$observed),
               predicted = unname(object$predicted), stringsAsFactors = FALSE)
  } else {
    .gpd_stop("Unsupported report object class.")
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (extension == "csv") {
    utils::write.csv(table, file, row.names = FALSE)
  } else if (extension == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) .gpd_stop("Install jsonlite to write JSON reports.")
    jsonlite::write_json(table, file, pretty = TRUE, dataframe = "rows", na = "null")
  } else if (extension == "html") {
    escape <- function(x) {
      x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      gsub(">", "&gt;", x, fixed = TRUE)
    }
    header <- paste(sprintf("<th>%s</th>", escape(names(table))), collapse = "")
    body <- apply(table, 1L, function(row) paste(sprintf("<td>%s</td>", escape(row)), collapse = ""))
    html <- c("<!doctype html><html><head><meta charset='utf-8'>",
              "<title>glycoPathDesign report</title>",
              "<style>body{font:16px Arial,sans-serif;max-width:1000px;margin:40px auto;color:#111}table{border-collapse:collapse;width:100%}th,td{border:1px solid #aaa;padding:8px;text-align:left}th{background:#eee}.note{color:#444;margin-top:20px}</style>",
              "</head><body><h1>glycoPathDesign report</h1><table><thead><tr>",
              header, "</tr></thead><tbody>", paste0("<tr>", body, "</tr>"),
              "</tbody></table><p class='note'>All conclusions are conditional on the supplied pathway graph, operating point, observation mapping, and error model.</p></body></html>")
    writeLines(html, file, useBytes = TRUE)
  } else {
    .gpd_stop("file must have extension .csv, .json, or .html.")
  }
  invisible(normalizePath(file, mustWork = TRUE))
}

#' @export
as.data.frame.glyco_panel_audit <- function(x, ...) .audit_table(x)

#' @export
as.data.frame.glyco_panel_set <- function(x, ...) x$metrics

#' @export
as.data.frame.glyco_robustness_audit <- function(x, ...) x$metrics

#' @export
as.data.frame.glyco_pathway_fit <- function(x, ...) {
  data.frame(glycoform = names(x$observed), observed = unname(x$observed),
             predicted = unname(x$predicted), stringsAsFactors = FALSE)
}
