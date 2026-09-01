`%||%` <- function(x, y) if (is.null(x)) y else x

.gpd_stop <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

.as_data_frame <- function(x, name) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    x <- utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(x)) .gpd_stop("%s must be a data.frame or a CSV path.", name)
  x
}

.require_columns <- function(x, columns, name) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    .gpd_stop("%s is missing required column(s): %s.", name, paste(missing, collapse = ", "))
  }
  invisible(x)
}

.normalise <- function(x, name = "values") {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < 0)) .gpd_stop("%s must be finite and non-negative.", name)
  total <- sum(x)
  if (!is.finite(total) || total <= 0) .gpd_stop("%s must have a positive sum.", name)
  x / total
}

.named_numeric <- function(x, expected, name, default = NULL) {
  if (is.null(x)) x <- default
  if (is.null(x)) .gpd_stop("%s is required.", name)
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- original_names
  if (is.null(names(x)) && length(x) == length(expected)) names(x) <- expected
  if (is.null(names(x))) .gpd_stop("%s must be named or have length %d.", name, length(expected))
  missing <- setdiff(expected, names(x))
  extra <- setdiff(names(x), expected)
  if (length(missing)) .gpd_stop("%s is missing: %s.", name, paste(missing, collapse = ", "))
  if (length(extra)) x <- x[setdiff(names(x), extra)]
  x <- x[expected]
  if (any(!is.finite(x))) .gpd_stop("%s must contain only finite values.", name)
  x
}

.matrix_rank <- function(x, tol = 1e-9) {
  if (!length(x) || !nrow(x) || !ncol(x)) return(0L)
  values <- svd(x, nu = 0, nv = 0)$d
  as.integer(sum(values > tol))
}

.safe_correlation <- function(x, y) {
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

.format_finite <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg"), "undefined")
}
