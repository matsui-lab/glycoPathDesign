#' Launch the glycoPathDesign Shiny application
#'
#' @param launch.browser Passed to `shiny::runApp`.
#' @param ... Additional arguments passed to `shiny::runApp`.
#' @export
run_app <- function(launch.browser = interactive(), ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) .gpd_stop("Install shiny to run the application.")
  if (!requireNamespace("DT", quietly = TRUE)) .gpd_stop("Install DT to use the pathway builder in the application.")
  app <- system.file("shiny", package = "glycoPathDesign")
  if (!nzchar(app)) .gpd_stop("The installed Shiny application could not be found.")
  shiny::runApp(app, launch.browser = launch.browser, ...)
}
