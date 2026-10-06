# Install glycoPathDesign and its Shiny dependencies before starting this entry point.
if (!requireNamespace("glycoPathDesign", quietly = TRUE)) {
  stop("Install the glycoPathDesign source package before launching this app.")
}
for (package in c("shiny", "DT")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(paste("Install", package, "before launching this app."))
}
app_dir <- system.file("shiny", package = "glycoPathDesign", mustWork = TRUE)
shiny::shinyAppDir(app_dir)
