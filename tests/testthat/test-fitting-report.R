test_that("synthetic composition is recovered by the forward model fit", {
  pathway <- canonical_pathway()
  truth <- setNames(c(0.2, -0.1, 0.15, -0.2, 0.1, 0.05, -0.05, 0.12),
                    pathway$rate_classes)
  observed <- predict_composition(pathway, truth)
  fit <- fit_pathway(pathway, observed, start = truth, n_starts = 1)
  expect_s3_class(fit, "glyco_pathway_fit")
  expect_lt(fit$objective, 1e-12)
  expect_gt(fit$fit_correlation, 0.999)
})

test_that("reports can be written in portable formats", {
  audit <- check_panel(canonical_pathway(),
    c("M9", "M8", "M5", "M5Gn", "M3Gn", "M3Gn2", "M3Gn3", "G1", "G2S1"))
  csv <- tempfile(fileext = ".csv")
  html <- tempfile(fileext = ".html")
  expect_true(file.exists(write_design_report(audit, csv)))
  expect_true(file.exists(write_design_report(audit, html)))
  expect_match(paste(readLines(html), collapse = " "), "conditional")
})
