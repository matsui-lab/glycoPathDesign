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

test_that("exact starting solutions converge without hiding real optimiser failures", {
  p <- canonical_pathway()
  for (loss in c("composition", "log_composition")) {
    fit <- fit_pathway(p, predict_composition(p), n_starts = 1, loss = loss)
    expect_equal(fit$objective, 0)
    expect_equal(fit$convergence, 0L)
    expect_equal(unname(fit$rates), rep(1, length(p$rate_classes)))
  }
  truth <- setNames(rep(.7, length(p$rate_classes)), p$rate_classes)
  failed <- fit_pathway(p, predict_composition(p, truth), n_starts = 1,
                        control = list(maxit = 0))
  expect_gt(failed$objective, 0)
  expect_true(failed$convergence != 0L)
})


test_that("CSV round-off at an exact starting solution does not cause false failure", {
  p <- canonical_pathway()
  f <- tempfile(fileext = ".csv")
  y <- predict_composition(p)
  write.csv(data.frame(glycoform = names(y), proportion = as.numeric(y)), f, row.names = FALSE)
  for (loss in c("composition", "log_composition")) {
    fit <- fit_pathway(p, read.csv(f), n_starts = 1, loss = loss)
    expect_equal(fit$convergence, 0L)
    expect_lt(fit$objective, 1e-25)
    expect_equal(unname(fit$rates), rep(1, length(p$rate_classes)))
  }
})
