test_that("full canonical Jacobian contains eight local rate directions", {
  pathway <- canonical_pathway()
  jacobian <- lograte_jacobian(pathway)
  expect_equal(dim(jacobian), c(22, 8))
  expect_equal(qr(jacobian, tol = 1e-8)$rank, 8)
  expect_lt(max(abs(colSums(jacobian))), 1e-8)
})

test_that("composition dimension prevents full rank below nine classes", {
  pathway <- canonical_pathway()
  audit <- check_panel(pathway, pathway$observed_classes[1:8])
  expect_false(audit$dimension_sufficient)
  expect_false(audit$full_rank)
})

test_that("published example nine-class panel is locally full rank", {
  pathway <- canonical_pathway()
  panel <- c("M9", "M8", "M5", "M5Gn", "M3Gn", "M3Gn2", "M3Gn3", "G1", "G2S1")
  audit <- check_panel(pathway, panel, target_sd = 0.25)
  expect_true(audit$dimension_sufficient)
  expect_true(audit$full_rank)
  expect_equal(audit$rank, 8)
  expect_true(is.finite(audit$condition_number))
  expect_true(is.finite(audit$amplification))
  expect_equal(audit$max_assay_cv, audit$target_sd / audit$amplification)
})

test_that("operating-point audit reports every requested point", {
  robust <- check_operating_points(canonical_pathway(),
    c("M9", "M8", "M5", "M5Gn", "M3Gn", "M3Gn2", "M3Gn3", "G1", "G2S1"),
    n = 5, seed = 2)
  expect_equal(nrow(robust$metrics), 5)
  expect_true(robust$full_rank_fraction >= 0 && robust$full_rank_fraction <= 1)
})

test_that("users can audit a declared subset of reaction classes", {
  pathway <- canonical_pathway()
  audit <- check_panel(pathway, c("M9", "M5", "M5Gn"),
                       rate_classes = c("MAN1", "MGAT1"))
  expect_equal(audit$n_rate_classes, 2)
  expect_named(audit$log_rates, pathway$rate_classes)
  expect_equal(colnames(audit$jacobian), c("MAN1", "MGAT1"))
})
