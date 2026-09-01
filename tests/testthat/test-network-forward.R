test_that("canonical pathway is valid and conserves observed composition", {
  pathway <- canonical_pathway()
  expect_s3_class(pathway, "glyco_pathway")
  expect_equal(nrow(pathway$nodes), 22)
  expect_equal(nrow(pathway$edges), 25)
  expect_equal(length(pathway$rate_classes), 8)
  predicted <- predict_composition(pathway)
  expect_equal(sum(predicted), 1, tolerance = 1e-12)
  expect_true(all(predicted >= 0))
})

test_that("observation mapping aggregates latent states", {
  pathway <- glyco_pathway(
    data.frame(id = c("A", "B", "C")),
    data.frame(from = c("A", "A"), to = c("B", "C"),
               rate_class = c("left", "right")),
    observations = data.frame(node = c("A", "B", "C"),
                              glycoform = c("precursor", "product", "product")),
    entry = "A"
  )
  predicted <- predict_composition(pathway)
  expect_named(predicted, c("precursor", "product"))
  expect_equal(sum(predicted), 1)
})

test_that("invalid graphs fail clearly", {
  expect_error(glyco_pathway(data.frame(id = c("A", "A")),
                             data.frame(from = "A", to = "A", rate_class = "r")),
               "unique")
  expect_error(glyco_pathway(data.frame(id = c("A", "B")),
                             data.frame(from = "A", to = "X", rate_class = "r")),
               "unknown")
})

test_that("pathway directories round-trip without losing boundary conditions", {
  original <- canonical_pathway()
  directory <- tempfile("pathway-")
  write_pathway(original, directory)
  restored <- read_pathway_directory(directory)
  expect_equal(restored$nodes, original$nodes)
  expect_equal(restored$edges, original$edges)
  expect_equal(restored$entry, original$entry)
  expect_equal(restored$secretion, original$secretion)
})
