test_that("panel enumeration distinguishes exhaustive and sampled searches", {
  pathway <- canonical_pathway()
  candidates <- pathway$observed_classes[1:10]
  exhaustive <- enumerate_panels(pathway, 9, candidates = candidates, max_panels = 100)
  expect_false(exhaustive$sampled)
  expect_equal(exhaustive$total_combinations, 10)
  expect_equal(nrow(exhaustive$metrics), 10)

  sampled <- enumerate_panels(pathway, 5, candidates = candidates, max_panels = 7, seed = 3)
  expect_true(sampled$sampled)
  expect_equal(sampled$evaluated, 7)
  expect_equal(length(unique(sampled$metrics$panel)), 7)
})

test_that("optimisation returns only full-rank candidates when requested", {
  pathway <- canonical_pathway()
  candidates <- c("M9", "M8", "M5", "M5Gn", "M3Gn", "M3Gn2", "M3Gn3", "G1", "G2S1", "G2")
  panels <- enumerate_panels(pathway, 9, candidates = candidates, max_panels = 100)
  best <- optimise_panels(panels, n = 20)
  expect_true(all(best$full_rank))
  if (nrow(best) > 1) expect_true(all(diff(best$amplification) >= 0))
})
