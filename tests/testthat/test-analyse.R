test_that("dd_analyse builds a well-formed daddy object with sensible defaults", {
  d <- dd_load_sample_data("scalar-pairwise")
  dd <- dd_analyse(d$x, t = d$t)

  expect_s3_class(dd, "daddy")
  expect_equal(dd$dt, dd$Dt)  # dt defaults to Dt
  expect_true(dd$Dt >= 1L)
  expect_equal(dd$bins, 20L)
  expect_equal(nrow(dd_drift(dd)), 20L)
  expect_equal(nrow(dd_diffusion(dd)), 20L)
})

test_that("dd_analyse rejects unknown/misspelled arguments instead of silently ignoring them", {
  d <- dd_load_sample_data("scalar-pairwise")
  expect_error(dd_analyse(d$x, t = d$t, binz = 20), "unused argument")
})

test_that("bins = 'auto' selects a bin count via Freedman-Diaconis", {
  d <- dd_load_sample_data("scalar-pairwise")
  dd <- dd_analyse(d$x, t = d$t, bins = "auto")
  expect_true(dd$bins >= 1L)
  expect_equal(nrow(dd_drift(dd)), dd$bins)
})

test_that("dd_fit attaches a poly1d fit without mutating the caller's copy", {
  d <- dd_load_sample_data("scalar-pairwise")
  dd <- dd_analyse(d$x, t = d$t)
  expect_null(dd$fits$drift)

  dd_fitted <- dd_fit(dd, "drift", degree = 3, threshold = 0)
  expect_s3_class(dd_fitted$fits$drift, "poly1d")
  expect_null(dd$fits$drift)  # original object untouched (R copy-on-modify)
})

test_that("dd_export_data returns a merged, sorted data frame", {
  d <- dd_load_sample_data("scalar-pairwise")
  dd <- dd_analyse(d$x, t = d$t)
  out <- dd_export_data(dd)
  expect_true(all(c("op", "avg_drift", "avg_diff", "n_drift", "n_diff") %in% names(out)))
  expect_equal(out$op, sort(out$op))
})

test_that("print/summary/plot methods run without error", {
  d <- dd_load_sample_data("scalar-pairwise")
  dd <- dd_analyse(d$x, t = d$t)
  dd <- dd_fit(dd, "drift", degree = 3, threshold = 0)
  dd <- dd_fit(dd, "diffusion", degree = 2, threshold = 0)

  expect_output(print(dd), "daddy")
  expect_output(summary(dd), "Summary of daddy")

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  p <- plot(dd)
  expect_s3_class(p, "patchwork")
})
