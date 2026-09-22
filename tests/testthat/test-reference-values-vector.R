# These tests check daddyR's vector estimator against reference numbers
# produced by actually running the original Python PyDaddy package (v1.1.1,
# patched only for its numpy>=1.24 incompatibility -- see
# DESIGN_DECISIONS.md) on the bundled vector_pairwise dataset, with
# bins = 20, Dt = dt = 1, t_int = 0.12 (matching PyDaddy's own default for
# this dataset). See tests/testthat/fixtures/ for the generated CSVs, and
# `data-raw/prepare_data.R` for how the bundled dataset itself was derived.
#
# As with the scalar reference tests, two kinds of quantity are compared
# differently:
#  - The raw (unbinned) drift series (F1/F2 inputs) do not depend on binning
#    at all and use an identical formula to PyDaddy's, so they -- and the
#    resulting fitted F1/F2 polynomials (threshold = 0, an exact OLS fit on
#    identical (x, y) pairs) -- should match to numerical precision.
#  - The diffusion matrix (G11/G22/G12) uses a deliberately different
#    residual formula from PyDaddy's (see DESIGN_DECISIONS.md and
#    dd_estimate_vector()'s documentation) which happens to coincide with
#    PyDaddy's own formula when Dt == 1 (as here) *except* for daddyR's
#    corrected bin edges feeding into the drift-subtraction step -- so these
#    are compared only loosely, as a sanity check that the algorithm is
#    computing something in the right ballpark, not a precise match.

fixture_path <- function(name) testthat::test_path("fixtures", name)

test_that("raw drift series matches PyDaddy exactly (bin-independent quantity)", {
  d <- dd_load_sample_data("vector-pairwise")
  params <- utils::read.csv(fixture_path("vec_params.csv"), header = FALSE)
  expected_t_int <- params[params[[1]] == "t_int", 2]
  expect_equal(d$t_int, expected_t_int, tolerance = 1e-8)

  est <- dd_estimate_vector(d$x1, d$x2, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  head_fixture <- utils::read.csv(fixture_path("vec_drift_head.csv"))
  n <- nrow(head_fixture)

  expect_equal(est$x1_for_drift_fit[1:n], head_fixture$x1, tolerance = 1e-6)
  expect_equal(est$x2_for_drift_fit[1:n], head_fixture$x2, tolerance = 1e-6)
  expect_equal(est$drift1_series[1:n], head_fixture$driftX, tolerance = 1e-4)
  expect_equal(est$drift2_series[1:n], head_fixture$driftY, tolerance = 1e-4)
})

test_that("F1/F2 polynomial fits (threshold = 0) match PyDaddy's ridge fit", {
  # threshold = 0 makes STLSQ a single plain OLS fit; since the (x, y) pairs
  # being fit are identical to PyDaddy's for F1/F2 (same drift formula, same
  # data), the fitted *function* should agree, independent of the two
  # implementations' different term orderings -- so this compares predicted
  # values on a grid, not raw coefficient vectors.
  d <- dd_load_sample_data("vector-pairwise")
  est <- dd_estimate_vector(d$x1, d$x2, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  F1 <- dd_fit_polynomial2d(est$x1_for_drift_fit, est$x2_for_drift_fit, est$drift1_series,
                             degree = 3, threshold = 0)
  F2 <- dd_fit_polynomial2d(est$x1_for_drift_fit, est$x2_for_drift_fit, est$drift2_series,
                             degree = 3, threshold = 0)

  grid <- utils::read.csv(fixture_path("vec_poly_grid_pairwise.csv"))
  F1_pred <- stats::predict(F1, grid$x1, grid$x2)
  F2_pred <- stats::predict(F2, grid$x1, grid$x2)

  expect_equal(F1_pred, grid$F1, tolerance = 1e-3)
  expect_equal(F2_pred, grid$F2, tolerance = 1e-3)
})

test_that("G11/G22/G12 polynomial fits are in the same ballpark as PyDaddy's", {
  # See file header: daddyR's diffusion-matrix formula deliberately differs
  # from PyDaddy's (fixing the Dt/dt conflation described in
  # dd_estimate_vector()'s docs), so exact agreement isn't expected here --
  # only that the fitted surfaces are positively correlated with PyDaddy's
  # and of a similar order of magnitude.
  d <- dd_load_sample_data("vector-pairwise")
  est <- dd_estimate_vector(d$x1, d$x2, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  G11 <- dd_fit_polynomial2d(est$x1_for_diff_fit, est$x2_for_diff_fit, est$diff11_series,
                              degree = 2, threshold = 0)
  G22 <- dd_fit_polynomial2d(est$x1_for_diff_fit, est$x2_for_diff_fit, est$diff22_series,
                              degree = 2, threshold = 0)
  G12 <- dd_fit_polynomial2d(est$x1_for_diff_fit, est$x2_for_diff_fit, est$diff12_series,
                              degree = 2, threshold = 0)

  grid <- utils::read.csv(fixture_path("vec_poly_grid_pairwise.csv"))
  G11_pred <- stats::predict(G11, grid$x1, grid$x2)
  G22_pred <- stats::predict(G22, grid$x1, grid$x2)
  G12_pred <- stats::predict(G12, grid$x1, grid$x2)

  expect_gt(stats::cor(G11_pred, grid$G11), 0.5)
  expect_gt(stats::cor(G22_pred, grid$G22), 0.5)
  # G12 is small in magnitude and noisier by nature (a cross term near zero
  # for this dataset); just check it's the same order of magnitude, not
  # tightly correlated.
  expect_lt(stats::sd(G12_pred), 10 * stats::sd(grid$G12) + 1e-6)
})
