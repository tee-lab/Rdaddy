test_that("dd_matrix_sqrt_psd reproduces a positive-definite matrix exactly", {
  G <- matrix(c(2, 0.5, 0.5, 1.5), nrow = 2)
  S <- dd_matrix_sqrt_psd(G)
  expect_equal(S %*% S, G, tolerance = 1e-8)
  expect_equal(S, t(S), tolerance = 1e-8)  # symmetric square root
})

test_that("dd_matrix_sqrt_psd clips negative eigenvalues with a warning, rather than erroring", {
  # An indefinite "diffusion matrix" (can happen if a fitted polynomial goes
  # negative off-support) -- PyDaddy's scipy.linalg.sqrtm() is not guaranteed
  # to raise here; daddyR always returns a valid real matrix instead.
  G_bad <- matrix(c(-1, 0, 0, 1), nrow = 2)
  expect_warning(S <- dd_matrix_sqrt_psd(G_bad), "negative eigenvalue")
  expect_true(all(is.finite(S)))
  # The clipped square root's square should be PSD and match on the
  # non-negative eigenspace (the y-direction here).
  expect_equal((S %*% S)[2, 2], 1, tolerance = 1e-8)
})

test_that("dd_simulate (scalar) produces a finite trajectory of the right length", {
  set.seed(10)
  drift <- new_poly1d(c(0, -0.5))       # f(x) = -0.5*x
  diffusion <- new_poly1d(c(0.04))      # g^2(x) = 0.04 (constant)
  dd <- structure(list(fits = list(drift = drift, diffusion = diffusion), vector = FALSE),
                   class = "daddy")
  x <- dd_simulate(dd, t_int = 0.01, timepoints = 500, x0 = 0.2)
  expect_length(x, 500)
  expect_true(all(is.finite(x)))
  expect_equal(x[1], 0.2)
})

test_that("dd_simulate errors clearly when required fits are missing", {
  dd_scalar <- structure(list(fits = list(drift = NULL, diffusion = NULL), vector = FALSE),
                          class = "daddy")
  expect_error(dd_simulate(dd_scalar, t_int = 0.01, timepoints = 10), "Fit both drift")

  dd_vec <- structure(list(fits = list(F1 = NULL, F2 = NULL, G11 = NULL, G22 = NULL, G12 = NULL),
                            vector = TRUE),
                       class = "daddy")
  expect_error(dd_simulate(dd_vec, t_int = 0.01, timepoints = 10), "Fit F1, F2, G11")
})

test_that("dd_simulate (vector) produces a finite trajectory with independent noise when G12 is unset", {
  set.seed(11)
  F1 <- new_poly2d(c(0, -0.5, 0), degree = 1)   # F1(x1,x2) = -0.5*x1
  F2 <- new_poly2d(c(0, 0, -0.8), degree = 1)   # F2(x1,x2) = -0.8*x2
  G11 <- new_poly2d(c(0.04, 0, 0), degree = 1)
  G22 <- new_poly2d(c(0.09, 0, 0), degree = 1)
  dd <- structure(
    list(fits = list(F1 = F1, F2 = F2, G11 = G11, G22 = G22, G12 = NULL), vector = TRUE),
    class = "daddy"
  )
  x <- dd_simulate(dd, t_int = 0.01, timepoints = 500, x0 = c(0.1, -0.1))
  expect_equal(dim(x), c(500L, 2L))
  expect_true(all(is.finite(x)))
  expect_equal(x[1, ], c(x1 = 0.1, x2 = -0.1))
})

test_that("dd_simulate (vector) with a cross-diffusion term runs and stays finite", {
  set.seed(12)
  F1 <- new_poly2d(c(0, -0.5, 0), degree = 1)
  F2 <- new_poly2d(c(0, 0, -0.8), degree = 1)
  G11 <- new_poly2d(c(0.04, 0, 0), degree = 1)
  G22 <- new_poly2d(c(0.09, 0, 0), degree = 1)
  G12 <- new_poly2d(c(0.02, 0, 0), degree = 1)
  dd <- structure(
    list(fits = list(F1 = F1, F2 = F2, G11 = G11, G22 = G22, G12 = G12), vector = TRUE),
    class = "daddy"
  )
  x <- dd_simulate(dd, t_int = 0.01, timepoints = 500, x0 = c(0, 0))
  expect_equal(dim(x), c(500L, 2L))
  expect_true(all(is.finite(x)))
})
