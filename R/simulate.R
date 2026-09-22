#' Symmetric positive-semidefinite matrix square root (via eigendecomposition)
#'
#' Used by [dd_simulate()]'s vector path to turn a fitted diffusion matrix
#' `G = [[G11, G12], [G12, G22]]` into a noise-coefficient matrix `S` with
#' `S %*% S == G`, for correlated 2D Wiener noise.
#'
#' @details
#' # Deviation from PyDaddy
#' PyDaddy's `simulate()` calls `scipy.linalg.sqrtm(G)` directly and only
#' reacts if that call happens to raise `LinAlgError`/`ValueError` -- but
#' `sqrtm` is not guaranteed to raise on an indefinite matrix; it can
#' silently return a matrix with a small (or not-so-small) imaginary part,
#' which then propagates into the simulated trajectory as `NaN`/complex
#' values with no clear error message. This can happen quite easily in
#' practice, because the fitted `G11`/`G22`/`G12` are unconstrained
#' polynomials that can evaluate to values outside a valid diffusion matrix
#' (e.g. negative "variance") once a simulated trajectory wanders outside
#' the range of the original training data.
#'
#' daddyR instead computes the symmetric matrix square root directly via
#' eigendecomposition and explicitly clips any negative eigenvalues to zero
#' before taking their square root -- so the result is always a valid,
#' finite, real matrix, at the cost of silently treating a locally-invalid
#' (indefinite) fitted diffusion matrix as the nearest valid one, rather than
#' erroring. A warning is issued when clipping actually changes the matrix,
#' so this doesn't happen silently.
#'
#' @param G A symmetric numeric matrix (need not be exactly symmetric; it is
#'   symmetrized first as `(G + t(G))/2`).
#' @return A symmetric matrix `S` such that `S %*% S` approximates `G` (exactly
#'   equal if `G` was already positive-semidefinite).
#' @keywords internal
dd_matrix_sqrt_psd <- function(G) {
  G <- (G + t(G)) / 2
  eig <- eigen(G, symmetric = TRUE)
  if (any(eig$values < -1e-8)) {
    warning("Diffusion matrix has a negative eigenvalue (", signif(min(eig$values), 3),
            "); clipping to 0. This usually means the simulated trajectory has ",
            "wandered outside the range the drift/diffusion polynomials were fit on.",
            call. = FALSE)
  }
  vals <- pmax(eig$values, 0)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(eig$vectors)
}

#' Simulate a scalar SDE trajectory from fitted drift/diffusion polynomials
#'
#' Uses the Euler-Maruyama scheme: `x[i+1] = x[i] + f(x[i])*h + g(x[i])*dW`,
#' with `dW ~ N(0, h)` and `g(x) = sqrt(max(diffusion(x), 0))`.
#'
#' @param drift,diffusion `poly1d` objects (e.g. `dd$fits$drift`,
#'   `dd$fits$diffusion`) giving `f(x)` and `g^2(x)`.
#' @param t_int Numeric, simulation time step.
#' @param timepoints Integer, number of time points to simulate (including
#'   the initial condition).
#' @param x0 Numeric, initial condition.
#' @return Numeric vector of length `timepoints`.
#' @keywords internal
dd_simulate_scalar <- function(drift, diffusion, t_int, timepoints, x0 = 0) {
  x <- numeric(timepoints)
  x[1] <- x0
  dW <- stats::rnorm(timepoints - 1) * sqrt(t_int)
  for (i in seq_len(timepoints - 1)) {
    f <- stats::predict(drift, x[i])
    g <- sqrt(max(stats::predict(diffusion, x[i]), 0))
    x[i + 1] <- x[i] + f * t_int + g * dW[i]
  }
  x
}

#' Simulate a vector (2D) SDE trajectory from fitted drift/diffusion polynomials
#'
#' Uses the Euler-Maruyama scheme, with the two Wiener increments correlated
#' through the diffusion matrix: `x[i+1,] = x[i,] + f(x[i,])*h + S(x[i,]) %*%
#' dW`, where `dW ~ N(0, h*I)` and `S` is the symmetric square root
#' ([dd_matrix_sqrt_psd()]) of `G = [[G11, G12], [G12, G22]]`.
#'
#' @details
#' # Deviation from PyDaddy
#' PyDaddy's vector `simulate()` uses plain Euler-Maruyama when a
#' cross-diffusion term is present but a diagonal-only shortcut otherwise
#' (`sqrt(abs(G11))`, `sqrt(abs(G22))`, independent noise) -- and, separately,
#' uses a *different, higher-order* integrator (`sdeint.itoSRI2`, a
#' stochastic Runge-Kutta scheme) for its scalar `simulate()`. daddyR uses
#' one integrator (Euler-Maruyama) and one code path (always build the full
#' 2x2 diffusion matrix and take its matrix square root, via
#' [dd_matrix_sqrt_psd()]) for both scalar and vector cases, rather than
#' reproducing this asymmetry. Euler-Maruyama has lower strong convergence
#' order than a Runge-Kutta scheme, so for high-precision work, use a
#' smaller `t_int` (finer time step) than you might otherwise need.
#'
#' @param F1,F2 `poly2d` objects giving the drift components.
#' @param G11,G22 `poly2d` objects giving the diagonal diffusion terms.
#' @param G12 A `poly2d` object giving the cross-diffusion term, or `NULL`
#'   (equivalent to a `poly2d` that is identically zero: independent noise).
#' @param t_int Numeric, simulation time step.
#' @param timepoints Integer, number of time points to simulate (including
#'   the initial condition).
#' @param x0 Numeric length-2 vector, initial condition `c(x1, x2)`.
#' @return A numeric matrix with `timepoints` rows and 2 columns (`x1`, `x2`).
#' @keywords internal
dd_simulate_vector <- function(F1, F2, G11, G22, G12, t_int, timepoints, x0 = c(0, 0)) {
  x <- matrix(NA_real_, timepoints, 2)
  colnames(x) <- c("x1", "x2")
  x[1, ] <- x0
  dW <- matrix(stats::rnorm(2 * (timepoints - 1)), ncol = 2) * sqrt(t_int)
  for (i in seq_len(timepoints - 1)) {
    xi1 <- x[i, 1]; xi2 <- x[i, 2]
    f <- c(stats::predict(F1, xi1, xi2), stats::predict(F2, xi1, xi2))
    g11 <- stats::predict(G11, xi1, xi2)
    g22 <- stats::predict(G22, xi1, xi2)
    g12 <- if (!is.null(G12)) stats::predict(G12, xi1, xi2) else 0
    G <- matrix(c(g11, g12, g12, g22), nrow = 2)
    S <- dd_matrix_sqrt_psd(G)
    x[i + 1, ] <- c(xi1, xi2) + f * t_int + as.vector(S %*% dW[i, ])
  }
  x
}

#' Simulate a trajectory from a fitted daddy model
#'
#' Generates a synthetic trajectory by numerically integrating the fitted
#' drift/diffusion (scalar) or `F1`/`F2`/`G11`/`G22`/`G12` (vector)
#' polynomials forward in time via Euler-Maruyama. Useful both as a
#' generative model in its own right and, compared back against the original
#' data (e.g. by re-running [dd_analyse()] on the simulated series), as a
#' self-consistency check on the fit -- PyDaddy's `model_diagnostics()` does
#' the same round-trip.
#'
#' @param dd A `"daddy"` object with the required fit(s) already set (see
#'   [dd_fit()]): `drift`/`diffusion` for a scalar object, or
#'   `F1`/`F2`/`G11`/`G22` (and, optionally, `G12`) for a vector object.
#' @param t_int Numeric, simulation time step (independent of the `t_int`
#'   the original data was observed at).
#' @param timepoints Integer, number of time points to simulate.
#' @param x0 Initial condition: a single number for a scalar object
#'   (default `0`), or a length-2 vector `c(x1, x2)` for a vector object
#'   (default `c(0, 0)`).
#' @return A numeric vector of length `timepoints` (scalar case), or a
#'   `timepoints` by 2 matrix with columns `x1`, `x2` (vector case).
#' @export
#' @examples
#' d <- dd_load_sample_data("scalar-pairwise")
#' dd <- dd_analyse(d$x, t = d$t_int)
#' dd <- dd_fit(dd, "drift", degree = 3, threshold = 0.01)
#' dd <- dd_fit(dd, "diffusion", degree = 2, threshold = 0.01)
#' sim <- dd_simulate(dd, t_int = d$t_int, timepoints = 500)
dd_simulate <- function(dd, t_int, timepoints, x0 = NULL) {
  stopifnot(inherits(dd, "daddy"))
  if (isTRUE(dd$vector)) {
    if (is.null(dd$fits$F1) || is.null(dd$fits$F2) ||
        is.null(dd$fits$G11) || is.null(dd$fits$G22)) {
      stop("Fit F1, F2, G11, and G22 (see dd_fit()) before calling dd_simulate(). ",
           "G12 (cross-diffusion) is optional; if unfit, independent noise is used.",
           call. = FALSE)
    }
    if (is.null(x0)) x0 <- c(0, 0)
    dd_simulate_vector(dd$fits$F1, dd$fits$F2, dd$fits$G11, dd$fits$G22, dd$fits$G12,
                        t_int = t_int, timepoints = timepoints, x0 = x0)
  } else {
    if (is.null(dd$fits$drift) || is.null(dd$fits$diffusion)) {
      stop("Fit both drift and diffusion (see dd_fit()) before calling dd_simulate().",
           call. = FALSE)
    }
    if (is.null(x0)) x0 <- 0
    dd_simulate_scalar(dd$fits$drift, dd$fits$diffusion,
                        t_int = t_int, timepoints = timepoints, x0 = x0)
  }
}
