#' Analyse a scalar time series as a stochastic differential equation
#'
#' Estimates binned drift and diffusion (Kramers-Moyal) coefficients from a
#' scalar time series `x`, following the method of Nabeel et al. (2025). This
#' is the main entry point of the package; the result is an object of class
#' `"daddy"` that [dd_drift()], [dd_diffusion()], [dd_fit()],
#' [dd_export_data()], and the `print`/`summary`/`plot` methods all operate
#' on.
#'
#' @details
#' # Deviations from PyDaddy's `Characterize()`
#'
#' * `Dt` defaults to `NULL`, meaning "auto-select from the autocorrelation
#'   time" (`ceiling(autocorr_time / 10)`). PyDaddy's equivalent function,
#'   `pydaddy.Characterize()`, hardcodes a default of `Dt = 1` even though its
#'   own internal auto-selection logic (reachable only through a different,
#'   undocumented code path) computes exactly this quantity. See
#'   `DESIGN_DECISIONS.md`.
#' * All parameters are named explicitly; there is no catch-all `...` that
#'   silently absorbs unrecognized arguments into internal state (as
#'   PyDaddy's constructors do via `self.__dict__.update(kwargs)`). A
#'   misspelled argument name here is a clean error, not a silently-ignored
#'   typo.
#' * No `show_summary` argument -- this is a pure function that does not
#'   print or plot as a side effect. Call `print()`, `summary()`, or `plot()`
#'   on the result explicitly.
#'
#' @param data Numeric vector: the time series to analyse.
#' @param t Either a single positive number (the time interval between
#'   observations) or a numeric vector of timestamps the same length as
#'   `data`.
#' @param Dt Integer lag for the drift estimate, or `NULL` (default) to
#'   auto-select from the autocorrelation time.
#' @param dt Integer lag for the diffusion estimate. Defaults to `Dt` (unlike
#'   PyDaddy, which defaults `dt` to a fixed `1` independent of `Dt`; letting
#'   diffusion and drift share a timescale by default avoids the
#'   index-alignment inconsistencies documented in `DESIGN_DECISIONS.md`
#'   that arise in PyDaddy's own code whenever `Dt != dt`).
#' @param bins Integer number of bins across the data range, or the string
#'   `"auto"` to choose via the Freedman-Diaconis rule. Default `20`.
#' @param inc Optional explicit bin width, used instead of `bins`.
#' @param op_range Optional length-2 vector giving the order-parameter range
#'   to bin over; defaults to the observed data range.
#'
#' @return An object of class `"daddy"`.
#' @export
#' @examples
#' d <- dd_load_sample_data("scalar-pairwise")
#' dd <- dd_analyse(d$x, t = d$t_int)
#' print(dd)
#' dd_drift(dd)
dd_analyse <- function(data, t = 1, Dt = NULL, dt = Dt, bins = 20L, inc = NULL,
                        op_range = NULL) {
  t_int <- dd_validate_series(data, t)
  x <- as.numeric(data)

  autocorr_time <- dd_autocorr_time(x)
  if (is.null(Dt)) Dt <- dd_auto_dt(x)
  if (is.null(dt)) dt <- Dt
  Dt <- as.integer(round(Dt)); dt <- as.integer(round(dt))

  if (identical(bins, "auto")) bins <- dd_auto_bins(x)
  bins <- as.integer(round(bins))

  est <- dd_estimate(x, t_int = t_int, Dt = Dt, dt = dt, bins = bins,
                      inc = inc, op_range = op_range)

  structure(
    list(
      x = x,
      t_int = t_int,
      Dt = Dt,
      dt = dt,
      bins = bins,
      autocorr_time = autocorr_time,
      estimate = est,
      fits = list(drift = NULL, diffusion = NULL),
      vector = FALSE
    ),
    class = "daddy"
  )
}

#' Analyse a bivariate (2D) time series as a vector stochastic differential
#' equation
#'
#' The vector-data analogue of [dd_analyse()]: estimates the drift vector
#' `(F1, F2)` and the diffusion matrix `(G11, G22, G12)` from two co-observed
#' series `x1`, `x2` (e.g. the two coordinates of a moving particle, or two
#' interacting population densities). Returns an object of the same
#' `"daddy"` class as [dd_analyse()] (with `$vector` set to `TRUE`), so
#' [dd_fit()], [dd_drift()], [dd_diffusion()], [dd_simulate()],
#' [dd_export_data()], and the `print`/`summary`/`plot` methods all work on
#' it the same way, just with vector-shaped inputs/outputs.
#'
#' @details
#' # Auto-selecting `Dt`
#' As in [dd_analyse()], `Dt = NULL` (the default) auto-selects a lag from
#' the autocorrelation time. For vector data there is no single series to
#' compute that from, so daddyR follows PyDaddy's own convention here:
#' the autocorrelation time is computed from `x1^2 + x2^2` (e.g. squared
#' distance from the origin, for position data) as a single representative
#' scalar summarizing the joint dynamics.
#'
#' See [dd_estimate_vector()] for the full list of ways the underlying
#' estimator deviates from (and fixes bugs in) PyDaddy's vector code path.
#'
#' @param x1,x2 Numeric vectors of equal length: the two co-observed series.
#' @param t Either a single positive number (the time interval between
#'   observations) or a numeric vector of timestamps the same length as
#'   `x1`/`x2`.
#' @param Dt Integer lag for the drift estimate, or `NULL` (default) to
#'   auto-select from the autocorrelation time of `x1^2 + x2^2`.
#' @param dt Integer lag for the diffusion estimate. Defaults to `Dt`.
#' @param bins Integer number of bins per axis. Default `20`.
#' @param inc_x,inc_y Optional explicit bin widths, used instead of `bins`
#'   on the corresponding axis.
#' @param op_range_x,op_range_y Optional length-2 vectors giving the
#'   order-parameter range to bin over; default to the observed range of
#'   `x1`/`x2`.
#'
#' @return An object of class `"daddy"` with `$vector` set to `TRUE`.
#' @export
#' @examples
#' d <- dd_load_sample_data("vector-pairwise")
#' dd <- dd_analyse_vector(d$x1, d$x2, t = d$t_int)
#' print(dd)
#' dd_drift(dd)
dd_analyse_vector <- function(x1, x2, t = 1, Dt = NULL, dt = Dt, bins = 20L,
                               inc_x = NULL, inc_y = NULL,
                               op_range_x = NULL, op_range_y = NULL) {
  t_int <- dd_validate_series(x1, t)
  if (!is.numeric(x2) || any(is.infinite(x2))) {
    stop("`x2` must be a finite numeric vector.", call. = FALSE)
  }
  if (length(x2) != length(x1)) {
    stop("`x1` and `x2` must be the same length.", call. = FALSE)
  }
  x1 <- as.numeric(x1); x2 <- as.numeric(x2)

  m_sq <- x1^2 + x2^2
  autocorr_time <- dd_autocorr_time(m_sq)
  if (is.null(Dt)) Dt <- dd_auto_dt(m_sq)
  if (is.null(dt)) dt <- Dt
  Dt <- as.integer(round(Dt)); dt <- as.integer(round(dt))

  bins <- as.integer(round(bins))

  est <- dd_estimate_vector(x1, x2, t_int = t_int, Dt = Dt, dt = dt, bins = bins,
                             inc_x = inc_x, inc_y = inc_y,
                             op_range_x = op_range_x, op_range_y = op_range_y)

  structure(
    list(
      x1 = x1, x2 = x2,
      t_int = t_int,
      Dt = Dt,
      dt = dt,
      bins = bins,
      autocorr_time = autocorr_time,
      estimate = est,
      fits = list(F1 = NULL, F2 = NULL, G11 = NULL, G22 = NULL, G12 = NULL),
      vector = TRUE
    ),
    class = "daddy"
  )
}

#' Binned drift estimate
#'
#' @param dd A `"daddy"` object from [dd_analyse()] or [dd_analyse_vector()].
#' @return For a scalar object, a data frame with columns `op`, `avg_drift`,
#'   `drift_se`, `n`. For a vector object, one row per `(x1, x2)` bin, with
#'   columns `x1`, `x2`, `avg_drift1`, `avg_drift2`, `n`.
#' @export
dd_drift <- function(dd) {
  stopifnot(inherits(dd, "daddy"))
  if (isTRUE(dd$vector)) return(dd_drift_vector(dd))
  e <- dd$estimate
  data.frame(op = e$op, avg_drift = e$avg_drift, drift_se = e$drift_ebar, n = e$drift_n)
}

#' @keywords internal
dd_drift_vector <- function(dd) {
  e <- dd$estimate
  grid <- expand.grid(x1 = e$op_x, x2 = e$op_y)
  data.frame(
    x1 = grid$x1, x2 = grid$x2,
    avg_drift1 = as.vector(e$avg_drift1), avg_drift2 = as.vector(e$avg_drift2),
    n = as.vector(e$drift1_n)
  )
}

#' Binned diffusion estimate
#'
#' @param dd A `"daddy"` object from [dd_analyse()] or [dd_analyse_vector()].
#' @return For a scalar object, a data frame with columns `op`, `avg_diff`,
#'   `diff_se`, `n`. For a vector object, one row per `(x1, x2)` bin, with
#'   columns `x1`, `x2`, `avg_diff11`, `avg_diff22`, `avg_diff12`, `n`.
#' @export
dd_diffusion <- function(dd) {
  stopifnot(inherits(dd, "daddy"))
  if (isTRUE(dd$vector)) return(dd_diffusion_vector(dd))
  e <- dd$estimate
  data.frame(op = e$op, avg_diff = e$avg_diff, diff_se = e$diff_ebar, n = e$diff_n)
}

#' @keywords internal
dd_diffusion_vector <- function(dd) {
  e <- dd$estimate
  grid <- expand.grid(x1 = e$op_x, x2 = e$op_y)
  data.frame(
    x1 = grid$x1, x2 = grid$x2,
    avg_diff11 = as.vector(e$avg_diff11), avg_diff22 = as.vector(e$avg_diff22),
    avg_diff12 = as.vector(e$avg_diff12), n = as.vector(e$diff11_n)
  )
}

#' Fit a polynomial to the drift or diffusion function
#'
#' Fits a sparse polynomial to the *unbinned* drift or diffusion series
#' (i.e. to every data point, not to the binned averages), matching
#' PyDaddy's `Daddy.fit()`. Unlike PyDaddy, this does not mutate `dd` in
#' place (R objects are not reference objects); it returns an updated copy,
#' which is the idiomatic R pattern for this kind of "add a computed result"
#' step: `dd <- dd_fit(dd, "drift", degree = 3)`.
#'
#' For a vector (`dd_analyse_vector()`) object, `which` is instead one of
#' `"F1"`, `"F2"` (drift components) or `"G11"`, `"G22"`, `"G12"` (diffusion
#' matrix entries), and the fitted result is a `poly2d` rather than a
#' `poly1d`. There is no separate `"G21"`: the diffusion matrix is symmetric
#' by construction, so `G12` is used for both off-diagonal entries wherever
#' the full matrix is needed (see [dd_simulate()]).
#'
#' @param dd A `"daddy"` object.
#' @param which One of `"drift"`/`"diffusion"` (scalar `dd`), or
#'   `"F1"`/`"F2"`/`"G11"`/`"G22"`/`"G12"` (vector `dd`).
#' @param degree Maximum polynomial degree.
#' @param threshold Sparsity threshold. Ignored if `tune = TRUE`.
#' @param alpha Ridge regularization strength.
#' @param tune If `TRUE`, choose `threshold` automatically via
#'   cross-validation (see [dd_tune_threshold()]/[dd_tune_threshold2d()]).
#' @param ... Additional arguments passed to the tuning function when
#'   `tune = TRUE`.
#' @return The `dd` object, with the corresponding entry of `dd$fits` set to
#'   the fitted `poly1d` (scalar) or `poly2d` (vector).
#' @export
dd_fit <- function(dd, which, degree = 3, threshold = 0, alpha = 0, tune = FALSE, ...) {
  stopifnot(inherits(dd, "daddy"))
  e <- dd$estimate

  if (isTRUE(dd$vector)) {
    which <- match.arg(which, c("F1", "F2", "G11", "G22", "G12"))
    xy <- switch(which,
      F1  = list(x1 = e$x1_for_drift_fit, x2 = e$x2_for_drift_fit, y = e$drift1_series),
      F2  = list(x1 = e$x1_for_drift_fit, x2 = e$x2_for_drift_fit, y = e$drift2_series),
      G11 = list(x1 = e$x1_for_diff_fit, x2 = e$x2_for_diff_fit, y = e$diff11_series),
      G22 = list(x1 = e$x1_for_diff_fit, x2 = e$x2_for_diff_fit, y = e$diff22_series),
      G12 = list(x1 = e$x1_for_diff_fit, x2 = e$x2_for_diff_fit, y = e$diff12_series)
    )
    fit <- if (tune) {
      dd_tune_threshold2d(xy$x1, xy$x2, xy$y, degree = degree, alpha = alpha, ...)
    } else {
      dd_fit_polynomial2d(xy$x1, xy$x2, xy$y, degree = degree, threshold = threshold, alpha = alpha)
    }
    dd$fits[[which]] <- fit
    return(dd)
  }

  which <- match.arg(which, c("drift", "diffusion"))
  if (which == "drift") {
    x <- e$x_for_drift_fit; y <- e$drift_series
  } else {
    x <- e$x_for_diff_fit; y <- e$diff_series
  }

  fit <- if (tune) {
    dd_tune_threshold(x, y, degree = degree, alpha = alpha, ...)
  } else {
    dd_fit_polynomial(x, y, degree = degree, threshold = threshold, alpha = alpha)
  }

  dd$fits[[which]] <- fit
  dd
}

#' Export drift/diffusion estimates as a data frame (and optionally a CSV)
#'
#' @param dd A `"daddy"` object.
#' @param file Optional file path; if given, the data frame is also written
#'   there via [utils::write.csv()].
#' @return A data frame combining `dd_drift()` and `dd_diffusion()` by order
#'   parameter (scalar `dd`), or by `(x1, x2)` bin (vector `dd`).
#' @export
dd_export_data <- function(dd, file = NULL) {
  stopifnot(inherits(dd, "daddy"))
  if (isTRUE(dd$vector)) {
    d1 <- dd_drift_vector(dd); d2 <- dd_diffusion_vector(dd)
    out <- data.frame(
      x1 = d1$x1, x2 = d1$x2,
      avg_drift1 = d1$avg_drift1, avg_drift2 = d1$avg_drift2,
      avg_diff11 = d2$avg_diff11, avg_diff22 = d2$avg_diff22, avg_diff12 = d2$avg_diff12,
      n = d1$n
    )
  } else {
    out <- merge(dd_drift(dd), dd_diffusion(dd), by = "op", suffixes = c("_drift", "_diff"))
    out <- out[order(out$op), ]
  }
  if (!is.null(file)) {
    utils::write.csv(out, file = file, row.names = FALSE)
  }
  out
}
