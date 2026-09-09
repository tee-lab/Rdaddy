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
#' * No `show_summary` argument — this is a pure function that does not
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
      fits = list(drift = NULL, diffusion = NULL)
    ),
    class = "daddy"
  )
}

#' Binned drift estimate
#' @param dd A `"daddy"` object from [dd_analyse()].
#' @return A data frame with columns `op`, `avg_drift`, `drift_se`, `n`.
#' @export
dd_drift <- function(dd) {
  stopifnot(inherits(dd, "daddy"))
  e <- dd$estimate
  data.frame(op = e$op, avg_drift = e$avg_drift, drift_se = e$drift_ebar, n = e$drift_n)
}

#' Binned diffusion estimate
#' @param dd A `"daddy"` object from [dd_analyse()].
#' @return A data frame with columns `op`, `avg_diff`, `diff_se`, `n`.
#' @export
dd_diffusion <- function(dd) {
  stopifnot(inherits(dd, "daddy"))
  e <- dd$estimate
  data.frame(op = e$op, avg_diff = e$avg_diff, diff_se = e$diff_ebar, n = e$diff_n)
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
#' @param dd A `"daddy"` object.
#' @param which One of `"drift"` or `"diffusion"`.
#' @param degree Maximum polynomial degree.
#' @param threshold Sparsity threshold. Ignored if `tune = TRUE`.
#' @param alpha Ridge regularization strength.
#' @param tune If `TRUE`, choose `threshold` automatically via
#'   cross-validation (see [dd_tune_threshold()]).
#' @param ... Additional arguments passed to [dd_tune_threshold()] when
#'   `tune = TRUE`.
#' @return The `dd` object, with `dd$fits$drift` or `dd$fits$diffusion` set
#'   to the fitted `poly1d`.
#' @export
dd_fit <- function(dd, which = c("drift", "diffusion"), degree = 3,
                    threshold = 0, alpha = 0, tune = FALSE, ...) {
  stopifnot(inherits(dd, "daddy"))
  which <- match.arg(which)
  e <- dd$estimate

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
#'   parameter.
#' @export
dd_export_data <- function(dd, file = NULL) {
  stopifnot(inherits(dd, "daddy"))
  out <- merge(dd_drift(dd), dd_diffusion(dd), by = "op", suffixes = c("_drift", "_diff"))
  out <- out[order(out$op), ]
  if (!is.null(file)) {
    utils::write.csv(out, file = file, row.names = FALSE)
  }
  out
}
