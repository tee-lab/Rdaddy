#' Load a bundled example dataset
#'
#' Loads one of the small time series datasets shipped with daddyR for
#' examples and testing, mirroring PyDaddy's `load_sample_dataset()`.
#' `"scalar-pairwise"`/`"scalar-ternary"` are scalar (1D) series, ready for
#' [dd_analyse()]; `"vector-pairwise"`/`"vector-ternary"` are bivariate (2D)
#' series, ready for [dd_analyse_vector()].
#'
#' @details
#' The vector datasets are a contiguous 20,000-point subsample (rounded to
#' 6 decimal places) of PyDaddy's own bundled
#' `model-data-vector-{pairwise,ternary}` datasets, which are simulated
#' trajectories from 2D pairwise- and ternary-interaction collective-motion
#' models (Jhawar & Guttal 2020) at a fixed observation interval of
#' `t_int = 0.12` -- subsampled (from an original ~600,000 points) purely to
#' keep package size reasonable; see `data-raw/prepare_data.R` for the exact
#' subsampling code.
#'
#' @param name One of `"scalar-pairwise"`, `"scalar-ternary"`,
#'   `"vector-pairwise"`, or `"vector-ternary"`.
#' @return For a scalar dataset, a list with `x` (the time series) and
#'   `t_int` (the time interval between observations), ready to pass to
#'   [dd_analyse()] as `dd_analyse(out$x, t = out$t_int)`. For a vector
#'   dataset, a list with `x1`, `x2` (the two co-observed series) and
#'   `t_int`, ready to pass to [dd_analyse_vector()] as
#'   `dd_analyse_vector(out$x1, out$x2, t = out$t_int)`.
#' @export
#' @examples
#' d <- dd_load_sample_data("scalar-pairwise")
#' str(d)
#' dv <- dd_load_sample_data("vector-pairwise")
#' str(dv)
dd_load_sample_data <- function(name = c("scalar-pairwise", "scalar-ternary",
                                          "vector-pairwise", "vector-ternary")) {
  name <- match.arg(name)
  if (startsWith(name, "vector")) {
    file <- switch(name,
      "vector-pairwise" = "vector_pairwise.csv",
      "vector-ternary" = "vector_ternary.csv"
    )
    path <- system.file("extdata", file, package = "daddyR")
    if (!nzchar(path)) stop("Could not find bundled dataset '", file, "'.", call. = FALSE)
    raw <- utils::read.csv(path, header = FALSE)
    return(list(x1 = raw[[1]], x2 = raw[[2]], t_int = 0.12))
  }
  file <- switch(name,
    "scalar-pairwise" = "scalar_pairwise.csv",
    "scalar-ternary" = "scalar_ternary.csv"
  )
  path <- system.file("extdata", file, package = "daddyR")
  if (!nzchar(path)) stop("Could not find bundled dataset '", file, "'.", call. = FALSE)
  raw <- utils::read.csv(path, header = FALSE)
  x <- raw[[1]]
  t <- raw[[2]]
  t_int <- (t[length(t)] - t[1]) / (length(t) - 1)
  list(x = x, t = t, t_int = t_int)
}
