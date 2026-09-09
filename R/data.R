#' Load a bundled example dataset
#'
#' Loads one of the small time series datasets shipped with daddyR for
#' examples and testing, mirroring PyDaddy's `load_sample_dataset()`.
#' Currently only scalar (1D) datasets are bundled; vector datasets (e.g.
#' the fish-schooling and cell-migration data from the original paper) will
#' be added once vector-data support is ported (see the package roadmap).
#'
#' @param name One of `"scalar-pairwise"` or `"scalar-ternary"`: simulated
#'   time series from simple pairwise- and ternary-interaction models
#'   (Jhawar & Guttal 2020), used throughout the PyDaddy paper's synthetic
#'   examples.
#' @return A list with `x` (the time series) and `t_int` (the time interval
#'   between observations), ready to pass to [dd_analyse()] as
#'   `dd_analyse(out$x, t = out$t_int)`.
#' @export
#' @examples
#' d <- dd_load_sample_data("scalar-pairwise")
#' str(d)
dd_load_sample_data <- function(name = c("scalar-pairwise", "scalar-ternary")) {
  name <- match.arg(name)
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
