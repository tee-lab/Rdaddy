#' Ridge regression (closed form)
#'
#' Equivalent to scikit-learn's `sklearn.linear_model.ridge_regression` with
#' default settings: no intercept centering (an intercept column must be
#' included in `X` explicitly, as the polynomial dictionaries here do), and
#' an optional `alpha` L2 penalty applied uniformly to every column,
#' including that intercept column, exactly as PyDaddy relies on.
#'
#' @param X Numeric design matrix (n x p).
#' @param y Numeric response vector (length n).
#' @param alpha Numeric >= 0, ridge penalty.
#' @param weights Optional numeric vector of sample weights (length n).
#' @return Numeric coefficient vector (length p).
#' @keywords internal
dd_ridge_fit <- function(X, y, alpha = 0, weights = NULL) {
  if (!is.null(weights)) {
    w <- sqrt(weights)
    X <- X * w
    y <- y * w
  }
  p <- ncol(X)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  if (alpha > 0) {
    XtX <- XtX + diag(alpha, p)
  }
  as.vector(solve(XtX, Xty))
}

#' Sequentially-thresholded ridge regression (STLSQ)
#'
#' Fits `y` as a sparse linear combination of the columns of `dictionary`:
#' ridge-fit all currently-active terms, zero out any coefficient with
#' magnitude at or below `threshold`, and repeat until the active set
#' stabilizes (or every term is eliminated). With `threshold = 0` this is
#' just ordinary (or ridge, if `alpha > 0`) regression.
#'
#' @details
#' # Deviation from PyDaddy
#' PyDaddy computes coefficient standard errors by inverting
#' `t(dictionary) %*% dictionary` for the *full* dictionary matrix, even
#' though most of its columns may have been thresholded to zero. Polynomial
#' dictionaries are highly collinear, so this matrix is frequently
#' ill-conditioned or exactly singular once the degree is more than a few,
#' occasionally raising a linear algebra error or silently returning inflated
#' standard errors for terms that were dropped. daddyR instead inverts only
#' the submatrix of *active* (kept) columns, which is what the standard error
#' formula is actually meant to describe, and falls back to `NA` standard
#' errors (with a warning) if that submatrix is itself singular or
#' under-determined, rather than erroring out.
#'
#' @param dictionary Numeric design matrix (n x p).
#' @param y Numeric response vector.
#' @param threshold Numeric >= 0, sparsity threshold.
#' @param alpha Numeric >= 0, ridge penalty.
#' @param weights Optional sample weights.
#' @return A list with `coeffs` (length p, zeroed for dropped terms) and
#'   `stderr` (length p, `NA` for dropped terms or if standard errors could
#'   not be computed).
#' @keywords internal
dd_stlsq <- function(dictionary, y, threshold = 0, alpha = 0, weights = NULL) {
  p <- ncol(dictionary)
  n <- nrow(dictionary)
  coeffs <- numeric(p)
  keep <- rep(TRUE, p)

  for (iter in seq_len(p)) {
    if (!any(keep)) {
      warning("Sparsity threshold eliminated all terms.", call. = FALSE)
      break
    }
    sub <- dd_ridge_fit(dictionary[, keep, drop = FALSE], y, alpha = alpha, weights = weights)
    coeffs[keep] <- sub
    new_keep <- abs(coeffs) > threshold
    coeffs[!new_keep] <- 0
    if (identical(new_keep, keep)) { keep <- new_keep; break }
    keep <- new_keep
  }

  stderr <- rep(NA_real_, p)
  n_active <- sum(keep)
  if (n_active > 0 && n > n_active) {
    Xk <- dictionary[, keep, drop = FALSE]
    yhat <- dictionary %*% coeffs
    rss <- sum((y - yhat)^2)
    sigma2 <- rss / (n - n_active)
    se_active <- tryCatch({
      XtX_inv <- solve(crossprod(Xk))
      sqrt(pmax(diag(XtX_inv) * sigma2, 0))
    }, error = function(e) {
      warning("Could not compute standard errors: the active-term design ",
              "matrix is singular (try a lower polynomial degree, or a ",
              "smaller threshold to keep more collinear terms together).",
              call. = FALSE)
      rep(NA_real_, n_active)
    })
    stderr[keep] <- se_active
  }

  list(coeffs = coeffs, stderr = stderr)
}

#' Build a polynomial design matrix `[1, x, x^2, ..., x^degree]`
#' @keywords internal
dd_poly_dictionary <- function(x, degree) {
  vapply(0:degree, function(d) x^d, numeric(length(x)))
}

#' A fitted 1D polynomial
#'
#' @param coeffs Numeric coefficient vector, `coeffs[i]` multiplies `x^(i-1)`.
#' @param stderr Optional standard errors, same length as `coeffs`.
#' @return An object of class `poly1d`.
#' @keywords internal
new_poly1d <- function(coeffs, stderr = NULL) {
  structure(list(coeffs = coeffs, stderr = stderr), class = "poly1d")
}

#' @export
predict.poly1d <- function(object, x, ...) {
  dict <- dd_poly_dictionary(x, length(object$coeffs) - 1L)
  as.vector(dict %*% object$coeffs)
}

#' @export
format.poly1d <- function(x, digits = 3, ...) {
  degree <- length(x$coeffs) - 1L
  terms <- character(0)
  for (d in degree:0) {
    c_val <- x$coeffs[d + 1]
    if (c_val == 0) next
    term <- if (d == 0) "" else if (d == 1) "x" else paste0("x^", d)
    if (!is.null(x$stderr) && is.finite(x$stderr[d + 1])) {
      coef_str <- sprintf("(%.*f \u00B1 %.*f)", digits, c_val, digits, x$stderr[d + 1])
    } else {
      coef_str <- sprintf("%.*f", digits, c_val)
    }
    terms <- c(terms, paste0(coef_str, term))
  }
  if (length(terms) == 0) return("0")
  paste(terms, collapse = " + ")
}

#' @export
print.poly1d <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Fit a sparse polynomial to (x, y) data
#'
#' @param x,y Numeric vectors of equal length (predictor and response).
#' @param degree Integer, maximum polynomial degree.
#' @param threshold Numeric >= 0, sparsity threshold (default 0: an ordinary,
#'   non-sparse fit).
#' @param alpha Numeric >= 0, ridge regularization strength.
#' @param weights Optional sample weights.
#' @return A `poly1d` object.
#' @export
dd_fit_polynomial <- function(x, y, degree, threshold = 0, alpha = 0, weights = NULL) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) <= degree) {
    stop("Not enough finite data points (", length(x), ") to fit a degree-",
         degree, " polynomial.", call. = FALSE)
  }
  dict <- dd_poly_dictionary(x, degree)
  fit <- dd_stlsq(dict, y, threshold = threshold, alpha = alpha, weights = weights)
  new_poly1d(fit$coeffs, fit$stderr)
}

#' K-fold cross-validation error for a polynomial fit
#' @keywords internal
dd_cv_error <- function(x, y, degree, threshold, alpha = 0, folds = 5) {
  n <- length(x)
  # Contiguous, unshuffled folds (matching sklearn's KFold(shuffle=False)):
  # the first n %% folds folds get one extra element.
  base_size <- n %/% folds
  extra <- n %% folds
  fold_sizes <- rep(base_size, folds) + c(rep(1L, extra), rep(0L, folds - extra))
  fold_id <- rep(seq_len(folds), times = fold_sizes)
  errs <- numeric(folds)
  for (k in seq_len(folds)) {
    test <- fold_id == k
    train <- !test
    if (sum(train) <= degree || sum(test) == 0) { errs[k] <- NA_real_; next }
    fit <- dd_fit_polynomial(x[train], y[train], degree, threshold = threshold, alpha = alpha)
    pred <- stats::predict(fit, x[test])
    errs[k] <- mean((y[test] - pred)^2)
  }
  mean(errs, na.rm = TRUE)
}

#' Choose a sparsity threshold by cross-validation
#'
#' Fits the polynomial across a range of candidate thresholds and picks the
#' one immediately before cross-validation error starts increasing sharply
#' (an elbow heuristic: as the threshold grows past the point where it starts
#' discarding real, non-noise terms, CV error should jump). This mirrors
#' PyDaddy's `tune_and_fit()`/`model_selection()` behavior.
#'
#' @inheritParams dd_fit_polynomial
#' @param thresholds Optional numeric vector of candidate thresholds. If
#'   `NULL`, a grid is chosen automatically from the magnitude of an
#'   unthresholded fit.
#' @param steps Number of thresholds to try when `thresholds` is `NULL`.
#' @param folds Number of cross-validation folds.
#' @return A `poly1d` object, fit at the chosen threshold.
#' @export
dd_tune_threshold <- function(x, y, degree, thresholds = NULL, alpha = 0,
                               steps = 20, folds = 5) {
  if (is.null(thresholds)) {
    fit0 <- dd_fit_polynomial(x, y, degree, threshold = 0, alpha = alpha)
    max_coef <- max(abs(fit0$coeffs))
    thresholds <- seq(0, max_coef, length.out = steps + 1)[seq_len(steps)]
  }
  cv_errors <- vapply(thresholds, function(th) {
    dd_cv_error(x, y, degree, threshold = th, alpha = alpha, folds = folds)
  }, numeric(1))

  delta <- diff(cv_errors)
  best_threshold <- if (all(is.na(delta))) thresholds[1] else thresholds[which.max(delta)]

  fit <- dd_fit_polynomial(x, y, degree, threshold = best_threshold, alpha = alpha)
  attr(fit, "threshold") <- best_threshold
  attr(fit, "cv_curve") <- data.frame(threshold = thresholds, cv_error = cv_errors)
  fit
}
