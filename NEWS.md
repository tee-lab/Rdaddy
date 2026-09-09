# daddyR 0.1.0

Initial native R port of PyDaddy (Phase 1: scalar time series only).

## New features

* `dd_analyse()`: estimate binned drift and diffusion (Kramers-Moyal
  coefficients) from a scalar time series.
* `dd_fit_polynomial()`: sparse polynomial fitting (sequentially-thresholded
  ridge regression) for the drift and diffusion functions, with optional
  cross-validated threshold selection via `dd_tune_threshold()`.
* `dd_autocorrelation()` and `dd_gaussianity_test()` diagnostics.
* `print.daddy()`, `summary.daddy()`, and `plot.daddy()` S3 methods.
* Bundled example dataset `scalar_pairwise` (a simulated pairwise-interaction
  model time series, from the original PyDaddy package).

## Relative to PyDaddy

This release covers PyDaddy's scalar (1D) workflow only. Vector (2D) data,
cross-diffusion, SDE simulation, and the full diagnostic/plotting suite are
planned for later releases. See `DESIGN_DECISIONS.md` for a full list of
behavioral differences from PyDaddy and the reasoning behind each one.
