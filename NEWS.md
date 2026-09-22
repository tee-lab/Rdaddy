# daddyR 0.2.0

Phase 2: vector (2D) time series support, SDE simulation, and an additional
diagnostic.

## New features

* `dd_analyse_vector()`: estimate the drift vector (`F1`, `F2`) and diffusion
  matrix (`G11`, `G22`, `G12`) from a bivariate (2D) time series, on the same
  `"daddy"` S3 class as `dd_analyse()` (with `$vector` set to `TRUE`).
  `dd_fit()`, `dd_drift()`, `dd_diffusion()`, `dd_export_data()`, and the
  `print`/`summary`/`plot` methods all dispatch automatically based on
  `dd$vector`, so the same function names work for both scalar and vector
  objects.
* `dd_fit_polynomial2d()` and `dd_tune_threshold2d()`: sparse 2-variable
  ("poly2d") polynomial fitting for the vector case, sharing the same
  sequentially-thresholded ridge regression machinery as the scalar
  `dd_fit_polynomial()`/`dd_tune_threshold()`.
* `dd_simulate()`: simulate a synthetic trajectory (Euler-Maruyama) from a
  fitted scalar or vector model -- useful as a generative model, and for
  self-consistency checks (compare a re-`dd_analyse()`d simulated series
  against the original fit).
* `dd_gaussianity_test_vector()`: convenience wrapper running
  `dd_gaussianity_test()` on both extracted noise components of a vector fit.
* `plot.daddy()` now produces a phase-portrait-plus-heatmap figure for
  vector objects (drift components `F1`/`F2` and diagonal diffusion terms
  `G11`/`G22`, binned over the `(x1, x2)` plane).
* Bundled example datasets `vector-pairwise`/`vector-ternary` (subsampled
  from PyDaddy's own bundled vector model data), loadable via
  `dd_load_sample_data()`.
* `dd_estimate()`'s (scalar) output gains a `noise_series` field: the
  pre-squared, drift-corrected residual series, which is the statistically
  correct input to `dd_gaussianity_test()` (as opposed to `diff_series`,
  which is already squared and so should not itself be tested for
  normality).

## Deviations from PyDaddy's vector implementation

PyDaddy's own vector code path (`SDE._vector_drift_diff`, `Daddy.simulate()`)
has several idiosyncrasies beyond the bin-edge issue already fixed for the
scalar case in 0.1.0 -- a closed-interval binning test that double-counts
points on shared bin edges, binned-average matrices that are filled
transposed relative to their declared shape (and read back with yet another,
inconsistent transpose convention in one diagnostic function), a diffusion
estimator that conflates the `Dt` and `dt` timescales, no positive-
semidefiniteness check before taking the diffusion matrix's square root in
`simulate()`, and an inconsistency where the scalar and vector simulators use
different-order numerical integrators. All are fixed in daddyR's vector
implementation; see `DESIGN_DECISIONS.md` and `dd_estimate_vector()`'s /
`dd_matrix_sqrt_psd()`'s documentation for the full details of each.

## Relative to PyDaddy

The full diagnostic/plotting suite (interactive sliders, `model_diagnostics()`
self-consistency plots) is still not ported; `dd_simulate()` covers the
generative/simulation half of that gap. See `DESIGN_DECISIONS.md` for the
complete list of behavioral differences from PyDaddy and the reasoning
behind each one.

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
