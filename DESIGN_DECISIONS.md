# Design decisions: how daddyR differs from PyDaddy, and why

This is a from-scratch, native R reimplementation of PyDaddy's scalar (1D)
workflow, not a line-by-line translation. Every deviation from PyDaddy's
behavior is listed here, with a citation to the PyDaddy source it changes and
the reasoning. References are to PyDaddy v1.1.1
(https://github.com/tee-lab/PyDaddy).

## Confirmed bugs in PyDaddy that this port does not reproduce

**1. `_remove_nan()` crashes on any numpy released since 2020.**
`pydaddy/preprocessing.py:38` calls `np.linspace(..., dtype=np.int)`.
`np.int` was deprecated in NumPy 1.20 and removed entirely in NumPy 1.24
(2023); this line raises `AttributeError: module 'numpy' has no attribute
'int'` on any current numpy install (confirmed against numpy 2.4.4 while
building this port).

**2. The Gaussianity test crashes on any numpy released since 2021.**
`pydaddy/analysis.py`'s `_get_critical_values()` calls
`np.histogram(kl_dist, normed=True)`. The `normed` argument was removed from
`numpy.histogram` in NumPy 1.24 (replaced by `density`); this raises
`TypeError` on current numpy.

**3. Bin edges used for grouping data don't match the reported bin
locations.** `pydaddy/characterize.py`'s `Main._preprocess()` (via
`_validate_inputs`) sets `inc = (max - min) / bins`, the width used to decide
which points fall in which bin (`sde.py`'s `_drift_and_diffusion`, condition
`X_ < (b + inc)`). But the bin *locations* it reports (`op`, used both as the
x-axis for every plot and as the fitting input) are
`np.linspace(min, max, bins)`, whose point-to-point spacing is
`(max - min) / (bins - 1)` — a different number whenever `bins > 1`. Since
`inc < spacing`, each bin of width `inc` starting at `op[i]` leaves a small
gap before `op[i+1]` starts; points that fall in that gap are silently
excluded from every bin's average. Confirmed on the bundled
`scalar_pairwise` dataset with the default `bins=20`: `inc = 0.1` but the
`op` spacing is `2/19 ≈ 0.1053`, a >5% gap compounding across the range.
daddyR's `dd_bin_edges()` constructs `bins + 1` edges that exactly tile the
data range and reports bin midpoints, so coverage is exact (see
`test-estimate.R`'s bin-coverage regression test). **Consequence:** daddyR's
binned drift/diffusion averages will not match PyDaddy's bit-for-bit (they're
computed over a very slightly different, corrected partition) — see
`test-reference-values.R` for how this is validated instead (the
bin-independent raw drift series and polynomial fit coefficients are checked
for exact agreement; the binned averages are checked for close numerical
agreement, not identity).

**4. The documented public API can't reach PyDaddy's own auto-timescale
logic.** `pydaddy/characterize.py`'s internal `Main` class computes
`Dt = ceiling(autocorrelation_time / 10)` whenever `Dt` is left `None` — this
is real, exercised code (used by the interactive timescale-slider feature).
But the actual public, documented entry point, `pydaddy.Characterize()`,
hardcodes `Dt=1` as its own default and always passes it through, so a user
calling the documented API can never trigger the auto-selection without
reaching into the undocumented `Main` class directly. daddyR's `dd_analyse()`
uses `Dt = NULL` (meaning "auto-select") as its actual default, since that
appears to be the intended behavior. Pass `Dt=` explicitly to get PyDaddy's
literal default of a fixed lag.

**5. `fit()` can silently misalign `x` and `y` when `Dt != dt`.**
`pydaddy/daddy.py`'s `Daddy.fit()` always slices the predictor as
`x = self._ddsde._X[:-self.Dt]`, even when fitting the diffusion function
`'G'`, whose target values (`self._ddsde._diffusion_`) were computed using
`dt`, not `Dt`. When `Dt == dt` (the overwhelmingly common case in the
tutorials and, we'd guess, most usage) this is harmless; when they differ,
the two arrays are misaligned by `|Dt - dt|` samples without any error or
warning. daddyR's `dd_fit()` always pairs each series with the x-values
`dd_estimate()` computed alongside it (`x_for_drift_fit` with
`drift_series`, `x_for_diff_fit` with `diff_series`), so this can't happen
regardless of whether `Dt` and `dt` differ. Relatedly, `dd_analyse()` now
defaults `dt` to whatever `Dt` resolves to (PyDaddy defaults them
independently, `Dt` from autocorrelation and `dt` to a flat `1`), which
avoids the whole class of `Dt != dt` edge cases by default; pass both
explicitly to control them independently.

## Statistical / numerical decisions

**6. The Gaussianity test is not actually a Kolmogorov-Smirnov-worthy KL
divergence, so it's been replaced with a standard test.**
`pydaddy/metrics.py`'s `_kl_divergence(p, q)` computes
`sum(p * log(abs((p + eps) / (q + eps))))` applied *elementwise to raw
samples* `p` and `q` — not to density or probability-mass estimates of `p`
and `q`, which is what the KL divergence formula actually requires. Applied
to raw draws from `N(0,1)` (which can be negative, hence the `abs()`), this
doesn't correspond to any standard information-theoretic quantity, and the
result depends on the arbitrary pairing/order of samples in `p` and `q`
rather than only on their distributions. (Item 2 above means this code path
also can't currently run at all on modern numpy.) daddyR's
`dd_gaussianity_test()` instead runs a one-sample Kolmogorov-Smirnov test of
the standardized residuals against `N(0,1)` — it answers the same question
("is this noise Gaussian?") with an established, well-calibrated procedure.
This is a genuine behavior change, not a refactor, and is worth your
independent sign-off given it affects a diagnostic result. Simple
alternatives if you'd prefer a bootstrap-flavored test instead of an
asymptotic one: an Anderson-Darling test, or a proper bootstrap using a
real (density-based) Jensen-Shannon divergence.

**7. Standard errors for sparse polynomial fits no longer risk crashing on
near-singular matrices.** `pydaddy/fitters.py`'s `PolyFitBase.fit()` computes
coefficient standard errors via
`np.linalg.inv(dictionary.T @ dictionary)`, inverting the *full* p x p
dictionary matrix even when most columns were thresholded to zero.
Polynomial dictionaries (`[1, x, x^2, x^3, ...]`) are highly collinear, so
this matrix is frequently ill-conditioned and occasionally exactly singular
past degree 3-4, either raising `LinAlgError` or returning wildly inflated
standard errors for terms that aren't even in the model. daddyR's
`dd_stlsq()` inverts only the submatrix of active (kept) terms — which is
what the standard-error formula is actually describing — and returns `NA`
with a warning (rather than an error) if that smaller matrix is itself
singular.

**8. Cross-validation folds now match scikit-learn's actual `KFold`
semantics.** An early draft of `dd_cv_error()` used interleaved fold
assignment; this was corrected to contiguous, unshuffled blocks (the first
`n %% folds` folds get one extra sample), matching
`sklearn.model_selection.KFold(shuffle=False)`, which is what PyDaddy's
`_get_cv_error()` uses.

## Removed dead code

**9. The automatic polynomial-order-selection machinery was never actually
wired up, so it isn't ported.** `pydaddy/preprocessing.py` defines
`_order()`, `_find_order()`, `_o1()`, `_o2()`, `_get_o1_o2()`, and
`_r2_vs_order_multi_dt()` — a fairly involved scheme for auto-detecting
polynomial degree from R²-vs-order curves at multiple timescales. But
`Main._preprocess()` (the only place preprocessing actually happens) calls
only `_validate_inputs()`; the entire order-detection block is commented
out. `_find_order()` itself also contains dead code after an unconditional
early `return` statement. None of this is reachable from any documented
PyDaddy API. daddyR does not port it. (A future phase could implement
degree selection properly via cross-validated BIC across degrees, which
would be a real feature rather than a resurrection of unreachable code.)

## API and structure cleanups (behavior-preserving)

**10. No silent keyword-argument absorption.** Every PyDaddy mixin
(`SDE`, `AutoCorrelation`, `Preprocessing`, `Metrics`, `Visualize`, `Daddy`,
`Main`, ...) initializes with `self.__dict__.update(kwargs)`, meaning any
keyword argument — including a misspelled one — is silently accepted and
stored as an attribute, whether or not anything reads it. daddyR's functions
take explicit named parameters with no catch-all `...` where it isn't
needed; passing an unrecognized argument is an immediate, clear R error
(`unused argument`) rather than a silent no-op (see
`test-analyse.R`'s "rejects unknown/misspelled arguments" test).

**11. No global warning suppression.** `pydaddy/characterize.py:13` runs
`warnings.filterwarnings("ignore")` at import time, suppressing every
warning the package (or anything it calls) might raise, package-wide, for
the rest of the Python process. daddyR raises targeted warnings only at the
specific points where something notable but non-fatal happened (e.g. the
singular-stderr fallback in item 7).

**12. No bare exception handlers.** Several places in PyDaddy
(`characterize.py:364`, `daddy.py:1110`, `visualize.py:1259,1510`) use bare
`except:` clauses, which catch and silently discard *everything*, including
`KeyboardInterrupt` and genuine programming errors, not just the specific
failure being anticipated. daddyR's few `tryCatch()` calls catch specific,
documented failure modes only.

**13. Flattened class hierarchy.** PyDaddy's classes form a nine-file,
multiple-inheritance mixin chain (`SDE` → `AutoCorrelation`/
`UnderlyingNoise` → `GaussianTest` → `Preprocessing` → `Main`/`Characterize`
→ `Daddy` → `Visualize`, with `Metrics` mixed into several of these), so
finding which class actually defines a given `self._method()` called from
deep in `Daddy` often means checking several files. daddyR is a flat set of
plain functions grouped by what they do (`estimate.R`, `poly_fit.R`,
`diagnostics.R`, `preprocess.R`), each independently callable and testable
without constructing a large stateful object first.

**14. Vectorized bin lookups.** PyDaddy assigns each point to a bin with a
Python `for` loop calling `np.argwhere(x < bins)[0][0]` once per data point
(`sde.py`, `analysis.py`). daddyR uses a single vectorized `findInterval()`
call. Besides being simpler, this matters in practice: the datasets this
package is meant for (animal-tracking and cell-migration time series) run to
10^5-10^6 points, where a per-point Python loop is a real performance cost.

**15. Functional (not mutating) `dd_fit()`.** `pydaddy/daddy.py`'s
`Daddy.fit()` mutates the object in place (`setattr(self, function_name,
res)`) as well as returning the fit. R data structures are not reference
objects by default, so daddyR's `dd_fit()` returns an updated *copy* of the
`daddy` object (`dd <- dd_fit(dd, "drift", degree = 3)`) — the idiomatic R
pattern for "compute a derived result and attach it," and one that doesn't
surprise a user by mutating something they didn't reassign.

**16. No console-script entry point.** PyDaddy ships a `pydaddy` console
script (`pydaddy/__console__.py`). This isn't a typical pattern for R
packages, whose users work through function calls in a script or console
session; it has been dropped in favor of the plain function-based API
(a Shiny app for interactive exploration is a possible future addition,
tracked in the porting roadmap, rather than a command-line tool).

## Phase 2: vector (2D) data and SDE simulation

PyDaddy's vector code path (mostly `sde.py`'s `SDE._vector_drift_diff` and
`daddy.py`'s `Daddy.simulate()`) has its own set of idiosyncrasies, beyond
the ones already fixed for the scalar case above.

**17. Vector binning has the scalar bin-edge bug (item 3), plus a
closed-interval double-counting bug of its own.** PyDaddy's 2D bin-edge
width (`inc_x`/`inc_y`) is computed as `(max-min)/bins` while the reported
bin locations use `np.linspace(..., bins)` — the same mismatch as item 3,
now confirmed on both axes. Separately, its 2D bin-membership test is
*closed on both ends* (`bin_x <= x & x <= bin_x + inc_x`), unlike its own
scalar test (`X_ < b+inc`, half-open) — so a point sitting exactly on a
shared edge between two bins is counted in **both**. `dd_estimate_vector()`
reuses [`dd_bin_index()`](R/estimate.R)'s half-open, exactly-tiling bins for
both axes, fixing both issues the same way item 3 fixed them for the scalar
case.

**18. Binned-average matrices are filled and read with inconsistent
transpose conventions.** PyDaddy declares `avgdriftX`/`avgdriftY`/etc. as
shape `(nx, ny)` but fills them `[y_bin, x_bin]` (the transpose of what the
declared shape suggests) — self-consistent within `sde.py` (the lookup
during diffusion estimation uses the same convention), but a *different*
diagnostic function elsewhere in PyDaddy (`Daddy.noise_diagnostics()`'s
residual computation) indexes the same arrays `[x_bin, y_bin]`, silently
reading the wrong bin's average drift whenever a point's x- and y-bin
indices differ. `dd_estimate_vector()` uses one explicit convention
everywhere — every matrix it returns is indexed `[bin_x, bin_y]` — so there
is no second, inconsistent reader to get out of sync.

**19. The diffusion matrix estimator conflates the drift and diffusion
timescales.** PyDaddy's vector diffusion formula lags the raw increment
by `Dt` (the *drift* timescale) but divides by `dt` (the *diffusion*
timescale) — structurally different from its own scalar fast-mode formula,
and only numerically equivalent to it when `Dt == dt`. (A related, purely
mechanical bug: `Daddy.fit()` slices its covariate matrix by `self.dt` when
fitting `G11`/`G22`/`G12`, but the target arrays it's pairing that with were
actually built with length `N - Dt`, not `N - dt` — a shape mismatch that
raises an error from the regression call whenever `Dt != dt`.)
`dd_estimate_vector()` instead extends the *scalar* residual formula (each
point's own bin's average drift, scaled by `t_int` only) to two dimensions,
then lags the resulting residual series by `dt` for `G11`, `G22`, *and*
`G12` — keeping `Dt` and `dt` independently meaningful, exactly as in the
scalar case, and avoiding the alignment bug entirely since every fit pairs
`x`/`y` arrays daddyR computed together (same fix as item 5). Note this
formula happens to numerically coincide with PyDaddy's own vector formula
whenever `Dt == 1` (the common case, and what the reference-value tests
use), since the two `Dt`-related terms cancel — the difference only shows up
once `Dt > 1`.

**20. No separate "G21".** Because a real diffusion matrix built this way is
symmetric by construction (`G12` and `G21` are the same
`residual1 * residual2` product), PyDaddy nonetheless carries "G21" as a
distinct fittable quantity throughout its API and simply forces it to equal
whatever `G12` was fit to (`daddy.py`'s `fit()`, `G21 = G12` whenever
`'G12'`/`'G21'` is requested) — a redundant concept that only creates a
chance to set one without the other. `dd_fit()` doesn't expose a `"G21"`
option at all; `G12` is used for both off-diagonal entries wherever the
full matrix is needed (e.g. [`dd_simulate()`]'s vector path).

**21. No safety check before taking the diffusion matrix's square root
during simulation.** `Daddy.simulate()` calls `scipy.linalg.sqrtm(G)`
directly, reacting only if that call happens to raise `LinAlgError` or
`ValueError` — but `sqrtm` isn't guaranteed to raise on an indefinite
matrix; it can silently return a matrix with a nonzero imaginary part, which
then propagates into the simulated trajectory as `NaN`/complex values with
no clear error. This is a real risk here, not a hypothetical one: the fitted
`G11`/`G22`/`G12` are unconstrained polynomials that can evaluate to an
invalid (indefinite) matrix once a simulated trajectory drifts outside the
range of the original training data. (PyDaddy also only clamps the diagonal
terms with `abs()` in its no-cross-term fallback branch, not in the
general `sqrtm` branch — an inconsistency of its own.) `dd_matrix_sqrt_psd()`
always computes the symmetric square root via eigendecomposition and
explicitly clips any negative eigenvalue to zero (with a warning when this
actually happens), so [`dd_simulate()`] always produces a finite, real
trajectory, uniformly for the diagonal-only and cross-diffusion cases alike.

**22. One integrator, used consistently.** PyDaddy's `simulate()` uses a
higher-order stochastic Runge-Kutta scheme (`sdeint.itoSRI2`) for scalar
data but plain Euler-Maruyama (`sdeint.itoEuler`) for vector data — an
asymmetry with no documented rationale. `dd_simulate()` uses Euler-Maruyama
for both, via the shared `dd_matrix_sqrt_psd()`-based noise construction;
this is simpler and consistent, at the cost of Euler-Maruyama's lower strong
convergence order (use a smaller `t_int` than you might otherwise need for
high-precision work).

**23. Standard errors for the 2D polynomial fits get item 7's fix "for
free."** PyDaddy's `PolyFitBase.fit()` (shared, unmodified, between its 1D
and 2D fitters) inverts the *full* candidate dictionary rather than the
sparsified active-term submatrix when computing coefficient standard
errors — the same bug as item 7, and if anything more likely to bite for
the vector case, since a 2D polynomial dictionary has more, more collinear
terms. Because daddyR's `dd_fit_polynomial2d()` reuses the same
`dd_stlsq()` used by the scalar `dd_fit_polynomial()`, it already inverts
only the active submatrix (with the same `NA`-and-warning fallback on
singularity) — no vector-specific fix was needed here.

**24. Gaussianity testing exposed for the vector case too (and a related
scalar fix).** `dd_estimate_vector()` returns `noise1`/`noise2`: the
pre-squared, drift-corrected residual series for each component — the
statistically appropriate input to [`dd_gaussianity_test()`] (as opposed to
`diff11_series`/`diff22_series`, which are already squared and so should not
themselves be tested for normality). `dd_gaussianity_test_vector()` is a
thin convenience wrapper running the test on both components.
`dd_estimate()`'s scalar output gained the analogous `noise_series` field at
the same time, for the same reason — this was a latent issue in the 0.1.0
release (nothing there was actually *wrong*, since `dd_gaussianity_test()`
takes whatever vector it's given, but the demo/documentation's suggested
usage tested the already-squared `diff_series`). Whether the *test itself*
(item 6's KS-test replacement) is the right choice is still open for your
review; this item is only about testing the right quantity, whichever test
is used.

## Explicitly out of scope

Cross-validated automatic degree selection (see item 9), the full
interactive diagnostic suite (`model_diagnostics()`'s simulate-and-re-
estimate self-consistency plots, the noise-autocorrelation and residual
QQ-plot panels, and the interactive timescale-slider views), and a
console/Shiny UI are not in this port yet. Nothing above should be read as
claiming feature parity; it's a list of deliberate differences within the
scope covered so far.
