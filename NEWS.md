# SGCTools 1.1.5

## Bug fixes

* `Kendall_mixed()` (and `conc_ties()`, used by `sgclm()`) undercounted concordant and
  discordant pairs. The insertion-sort step of the merge-sort pair counter, used on
  chunks of fewer than 10 observations, credited only the moving element and never the
  element it passed. The estimate was pulled toward zero, by about 15% at 10-15
  pairwise-complete observations, 3-6% at 24-30, and under 1% at 100 or more. With 10
  perfectly concordant observations it returned 0.778 instead of 1. It now equals
  Kendall's tau_a exactly (ties count 0 and stay in the denominator), matching a
  brute-force computation and `mixedCCA::KendallTau()`.
* `bridgeF_tt()` and `fromKtoR_mixed()` have explicit arguments again. The
  `function(...)` wrappers introduced in 1.1.4 broke
  `Vectorize(bridgeF_tt, vectorize.args = c("r", "zratio1", "zratio2"))`, so
  `fpca.sgc.lat()` failed on truncated components.
* New tests: `tests/testthat/test-kendall.R`.

# SGCTools 1.1.4

## Bug fixes

* `getLatentPreds()` now builds its cutoffs and marginal transformations with the
  `type` vector you supply. Previously the internal `preprocess_data()` was called
  without `type` and re-guessed each column's type from the data, which could
  disagree with `type` for sparsely observed columns.
* Ordinal cutoffs are now built over the full level set `0, 1, ..., K` (K = largest
  observed ordinal code), not only over the levels observed in each column. When a
  sparse ordinal column missed a level, every later cutoff shifted and the
  truncated-normal step received an `NA` bound ("'upper' not specified or contains
  NA"), so the whole call failed. An absent level now gets a zero-width interval.
  Wherever the previous code ran and every level 0..K was observed in a column,
  the results are identical. A column that lacked level 0 could previously receive
  the bounds of the wrong category; this is now correct.
* A column with no observed value no longer triggers an error in the marginal
  transformation estimate.
* Input checks: `type` must have length `ncol(X)` and use "cont", "trunc", "ord",
  "bin"; ordinal codes must be non-negative integers; binary codes must be 0/1.

## New

* `getLatentPreds()` gains `pmv_maxpts` and `pmv_abseps` (Monte-Carlo controls of the
  truncated-normal conditional mean). The defaults reproduce the previous precision.
* Unit tests (`tests/testthat`).

# SGCTools 1.1.3

* See `CLEANUP_NOTES.md`.
