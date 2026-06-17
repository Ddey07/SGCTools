# SGCTools cleanup notes

Development copy of `SGCTools` (forked from `Ddey07/SGCTools` @ `594f28d`).
Scope agreed with maintainer: **bug fixes + safe cleanup only** — the public API
(`Kendall_mixed`, `fpca.sgc.lat`, `fromXtoRMixed`, `getLatentPreds`, `sgclm`) and
the numerical estimation paths are unchanged. Every removed object was verified to
have **zero call sites** in `R/` before deletion.

## Bug fixes

1. **Duplicate `Kendall_mixed` (real bug).** A second, internal copy of
   `Kendall_mixed` lived in `R/auxil-functions.R`. Because R collates `R/` files
   in C locale, `auxil-functions.R` loaded *after* `Kendall_mixed.R`, so the
   internal copy **silently overrode the exported, documented one** — meaning the
   shipped `Kendall_mixed()` actually ran without `diag(tau) <- 1` and without the
   input validation promised in the docs. Removed the internal copy; the canonical
   exported version in `R/Kendall_mixed.R` is now the only definition.

2. **Duplicate `lo.dupl`.** Defined twice in `auxil-functions.R`; the first
   definition was dead (immediately overwritten by the second). Removed the first;
   kept the effective second one (used by `sgclm`). No behavior change.

3. **README install command typo.** `devtools::install_github("DDey07/SGCTools")`
   → `"Ddey07/SGCTools"` (wrong capitalization would 404). Fixed in `README.Rmd`
   and `Readme.md`.

## Safe cleanup

4. **Console spam silenced in `fpca.sgc.lat`.** The `nls()`/`nlsLM()` control lists
   had `printEval=TRUE` and `trace=TRUE`, dumping per-iteration output on every call.
   Set both to `FALSE`. **Output only — the fit and all returned values are identical.**

5. **Pruned unused dependencies** (DESCRIPTION / NAMESPACE / `import_package.R`):
   - Removed `tidyverse` (a meta-package; never imported by code and not allowed as a
     hard dependency on CRAN), `refund`, `pcaPP`, and `mnormt` — none are referenced
     anywhere in `R/`.
   - Moved `clusterGeneration` to **Suggests** (used only via `::` in examples).
   - Kept everything actually used: `Matrix`, `fMultivar`, `fda`, `minpack.lm`,
     `mixedCCA`, `mvtnorm`, `numDeriv`, `tmvtnorm`, `splines`, `stats`, `utils`, `Rcpp`.

6. **Removed dead/unused internal helpers** from `auxil-functions.R` (no call sites):
   `g`, `gprime`, `truncate`, `Kendall_bin`, `ties`, `bridge_bc`, `fb`, `fb.sp.uni`,
   `entropy`, `beta.elim2`, `deriv.b2`. Also removed an unused `ids` assignment and a
   stale commented line inside `conc_ties`.
   - Note: `beta.elim2` / `deriv.b2` were stubs for multivariate-coefficient
     asymptotics. They are unused here and recoverable from git history if the new
     multivariate package needs them.

7. **Robustness:** `fpca.sgc.lat` now validates `type` up front
   (`"cont"`/`"bin"`/`"ord"`/`"trunc"`); previously a typo produced a confusing
   "object 'eunsc' not found" error deep in the function. Valid calls are unaffected.

## Not changed (deliberately)

- The random NLS start `init <- runif(...)` is left as-is (adding a seed would change
  numerical output). Consider exposing an optional `seed=` argument later.
- `getFromNamespace("fromKtoR_mixed"/"bridgeF_tt", "mixedCCA")` is kept but is
  fragile across mixedCCA versions — flagged for a future hardening pass.
- `man/*.Rd` were left untouched; regenerate them with `devtools::document()`
  (see below), which will also rewrite `NAMESPACE` from the roxygen tags.

## Verify locally (R + devtools)

```r
# from the SGCTools_dev/ directory
devtools::document()      # regenerate NAMESPACE + man/ from roxygen (recompiles Rcpp)
devtools::load_all()      # quick load to catch parse/namespace errors
devtools::check()         # full R CMD check
```

Static checks already performed here: all five exported functions resolve to a single
definition, every helper they call is present, brace/paren/bracket balance is preserved
in all edited files, and no removed name has a remaining call site.
