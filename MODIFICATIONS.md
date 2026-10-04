# Pooled SDR and TMLE development branch

This modification starts from Herb Susmann's current GitHub `curve` branch,
commit `f33c53174b29d881916022df178029b30e7275f6` (2025-12-15).
The package version is `1.6.0.9001`; the local branch is
`pooled-tmle-upstream-fixes`.

## Use

```r
fit <- lmtp_curve(
  data, trt, outcome, baseline = baseline, shift = policy,
  estimator = "tmle", folds = 5
)
```

`estimator = "sdr"` remains the default. The pooled learner backend now follows
upstream `SuperLearner`, with `SL.*` candidate names. The common Herb aliases
`glm`, `mean`, `ranger`, and `cv_glmnet` are translated to their SuperLearner
counterparts. Other mlr3 candidate names require a corresponding SL wrapper.

## Pooled TMLE

Both methods use one pooled initial nuisance learner per recursive depth, with
time as a predictor. For TMLE, the training recursion uses the next shifted
regression prediction as its response. The entire initial recursion finishes
inside the outer training set before targeting begins.

For outcome horizon `s`, regressions are targeted backward over `r = s, ..., 1`.
Each target uses a weighted logistic intercept fluctuation on the outer validation
fold, with response equal to the next targeted shifted regression (the observed
outcome at the terminal step), offset equal to the current natural regression,
and weights equal to the cumulative treatment/censoring ratio times sampling
weights. The terminal step also includes the sporadic-measurement weight.

Separate intercepts for each `(s, r)` solve the score equation for every curve
point. A single intercept across all outcome horizons would only solve a pooled
equation. The intercept is solved directly by the same logistic score used by
upstream TMLE; boundary responses use its limiting values. The TMLE estimate is
the mean of the targeted first shifted regression. `fit$targeting` records
coefficients and residual scores.

The held-out fluctuation equation follows upstream commit `9996b04`. Unlike
interleaving validation targeting and subsequent initial regressions, the two
stages here prevent validation responses from entering later nuisance training.
Training predictions used internally by the initial outcome recursion come from
fits within the outer training set; reported nuisance predictions are out of
sample with respect to the outer validation fold.

As in Herb's branch, survival estimates receive a monotonicity projection.
`unprojected_estimates` records the pre-projection estimates to which the TMLE
score equations apply. The existing SDR isotonic pseudo-outcome calibration is
retained for SDR and is not applied to TMLE.

## Shared splits

`folds` accepts either a count or a list of
`list(training_set = ..., validation_set = ...)` entries using original wide-row
indices. Training sets must complement validation sets; validation sets must
partition all rows. A cluster cannot appear on both sides of an outer split.

`learner_folds` accepts one positive-integer label vector per outer fold, in that
fold's training-row order. All rows for a cluster must have the same label.
These labels control treatment, sporadic-measurement, and outcome learner CV,
including repeated copies in density-ratio stacks and pooled person-time rows.
Each fit restricts the supplied labels to its eligible rows. At least two labels
must remain. Without supplied labels, a common inner split is generated using the
larger learner fold count, with a minimum of two.

The fitted object returns the outer and inner plans for reuse. `folds = 1` is
explicitly an in-sample analysis. SDR training pseudo-outcome weights are refit
with inner cross-fitting wholly within each outer training set, rather than
borrowing global out-of-fold fits that may have trained on the outer validation
fold. SuperLearner CV inside those additional inner nuisance fits is restricted
to their own training data.

Continuous responses use explicit `bounds` when supplied, including the same
scaling for training and validation. Bridge-transformed responses should be
treated as continuous and require appropriate bounds. The bridge, adjoint
estimation, and combined one-step correction are separate integration work.

## Upstream fixes carried over

Sources were fetched from current upstream `master`
(`f04a2b47f46debc515ce4ae778e05ebfde922c44`) and `devel`
(`0631f120216e875ec41e6a24b6305a1778a66804`).

| Source commit | Fix or update |
| --- | --- |
| `49708e9` | History-dependent shifts read natural treatment history. Cherry-picked. |
| `aac4880` | Multivariate policy-following indicators are calculated by subject. Cherry-picked. |
| `e68d3a7` | Supplied shifted-data validation flattens multivariate treatment names. |
| `48ddf8e` | Explicit outcome bounds are passed to and used by the task. |
| `d56fec0` | Correct internal cluster ID column for outer folds; stratification uses the final outcome in the pooled branch. |
| `9996b04` | Ordinary TMLE targeting uses validation responses and validation weights. |
| `7610633` | Local fold generator replaces origami, with additional singleton-stratum and empty-fold safeguards. |
| `f04a2b4` | Learner CV fold count cannot exceed the number of training clusters. |
| `e48107b` | Risk-ratio and odds-ratio inference uses the log scale. Cherry-picked. |
| `6f0127a` | `mtp = TRUE` default. |
| Upstream `R/sl.R` | SuperLearner adapter, `onlySL = TRUE`, and prediction feature selection. |
| `0631f12` | Point treatment can have a longer time-varying censoring horizon. |

Additional integration repairs declare the previously missing `OrdMonReg`
dependency, require `ife >= 0.2.0`, correct the prediction-bound default, encode
sporadically missing lagged outcomes with missingness indicators, avoid treating
a repeated terminal outcome as an observed lag, and handle final sporadic
measurement and absorbed survival states. Cumulative EIF weights retain their
observation dimension even in small validation folds. A point-treatment policy
is applied only at its treatment time, rather than being applied repeatedly
during later curve and fixed-horizon regression predictions.

## Validation

The tests cover known-truth randomized treatment curves, agreement between pooled
and ordinary TMLE for a single outcome time, every held-out score equation with
sampling weights, censored survival curves, multivariate treatments, missing
outcomes, explicit continuous bounds, grouped sample splits, and invariance of a
fold's initial nuisance fits to changes in its outer validation data. Legacy
estimator fidelity and argument tests are retained.

The missing isotonic dependency was installed into `.library`, a private project
library. The existing installed lmtp package was not used as an implementation
base or replaced. Build and check outputs are in the ignored `validation/`
directory. Final validation on 2026-10-03 with R 4.5.2:

- Source test suite: passed, with no warnings.
- Built-package tests: **275 passed, 0 failed, 0 warnings, 0 skipped**.
- `R CMD check --no-manual --no-build-vignettes`: **Status: OK**, with no errors,
  warnings, or notes. Manual and vignette builds were excluded.
- Weighted TMLE residual score tests use a `1e-9` tolerance at every tested
  outcome horizon and regression time.
- Seeded known-truth tests cover randomized binary curves (maximum error below
  `0.035`) and a continuous point-treatment MTP (maximum error below `0.025`
  across both estimators and all horizons).

Reproduce from the repository root:

```sh
R_LIBS="$PWD/.library" Rscript -e 'devtools::test(".", reporter="summary", stop_on_failure=TRUE)'
R_LIBS="$PWD/.library" R CMD build . --no-build-vignettes --no-manual
NOT_CRAN=true R_LIBS="$PWD/.library" _R_CHECK_FORCE_SUGGESTS_=false R CMD check lmtp_1.6.0.9001.tar.gz --no-manual --no-build-vignettes --output=validation
```

The installable source archive is `lmtp_1.6.0.9001.tar.gz`.
The check log is `validation/lmtp.Rcheck/00check.log`; the built-package test
summary is `validation/lmtp.Rcheck/tests/testthat.Rout`.
