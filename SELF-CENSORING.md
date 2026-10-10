# lmtp self-censoring branch

Repository: https://github.com/idiazst/lmtp; branch: `self-censoring`.
Upstream and Herb's Git histories remain in the branch ancestry. The standard
upstream branches are unchanged. This organization does not port any additional
upstream changes or alter the estimators.

```r
remotes::install_github("idiazst/cmbridge")
remotes::install_github("idiazst/lmtp", ref = "self-censoring")
```

- `R/bridge_curve.R`: observed-history bridge/adjoint fits, `lmtp_bridge()` and
  the one-step correction for either pooled SDR or logistic TMLE.
- `R/bridge_response.R`: sequential responses rebuilt inside training samples.
- `R/sample_splits.R`, `R/make_folds.R`: supplied outer and learner splits,
  with grouped copies of each person kept together.
- `R/sl.R`, `R/saturated_learner.R`: SuperLearner integration and training-only
  penalty selection; `R/curve.R` and `R/curve_tmle.R`: pooled sequential fits.
- `tests/testthat`: upstream tests and integration tests covering targeting,
  out-of-sample nuisance fits, supplied splits and shared penalty selection.

Run `testthat::test_local(stop_on_failure = TRUE)` from this package directory.
Generic moment estimators and ensembles are tested in cmbridge. The simulation
repository pins the exact versions used for the saved study; it does not silently
replace those dependencies with the newest branch head.
