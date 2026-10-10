# Self-censoring branch of lmtp

The package remains **lmtp**, version 1.6.0.9021, on branch `self-censoring`
in the idiazst/lmtp fork. Generic moment estimators are supplied by cmbridge.
Simulation code and latest results live in lmtp_selfcensor_sim.

## Foundation

This implementation descends from Herb Susmann's pooled `curve` branch,
commit f33c53174b29d881916022df178029b30e7275f6. Its Git ancestry is retained.
The previous implementation work ported selected upstream lmtp fixes concerning
natural treatment histories, multivariate shifts, censoring, bounds, grouped
folds, outcome dimensions and log-scale relative-risk/odds-ratio inference.
The repository reorganization does not merge additional upstream changes.

## Sequential estimators and correction

`lmtp_curve(estimator = "sdr")` and `lmtp_curve(estimator = "tmle")` use
SuperLearner and pooled outcome regressions. Logistic TMLE first completes the
initial regression recursion within the outer training sample, then targets
separate held-out scores for each endpoint and recursive time. It records the
fluctuation coefficients and score residuals. SDR retains its recursive update.

`lmtp_bridge()` adds observed-history self-censoring bridge/adjoint estimation
and the one-step correction to either sequential estimator. All three bridge
candidates use inverse-expit parameterization during estimation. Adjoint
candidates retain unrestricted identity links. The candidate libraries are sieve
minimum distance, Landweber and PMMR; saturated L1 remains a SuperLearner
regression/classification candidate, rather than a conditional-moment candidate.

## Shared training samples

Outer `folds` accepts a count or supplied complementary training/validation sets.
`learner_folds` supplies person-level labels within each outer training sample.
Grouped copies in classifier stacks and pooled person-time rows stay together.
Nested training callbacks retain the shared person assignments for SuperLearner,
cmbridge ensemble selection, sequential responses and penalty selection.

`R/bridge_response.R` rebuilds upstream nuisances and backward responses inside
the training sample needed for each fit, excluding its validation outcomes.
SDR/TMLE can reuse bridge/treatment fit caches only for matching data, libraries,
splits and kernel; regression responses remain estimator specific.

## Scoring and penalties

Conditional-cell scoring uses cmbridge's distinct-pair U-statistic. Gaussian and
linear scores also remove self-products. The raw ensemble Gram is projected to
positive semidefinite before simplex optimization; penalty losses remain raw.
Positive sieve/L1 penalty grids are evaluated on training-only CV splits.
Boundary minima are recorded, grids extend, and unresolved boundaries fail.
Incomplete/nonfinite trials cannot be selected; candidate failures are retained.
Treatment predictors remain in the penalty designs.

cmbridge fitted functions are not clipped after optimization. The integration's
approved wider application interval is [-100,100]. Logistic TMLE is unchanged;
no Gaussian mean update is introduced by this branch.

## Verification

The package retains upstream unit tests and adds pooled-TMLE, upstream-fix,
bridge-curve and sequential-response integration tests. They check score
reconstruction, validation exclusion, shared person splits, penalty CV and
inverse-expit bridges/unrestricted adjoints. See SELF-CENSORING.md for commands.
The latest simulation archives pin this numerical implementation; development
changes to package documentation do not silently replace those builds.
