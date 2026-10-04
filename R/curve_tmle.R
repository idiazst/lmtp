# Pool the initial regressions by recursive depth, then target each outcome
# time separately. A single pooled intercept would not solve every curve point's
# score equation. No targeted validation responses enter a nuisance training fit.
target_curve <- function(task, fold, ratios, sporadic_weights, natural, shifted, tol) {
  data <- task$natural[task$folds[[fold]]$validation_set, , drop = FALSE]
  sampling_weights <- task$weights[task$folds[[fold]]$validation_set]
  targeting <- vector("list", task$tau)
  influence_functions <- matrix(NA_real_, nrow(data), task$tau)

  for (s in seq_len(task$tau)) {
    targeting[[s]] <- data.frame(time = rev(seq_len(s)), epsilon = NA_real_, score = NA_real_)
    weights <- compute_weights(ratios, sporadic_weights, 1, s)
    for (r in rev(seq_len(s))) {
      at_risk <- task$at_risk_N(data, r - 1) & task$at_risk_D(data, r - 1)
      observed <- task$observed(data, r)
      y <- shifted[[s]][, r + 1L]
      q <- natural[[s]][, r]
      w <- weights[, r] * sampling_weights
      eligible <- at_risk & observed & is.finite(y) & is.finite(q) & is.finite(w) & w > 0
      if (!any(eligible)) stop("No positive-weight validation observations for TMLE targeting at outcome ", s, ", time ", r, ".")
      if (any(y[eligible] < 0 | y[eligible] > 1)) stop("TMLE targeting outcomes must be scaled to [0, 1].")

      # Same weighted logistic intercept fluctuation as upstream's held-out TMLE.
      epsilon <- curve_fluctuation(y[eligible], q[eligible], w[eligible], tol)
      update_rows <- at_risk & task$observed(data, r - 1) & is.finite(q)
      natural[[s]][update_rows, r] <- plogis(qlogis(bound(q[update_rows], tol)) + epsilon)
      shifted[[s]][update_rows, r] <- plogis(qlogis(bound(shifted[[s]][update_rows, r], tol)) + epsilon)
      i <- s - r + 1L
      targeting[[s]]$epsilon[i] <- epsilon
      targeting[[s]]$score[i] <- weighted.mean(y[eligible] - natural[[s]][eligible, r], w[eligible])
    }
    influence_functions[, s] <- eif(ratios, sporadic_weights, shifted[[s]], natural[[s]], tau = s)
  }
  list(natural = natural, shifted = shifted, influence_functions = influence_functions,
       targeting = targeting)
}

# Solve the weighted intercept score directly, avoiding GLM convergence failure
# for separated folds. +/-Inf give the appropriate boundary submodel limits.
curve_fluctuation <- function(y, q, weights, tol) {
  if (all(y == 0)) return(-Inf)
  if (all(y == 1)) return(Inf)
  offset <- qlogis(bound(q, tol))
  weights <- weights / max(weights)
  weights <- weights / sum(weights)
  score <- function(epsilon) sum(weights * (y - plogis(offset + epsilon)))
  stats::uniroot(score, c(-100, 100), tol = 1e-12)$root
}
