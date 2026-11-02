#' Full-history counterfactual curves under informative measurement
#'
#' Implements the full-history identification and gradient recursion, including
#' a separate bridge at each transition and target-specific adjoint corrections.
#' Regressions are pooled by backward depth. Measurement is intermittent, not
#' censoring; only complete retained states enter conditional fits.
#'
#' @param data Observed wide data.frame.
#' @param trt Treatment column names, or one character vector per time.
#' @param outcome Outcome columns at times 2 through tau+1.
#' @param measurement Binary measurement columns at times 1 through tau+1.
#'   Baseline measurement must equal one.
#' @param health List of health-column vectors at times 1 through tau+1.
#'   Corresponding components must have the same interpretation across times.
#'   Each follow-up health vector must include its outcome.
#' @param alpha Known no-visit treatment rule, function(data, trt), or a list
#'   of such functions, one per treatment time. It receives pre-treatment
#'   history with the current treatment columns removed.
#' @param d Visit-time policy, function(data, trt), or a list of such functions.
#' @param t No-visit policy, function(data, trt, epsilon), or a list. The
#'   treatment columns supplied to this function already contain alpha(data,trt).
#'   Each policy returns a numeric vector for scalar treatment, or a matrix or
#'   data.frame with one column per treatment component.
#' @param epsilon Equal-weight integration nodes for the known randomizer law,
#'   or a list of nodes per time. Defaults to 20 uniform midpoints. For other
#'   laws supply matching quantiles or draws. Increase nodes to check integration
#'   accuracy. These nodes integrate regression predictions only.
#' @param randomizer Function(n, time) drawing n independent randomizers from
#'   the known law Q at a treatment time. Used to generate continuation-policy
#'   classifier samples. Defaults to independent Uniform(0,1) draws. This law
#'   must match epsilon; continuous policies require a continuous randomizer.
#' @param ratio_draws Positive number of independent continuation samples per
#'   training person. The classifier odds are divided by this class-size ratio.
#' @param mtp Whether d is a modified policy (TRUE) or a static policy (FALSE).
#'   Static policies assign exact zero density ratios to non-followers.
#' @param baseline Always-observed baseline columns retained at every time.
#' @param k Retained health-history length, a positive integer or Inf. Alpha
#'   must be determined by this state and the retained observed history.
#' @param estimator Either "sdr" or "tmle" (weighted logistic fluctuation).
#' @param folds Outer fold count (at least two), or explicit person splits.
#' @param learner_folds Optional inner labels in outer training-row order.
#' @param learners_outcome,learners_trt SuperLearner libraries.
#' @param bridge_library,adjoint_library Conditional-moment libraries, as in
#'   [lmtp_bridge()]. Positive-penalty tuning uses only training observations.
#' @param bounds Fixed interval containing all recursive pseudo-outcomes,
#'   including zero. Values outside it cause an error, never silent clipping.
#' @param control Controls from [lmtp_control()]. Density-ratio trimming is
#'   calibrated on training predictions only; use .trim=1 to disable it.
#' @param id Optional cluster ID column.
#' @param nuisance_cache Optional environment shared by SDR/TMLE calls with the
#'   same data, policies and splits. Only identical training problems are reused;
#'   altered training data or tuning specifications create separate cache entries.
#' @return An lmtp_curve with estimates, combined influence functions, separate
#'   sequential and bridge components, propagation weights and nuisance fits.
#' @export
lmtp_full_history <- function(data, trt, outcome, measurement, health,
                              alpha, d, t, epsilon = (seq_len(20) - .5) / 20,
                              randomizer = function(n, time) stats::runif(n),
                              ratio_draws = 5L, mtp = TRUE,
                              baseline = NULL, k = Inf,
                              estimator = c("sdr", "tmle"), folds = 3,
                              learner_folds = NULL,
                              learners_outcome = c("SL.saturated_l1_cv", "SL.mean"),
                              learners_trt = learners_outcome,
                              bridge_library = NULL, adjoint_library = NULL,
                              bounds = c(-100, 100), control = lmtp_control(), id = NULL,
                              nuisance_cache = NULL) {
  estimator <- match.arg(estimator)
  if (is.null(nuisance_cache)) nuisance_cache <- new.env(parent = emptyenv())
  if (!is.environment(nuisance_cache)) stop("nuisance_cache must be an environment.")
  z <- fh_setup(data, trt, outcome, measurement, health, baseline, k, alpha, d, t, epsilon, id)
  if (!is.function(randomizer) || length(ratio_draws) != 1L || !is.finite(ratio_draws) ||
      ratio_draws < 1L || ratio_draws != as.integer(ratio_draws))
    stop("randomizer must be a function and ratio_draws a positive integer.")
  if (length(mtp) != 1L || !is.logical(mtp) || is.na(mtp)) stop("mtp must be TRUE or FALSE.")
  z$randomizer <- randomizer; z$ratio_draws <- as.integer(ratio_draws); z$mtp <- mtp
  if (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1] >= bounds[2] ||
      bounds[1] > 0 || bounds[2] < 0) stop("bounds must be a finite increasing interval containing zero.")
  if (!is.list(folds)) {
    if (length(folds) != 1L || !is.finite(folds) || folds < 2 || folds != as.integer(folds))
      stop("Full-history estimation requires at least two outer folds.")
    folds <- make_folds(z$data, V = folds, cluster_ids = z$id)
  }
  if (length(folds) < 2L) stop("Full-history estimation requires at least two outer folds.")
  folds <- validate_folds(folds, z$n, z$id)
  learner_folds <- validate_learner_folds(learner_folds, folds, z$id)
  libraries <- fh_libraries(bridge_library, adjoint_library)
  seq_score <- br_score <- omega <- matrix(0, z$n, z$tau)
  ratios_d <- ratios_t <- bridge_values <- matrix(NA_real_, z$n, z$tau)
  fits <- targeting <- vector("list", length(folds))
  for (j in seq_along(folds)) {
    train <- folds[[j]]$training_set; valid <- folds[[j]]$validation_set
    labels <- if (is.null(learner_folds)) NULL else learner_folds[[j]]
    engine <- fh_engine(z, train, labels, learners_outcome, learners_trt,
                        libraries, bounds, control, estimator, nuisance_cache)
    nu <- engine$nuisance(train, valid)
    omega[valid, ] <- nu$omega
    ratios_d[valid, ] <- nu$d; ratios_t[valid, ] <- nu$t; bridge_values[valid, ] <- nu$beta
    models <- lapply(0:(z$tau - 1L), function(depth) engine$model(train, depth))
    predictions <- lapply(seq_len(z$tau), function(h) lapply(seq_len(h), function(s)
      engine$predictions(models[[h - s + 1L]], valid, s)))
    target_j <- vector("list", z$tau)
    adjoints <- vector("list", z$tau)
    for (h in seq_len(z$tau)) {
      p <- predictions[[h]]
      eta <- engine$terminal(train, valid, h)
      diagnostics <- list()
      for (s in rev(seq_len(h))) {
        if (estimator == "tmle") {
          ii <- which(z$J[valid, s] & nu$omega[, s] > 0)
          if (!length(ii)) stop("No positive-weight validation observations for full-history TMLE.")
          scaled <- (eta[ii] - bounds[1]) / diff(bounds)
          fh_check_bounds(eta[ii], bounds)
          eps <- curve_fluctuation(scaled, (p[[s]]$natural[ii] - bounds[1]) / diff(bounds),
                                   nu$omega[ii, s], control$.bound)
          for (name in c("natural", "d", "t")) {
            # Integrate after fluctuation: logit updates do not commute with
            # averaging over the continuation randomizer.
            if (name == "t") {
              values <- p[[s]]$t_nodes
              values <- bounds[1] + diff(bounds) * plogis(qlogis(bound(
                (values - bounds[1]) / diff(bounds), control$.bound)) + eps)
              p[[s]]$t <- rowMeans(values)
            } else p[[s]][[name]] <- bounds[1] + diff(bounds) * plogis(qlogis(bound(
              (p[[s]][[name]] - bounds[1]) / diff(bounds), control$.bound)) + eps)
          }
          diagnostics[[length(diagnostics) + 1L]] <- data.frame(time = s, epsilon = eps,
            score = weighted.mean(eta[ii] - p[[s]]$natural[ii], nu$omega[ii, s]))
        }
        if (s > 1L) eta <- fh_q(z, valid, s, p[[s]], nu$beta[, s - 1L])
      }
      score <- p[[1]]$d
      for (s in seq_len(h)) {
        eta <- if (s == h) engine$terminal(train, valid, h) else
          fh_q(z, valid, s + 1L, p[[s + 1L]], nu$beta[, s])
        score <- score + nu$omega[, s] * (eta - p[[s]]$natural)
        adj <- engine$adjoint(train, valid, s, h)
        residual <- rep(-1, length(valid))
        measured <- which(z$data[[z$measurement[s + 1L]]][valid] == 1 & z$J[valid, s])
        residual[measured] <- nu$beta[measured, s] - 1
        br_score[valid, h] <- br_score[valid, h] - z$J[valid, s] * adj$value * residual
        adjoints[[h]][[s]] <- adj$model
      }
      seq_score[valid, h] <- if (estimator == "tmle") p[[1]]$d else score
      # The empirical gradient retains the sequential residual for TMLE too.
      predictions[[h]] <- list(predictions = p, score = score)
      target_j[[h]] <- if (length(diagnostics)) do.call(rbind, diagnostics) else NULL
    }
    fits[[j]] <- list(bridges = nu$fits, regressions = models, adjoints = adjoints,
                      beta_validation = nu$beta, ratios_d = nu$d, ratios_t = nu$t,
                      training_rows = train, validation_rows = valid, learner_folds = engine$splits(train),
                      sequential = predictions, cache_counts = engine$counts())
    targeting[[j]] <- target_j
  }
  points <- colMeans(seq_score + br_score)
  gradients <- seq_score
  if (estimator == "tmle") for (j in seq_along(folds)) for (h in seq_len(z$tau))
    gradients[folds[[j]]$validation_set, h] <- fits[[j]]$sequential[[h]]$score
  gradients <- gradients + br_score
  out <- list(estimator = paste(toupper(estimator), "full-history bridge curve"),
    estimates = lapply(seq_len(z$tau), function(h) ife::ife(points[h], gradients[, h],
      rep(1, z$n), as.character(z$id))),
    folds = folds, learner_folds = lapply(fits, function(x) x$learner_folds),
    shift = paste(deparse(substitute(d)), collapse = ""), outcome_type = "continuous",
    unprojected_estimates = points, sequential_estimates = colMeans(seq_score),
    sequential_score = gradients - br_score, bridge_correction = br_score,
    influence_functions = sweep(gradients, 2L, points), propagation_weights = omega,
    density_ratios_d = ratios_d, density_ratios_t = ratios_t, bridge_values = bridge_values,
    bridge_fits = fits, targeting = targeting, measurement = measurement, k = k,
    integration_nodes = z$epsilon)
  class(out) <- "lmtp_curve"
  out
}

fh_libraries <- function(bridge, adjoint) {
  if (is.null(bridge)) bridge <- list(
    sieve_md = list(method = "sieve_md", control = list(target_basis = "cell",
      instrument_basis = "cell", link = "inverse_logit", penalty_scales = 10^seq(-10, 0, by = .25))),
    landweber = list(method = "landweber", control = list(link = "inverse_logit")),
    pmmr = list(method = "pmmr", control = list(link = "inverse_logit")))
  if (is.null(adjoint)) adjoint <- lapply(bridge, function(candidate) {
    if (is.character(candidate)) candidate <- list(method = candidate)
    if (is.null(candidate$control)) candidate$control <- list()
    candidate$control$link <- "identity"
    candidate
  })
  list(bridge = bridge, adjoint = adjoint)
}

fh_setup <- function(data, trt, outcome, measurement, health, baseline, k, alpha, d, t, epsilon, id) {
  data <- as.data.frame(data); n <- nrow(data); tau <- length(outcome)
  if (!tau || n < 4L) stop("At least one outcome and four observations are required.")
  if (!is.list(trt)) trt <- as.list(trt)
  if (length(trt) != tau || length(health) != tau + 1L || length(measurement) != tau + 1L)
    stop("trt/outcome need tau entries; health/measurement need tau+1 including baseline.")
  if (any(lengths(health) != length(health[[1L]])) || !length(health[[1L]]) ||
      any(lengths(trt) != length(trt[[1L]])) || !length(trt[[1L]]))
    stop("Health and treatment vectors must have consistent positive dimensions across times.")
  if (length(k) != 1L || is.na(k) || k < 1 || (!is.infinite(k) && k != as.integer(k)))
    stop("k must be a positive integer or Inf.")
  columns <- unique(c(unlist(trt), outcome, measurement, unlist(health), baseline, id))
  if (!all(columns %in% names(data))) stop("Some requested columns are absent from data.")
  if (anyDuplicated(c(unlist(trt), unlist(health), measurement)))
    stop("Treatment, health and measurement columns must be distinct across times.")
  if (!all(unlist(data[measurement]) %in% c(0, 1))) stop("Measurement must be observed and binary.")
  if (any(data[[measurement[1L]]] != 1)) stop("Baseline health must be measured for everyone.")
  for (s in seq_len(tau + 1L)) {
    observed <- data[[measurement[s]]] == 1
    if (!all(vapply(data[health[[s]]], is.numeric, TRUE)) ||
        any(!is.finite(as.matrix(data[observed, health[[s]], drop = FALSE]))))
      stop("Measured health must be numeric and finite.")
    if (s > 1L && !outcome[s - 1L] %in% health[[s]]) stop("Each outcome must belong to its health vector.")
    data[!observed, health[[s]]] <- NA_real_
  }
  if (!all(vapply(data[c(unlist(trt), baseline)], is.numeric, TRUE)) ||
      any(!is.finite(as.matrix(data[c(unlist(trt), baseline)]))))
    stop("Treatment and baseline predictors must be observed, numeric and finite.")
  ids <- if (is.null(id)) seq_len(n) else data[[id]]
  if (anyNA(ids)) stop("Cluster IDs cannot be missing.")
  rules <- lapply(list(alpha = alpha, d = d, t = t), function(rule) {
    if (is.function(rule)) rule <- rep(list(rule), tau)
    if (length(rule) != tau || !all(vapply(rule, is.function, TRUE)))
      stop("alpha, d and t must be functions or lists of tau functions.")
    rule
  })
  if (!is.list(epsilon)) epsilon <- rep(list(epsilon), tau)
  if (length(epsilon) != tau || any(!lengths(epsilon)) ||
      !all(vapply(epsilon, function(x) is.numeric(x) && all(is.finite(x)), TRUE)))
    stop("epsilon must supply nonempty finite numeric integration nodes for every time.")
  J <- matrix(FALSE, n, tau)
  for (s in seq_len(tau)) {
    first <- if (is.infinite(k)) 1L else max(1L, s - k + 1L)
    J[, s] <- rowSums(data[measurement[first:s]]) == s - first + 1L
  }
  z <- list(data = data, n = n, tau = tau, trt = trt, outcome = outcome,
    measurement = measurement, health = health, baseline = baseline, k = k,
    id = ids, epsilon = epsilon, rules = rules, J = J)
  # Verify the known no-visit rule without exposing latent health or outcomes.
  for (s in seq_len(tau)) {
    rows <- which(data[[measurement[s]]] == 0)
    if (length(rows)) {
      value <- fh_policy(z, rows, s, "alpha")
      if (any(abs(value - as.matrix(data[rows, trt[[s]], drop = FALSE])) > 1e-8))
        stop("Observed no-visit treatment disagrees with alpha at time ", s, ".")
    }
  }
  z
}

# Policy callbacks receive only the retained observed history, never future
# outcomes or unrecorded covariates accidentally supplied by a simulator.
fh_policy <- function(z, rows, s, branch, epsilon = NULL) {
  first <- if (is.infinite(z$k)) 1L else max(1L, s - z$k + 1L)
  columns <- unique(c(z$baseline, unlist(z$health[first:s]),
    unlist(z$trt[first:s]), z$measurement[first:s]))
  data <- z$data[rows, columns, drop = FALSE]; a <- z$trt[[s]]
  if (branch == "alpha") data <- data[, setdiff(names(data), a), drop = FALSE]
  if (branch == "t") data[, a] <- fh_policy(z, rows, s, "alpha")
  args <- list(data, a)
  if (branch == "t") args$epsilon <- epsilon
  value <- do.call(z$rules[[branch]][[s]], args)
  if (is.null(dim(value))) {
    if (length(a) != 1L || length(value) != length(rows)) stop("Invalid policy output dimensions.")
    value <- matrix(value, ncol = 1L)
  }
  value <- as.matrix(value)
  if (!identical(dim(value), c(length(rows), length(a))) || any(!is.finite(value)))
    stop("Policy output must be finite with one row per observation and one column per treatment.")
  value
}

# Canonical lag names allow pooling across time. Padding carries an explicit
# indicator and cannot be confused with an observed zero health value.
fh_design <- function(z, rows, s, treatment = NULL, bridge = FALSE) {
  depth <- if (is.infinite(z$k)) z$tau else min(z$k, z$tau)
  x <- z$data[rows, z$baseline, drop = FALSE]
  for (lag in 0:(depth - 1L)) {
    at <- s - lag
    for (j in seq_along(z$health[[1L]]))
      x[[paste0("health", j, "_lag", lag)]] <- if (at < 1L) 0 else z$data[[z$health[[at]][j]]][rows]
    x[[paste0("available_lag", lag)]] <- as.numeric(at >= 1L)
    if (lag > 0L) for (j in seq_along(z$trt[[1L]]))
      x[[paste0("trt", j, "_lag", lag)]] <- if (at < 1L) 0 else z$data[[z$trt[[at]][j]]][rows]
  }
  if (bridge) {
    for (j in seq_along(z$health[[1L]]))
      x[[paste0("next_health", j)]] <- z$data[[z$health[[s + 1L]][j]]][rows]
  } else {
    if (is.null(treatment)) treatment <- as.matrix(z$data[rows, z$trt[[s]], drop = FALSE])
    for (j in seq_along(z$trt[[s]])) x[[paste0("trt", j)]] <- treatment[, j]
  }
  x$time <- rep(s, length(rows))
  x
}

fh_check_bounds <- function(y, bounds) {
  if (any(!is.finite(y)) || any(y < bounds[1] | y > bounds[2]))
    stop("Recursive pseudo-outcomes exceed bounds; supply a wider fixed interval.")
  invisible(y)
}

fh_q <- function(z, rows, s, p, beta) {
  value <- numeric(length(rows))
  observed <- which(z$data[[z$measurement[s]]][rows] == 1 & z$J[rows, s])
  value[observed] <- p$d[observed] + (beta[observed] - 1) * p$t[observed]
  value
}
