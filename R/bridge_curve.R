#' Observed-history counterfactual curves under informative measurement
#'
#' Fits fold-specific bridges and adjoints, calls the pooled SDR or TMLE
#' sequential estimator, and adds the adjoint correction as a one-step update.
#' Measurement is not censoring. Recorded outcomes remain history variables;
#' transformed terminal responses are substituted only inside each outer fold.
#'
#' @param data A data.frame containing observed data.
#' @param trt Treatment columns, or one vector of columns per time.
#' @param outcome Outcome columns, one per horizon.
#' @param measurement Binary measurement columns, one per outcome horizon.
#' @param health Columns in the measured health vector, one vector per horizon.
#' @param baseline Baseline predictor columns.
#' @param time_vary Recorded current covariate columns, one vector per treatment time.
#' @param shift Modified treatment policy.
#' @param estimator Either "sdr" or "tmle".
#' @param folds Outer split count or an explicit split list.
#' @param learner_folds Inner labels in outer training-row order.
#' @param learners_outcome,learners_trt SuperLearner libraries.
#' @param bridge_library,adjoint_library Conditional-moment libraries. The
#'   default adjoint copies bridge learners with lower and upper bounds removed.
#' @param bounds Fixed bounds for the transformed outcome, including zero.
#' @param control lmtp controls.
#' @param id Optional cluster ID column.
#' @param mtp Whether the policy is a modified rather than deterministic policy.
#' @param nuisance_fit Optional callback (training_data, validation_data,
#'   inner_folds, task, fold). Returns beta_training, beta_validation,
#'   lambda_validation matrices; optionally ratios_training, ratios_validation
#'   and oracle_regressions for controlled validation. All values use raw scales.
#' @return An lmtp_curve object with the combined gradient and bridge diagnostics.
#' @export
lmtp_bridge <- function(data, trt, outcome, measurement, health,
                        baseline = NULL, time_vary = NULL, shift,
                        estimator = c("sdr", "tmle"), folds = 3,
                        learner_folds = NULL,
                        learners_outcome = c("SL.saturated_l1_small", "SL.saturated_l1", "SL.saturated_l1_large"),
                        learners_trt = learners_outcome,
                        bridge_library = NULL, adjoint_library = NULL,
                        bounds = c(0, 6), control = lmtp_control(), id = NULL,
                        mtp = TRUE, nuisance_fit = NULL) {
  estimator <- match.arg(estimator)
  if (length(measurement) != length(outcome) || length(health) != length(outcome)) {
    stop("measurement and health must have one entry per outcome horizon.")
  }
  if (bounds[1] != 0) stop("Bridge outcome bounds must start at zero.")
  original <- as.data.frame(data)
  if (!all(unlist(original[measurement]) %in% c(0, 1))) stop("Measurement indicators must be binary and observed.")
  for (s in seq_along(outcome)) {
    if (!outcome[s] %in% health[[s]]) stop("Each health vector must include its outcome.")
    measured <- original[[measurement[s]]] == 1
    if (any(!is.finite(as.matrix(original[measured, health[[s]], drop = FALSE])))) {
      stop("Measured health must be numeric, finite and observed.")
    }
    # An unmeasured value cannot enter recorded history, even if a simulation
    # passes a complete latent-data frame alongside measurement indicators.
    original[!measured, health[[s]]] <- NA_real_
  }
  encoded <- original
  encoded$..bridge_row <- seq_len(nrow(encoded))
  original$..bridge_row <- seq_len(nrow(original))
  if (!is.null(time_vary)) for (s in seq_along(time_vary)) {
    columns <- character()
    for (name in time_vary[[s]]) {
      value <- original[[name]]
      if (!is.numeric(value)) stop("Time-varying bridge covariates must have numeric finite-level coding.")
      copy <- paste0("..bridge_", s, "_", name)
      missing <- paste0(copy, "_missing")
      encoded[[copy]] <- ifelse(is.na(value), 0, value)
      encoded[[missing]] <- as.numeric(is.na(value))
      columns <- c(columns, copy, missing)
    }
    time_vary[[s]] <- columns
  }
  task <- LmtpTask$new(data = encoded, shifted = make_shifted(encoded, trt, NULL, shift, NULL),
                       A = trt, Y = outcome, L = time_vary, W = baseline,
                       C = NULL, D = NULL, k = Inf, id = id, outcome_type = "continuous",
                       folds = folds, weights = NULL, bounds = bounds, learner_folds = learner_folds)
  prepare_curve_splits(task, control)
  if (is.null(bridge_library)) {
    bridge_library <- list(
      small = list(method = "saturated_l1", control = list(penalty = .001 / sqrt(task$n), lower = 1, upper = bounds[2])),
      medium = list(method = "saturated_l1", control = list(penalty = .01 / sqrt(task$n), lower = 1, upper = bounds[2])),
      large = list(method = "saturated_l1", control = list(penalty = .1 / sqrt(task$n), lower = 1, upper = bounds[2])))
  }
  if (is.null(adjoint_library)) {
    adjoint_library <- lapply(bridge_library, function(x) { x$control$lower <- -Inf; x$control$upper <- Inf; x })
  }
  fitted <- vector("list", length(task$folds))
  ratio_matrix <- correction <- matrix(NA_real_, task$n, task$tau)
  for (j in seq_along(task$folds)) {
    train <- task$folds[[j]]$training_set
    valid <- task$folds[[j]]$validation_set
    got <- if (is.null(nuisance_fit)) {
      bridge_nuisance_fit(original, encoded, task, j, measurement, health,
                           bridge_library, adjoint_library, learners_trt, mtp, control)
    } else nuisance_fit(original[train, , drop = FALSE], original[valid, , drop = FALSE],
                         task$learner_folds[[j]], task, j)
    for (name in c("beta_training", "beta_validation", "lambda_validation")) {
      nr <- if (name == "beta_training") length(train) else length(valid)
      if (!is.matrix(got[[name]]) || !identical(dim(got[[name]]), c(nr, task$tau))) stop("Invalid nuisance matrix: ", name)
    }
    terminal <- function(rows, beta) {
      ans <- matrix(0, length(rows), task$tau)
      for (s in seq_len(task$tau)) {
        measured <- original[[measurement[s]]][rows] == 1
        ans[measured, s] <- original[[outcome[s]]][rows[measured]] * beta[measured, s]
      }
      if (any(!is.finite(ans)) || any(ans < bounds[1] | ans > bounds[2])) stop("Bridge outcomes outside supplied bounds.")
      ans
    }
    got$terminal_training <- terminal(train, got$beta_training)
    got$terminal_validation <- terminal(valid, got$beta_validation)
    for (s in seq_len(task$tau)) {
      m <- original[[measurement[s]]][valid]
      residual <- rep(-1, length(valid))
      residual[m == 1] <- got$beta_validation[m == 1, s] - 1
      correction[valid, s] <- -got$lambda_validation[, s] * residual
    }
    if (!is.null(got$ratios_validation)) ratio_matrix[valid, ] <- got$ratios_validation
    fitted[[j]] <- got
  }
  ratios <- if (all(is.finite(ratio_matrix))) list(ratios = ratio_matrix, fits = NULL) else
    cf_r(task, learners_trt, mtp, control, function(...) NULL)
  sporadic <- matrix(1, task$n, task$tau)
  curve <- cf_curve(task, ratios$ratios, sporadic, learners_outcome, control,
                     function(...) NULL, estimator, learners_trt, mtp, fitted)
  sequential <- theta_curve(task, curve$influence_functions, curve$shifted,
                             ratios$ratios, sporadic, curve$fits, ratios$fits, NULL,
                             deparse(substitute(shift)), estimator, curve$targeting)
  out <- sequential
  out$sequential_estimates <- sequential$estimates
  out$bridge_correction <- correction
  out$bridge_fits <- fitted
  out$measurement <- measurement
  for (s in seq_len(task$tau)) {
    point <- sequential$estimates[[s]]
    out$estimates[[s]] <- ife::ife(point@x + mean(correction[, s]), point@eif + correction[, s],
                                  task$weights, as.character(task$id))
  }
  out$unprojected_estimates <- vapply(out$estimates, function(x) x@x, numeric(1))
  out$estimator <- paste(toupper(estimator), "observed-history bridge curve")
  out
}

bridge_nuisance_fit <- function(original, encoded, task, fold, measurement, health,
                                bridge_library, adjoint_library, learners, mtp, control) {
  rows <- task$folds[[fold]]$training_set
  valid <- task$folds[[fold]]$validation_set
  labels <- task$learner_folds[[fold]]
  ratio_cache <- new.env(parent = emptyenv())
  split_cache <- new.env(parent = emptyenv())
  nested_labels <- function(local_train) {
    key <- digest::digest(local_train)
    if (!exists(key, split_cache, inherits = FALSE)) {
      if (identical(local_train, seq_along(rows))) value <- labels else {
        plan <- make_folds(encoded[rows[local_train], , drop = FALSE], V = 2,
                           cluster_ids = task$id[rows[local_train]])
        value <- integer(length(local_train))
        for (k in seq_along(plan)) value[plan[[k]]$validation_set] <- k
      }
      assign(key, value, split_cache)
    }
    get(key, split_cache, inherits = FALSE)
  }
  ratio_fit <- function(local_train, local_valid) {
    # Hash row indices: environment binding names have a 10,000-byte limit.
    key <- digest::digest(local_train)
    if (!exists(key, ratio_cache, inherits = FALSE)) {
      sub <- task$clone(deep = TRUE)
      sub$folds <- list(list(training_set = rows[local_train], validation_set = c(rows, valid)))
      sub$learner_folds <- list(nested_labels(local_train))
      value <- estimate_r(sub, 1, learners, mtp, control, function(...) NULL)$ratios
      assign(key, value, ratio_cache)
    }
    get(key, ratio_cache, inherits = FALSE)
  }
  final_ratios <- ratio_fit(seq_along(rows), integer())
  training_ratios <- matrix(NA_real_, length(rows), task$tau)
  for (label in unique(labels)) {
    iv <- which(labels == label)
    training_ratios[iv, ] <- ratio_fit(which(labels != label), iv)[iv, , drop = FALSE]
  }
  oof_cache <- new.env(parent = emptyenv())
  training_ratio_fit <- function(train) {
    if (identical(train, seq_along(rows))) return(training_ratios)
    key <- digest::digest(train)
    if (!exists(key, oof_cache, inherits = FALSE)) {
      inner <- nested_labels(train)
      value <- matrix(NA_real_, length(train), task$tau)
      for (label in unique(inner)) {
        vi <- which(inner == label)
        fit <- ratio_fit(train[inner != label], train[vi])
        value[vi, ] <- fit[train[vi], , drop = FALSE]
      }
      assign(key, value, oof_cache)
    }
    get(key, oof_cache, inherits = FALSE)
  }
  bt <- matrix(1, length(rows), task$tau)
  bv <- lv <- matrix(1, length(valid), task$tau)
  beta_fits <- adjoint_fits <- vector("list", task$tau)
  for (s in seq_len(task$tau)) {
    history <- task$vars$history("A", s)
    A <- if (length(task$vars$A) > 1L) task$vars$A[[s]] else task$vars$A[[1]]
    B <- as.matrix(encoded[c(rows, valid), c(history, A), drop = FALSE])
    V <- as.matrix(cbind(encoded[c(rows, valid), history, drop = FALSE],
                         original[c(rows, valid), health[[s]], drop = FALSE]))
    M <- original[[measurement[s]]][c(rows, valid)]
    ntrain <- length(rows)
    beta <- cmbridge::fit_bridge_ensemble(B[seq_len(ntrain), , drop = FALSE], V[seq_len(ntrain), , drop = FALSE],
                                          M[seq_len(ntrain)], library = bridge_library,
                                          fold_id = labels, kernel = "cell")
    beta_fits[[s]] <- beta
    observed <- which(M[ntrain + seq_along(valid)] == 1)
    bv[observed, s] <- predict(beta, V[ntrain + observed, , drop = FALSE])
    # Inner out-of-fold bridge responses for sequential-regression training.
    for (label in unique(labels)) {
      ti <- which(labels != label); vi <- which(labels == label & M[seq_len(ntrain)] == 1)
      nested <- nested_labels(ti)
      inner_beta <- cmbridge::fit_bridge_ensemble(B[ti, , drop = FALSE], V[ti, , drop = FALSE], M[ti],
                                                  library = bridge_library, fold_id = nested, kernel = "cell")
      bt[vi, s] <- predict(inner_beta, V[vi, , drop = FALSE])
    }
    loading <- function(train, validation) {
      ratios <- ratio_fit(train, validation)
      omega <- apply(ratios[seq_len(ntrain), seq_len(s), drop = FALSE], 1, prod)
      y <- original[[task$vars$Y[s]]][rows]
      y[is.na(y)] <- 0
      value <- omega * y
      # The loading response on training rows also uses out-of-fold ratios.
      oof <- training_ratio_fit(train)
      train_value <- apply(oof[, seq_len(s), drop = FALSE], 1, prod) * y[train]
      list(train = train_value, validation = value[validation])
    }
    adjoint <- cmbridge::fit_adjoint_ensemble(B[seq_len(ntrain), , drop = FALSE], V[seq_len(ntrain), , drop = FALSE],
                                             M[seq_len(ntrain)], loading, library = adjoint_library,
                                             fold_id = labels, kernel = "cell")
    adjoint_fits[[s]] <- adjoint
    lv[, s] <- predict(adjoint, B[ntrain + seq_along(valid), , drop = FALSE])
  }
  list(beta_training = bt, beta_validation = bv, lambda_validation = lv,
       ratios_training = training_ratios,
       ratios_validation = final_ratios[length(rows) + seq_along(valid), , drop = FALSE],
       beta_fits = beta_fits, adjoint_fits = adjoint_fits)
}
