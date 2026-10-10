cf_curve <- function(task, ratios, sporadic_weights, learners, control, pb,
                     estimator = "sdr", learners_trt = learners, mtp = TRUE,
                     bridge_folds = NULL) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      supplied <- if (is.null(bridge_folds)) NULL else bridge_folds[[fold]]
      training_weights <- if (!is.null(supplied$ratios_training)) {
        list(ratios = supplied$ratios_training,
             sporadic = matrix(1, nrow(supplied$ratios_training), task$tau))
      } else if (estimator == "sdr") {
        curve_training_weights(task, fold, ratios, sporadic_weights, learners_trt, mtp, control,
                               bridge = !is.null(supplied))
      } else NULL
      estimate_curve(task, fold, ratios, sporadic_weights, learners, control, pb,
                     estimator, training_weights, supplied)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(influence_functions = recombine(rbind_depth(ans, "influence_functions"), task$folds),
       shifted = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "shifted", t), task$folds)),
       natural = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "natural", t), task$folds)),
       pseudo = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "pseudo_valid", t), task$folds)),
       fits = lapply(ans, \(x) x[["fits"]]),
       targeting = lapply(ans, \(x) x[["targeting"]]))
}

estimate_curve <- function(task, fold, ratios, sporadic_weights, learners, control, pb,
                           estimator, training_weights = NULL, supplied = NULL) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  ratios <- get_folded_data(ratios, task$folds, fold)
  sporadic_weights <- get_folded_data(sporadic_weights, task$folds, fold)
  if (!is.null(training_weights)) {
    ratios$train <- training_weights$ratios
    sporadic_weights$train <- training_weights$sporadic
  }
  tol <- if (estimator == "tmle") max(control$.bound, .Machine$double.eps) else NULL

  # Matrix to store predictions under the intervention for the training data
  mst <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau + 1))
  # Matrix to store predictions under the observed for the training data
  mnt <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau))
  # Matrix to store predictions under the intervention for the validation data
  msv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau + 1))
  # Matrix to store predictions under the observed for the validation data
  mnv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau))

  influence_functions <- matrix(nrow = nrow(natural$valid), ncol = task$tau)

  pseudo_valid <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau + 1))

  # Pivot into long format
  natural <- lapply(natural, \(x) pivot(x, task$vars))
  shifted <- lapply(shifted, \(x) pivot(x, task$vars))

  if (!is.null(supplied)) {
    # Replace terminal responses only AFTER pivot creates recorded-history
    # covariates. A bridge pseudo-outcome is never an observed lagged outcome.
    natural$train$..i..Y_1 <- as.vector(task$scale(supplied$terminal_training))
    natural$valid$..i..Y_1 <- as.vector(task$scale(supplied$terminal_validation))
    natural$train$..i..R <- natural$valid$..i..R <- 1
  }

  # Keep lagged variables as natural values for predictions
  # Under shift training
  ust <- natural$train
  As <- grep("^..i..A", names(ust), value = TRUE)
  As <- grep("_lag", As, value = TRUE, invert = TRUE)
  shift_train <- if (length(task$vars$A) == 1L) ust$time == 1L else rep(TRUE, nrow(ust))
  ust[shift_train, As] <- shifted$train[shift_train, As]

  # Under shift valid
  usv <- natural$valid
  shift_valid <- if (length(task$vars$A) == 1L) usv$time == 1L else rep(TRUE, nrow(usv))
  usv[shift_valid, As] <- shifted$valid[shift_valid, As]

  taus <- order(unique(natural$train$time))
  for (t in seq_along(taus)) {
    mst[[t]][, t + 1] <- subset(natural$train, time == taus[t])$..i..Y_1
    msv[[t]][, t + 1] <- subset(natural$valid, time == taus[t])$..i..Y_1

    pseudo_valid[[t]][, t + 1] <- subset(natural$valid, time == taus[t])$..i..Y_1
  }

  fits <- vector("list", length = task$tau)
  for (t in 1:task$tau) {
    if (isTRUE(task$survival)) {
      # Indicator for not experiencing the outcome or the competing risk already
      at_risk <- natural$train$..i..N == 1 & natural$train$..i..D_1 == 0
      at_risk[is.na(at_risk)] <- FALSE
    } else {
      at_risk <- rep(TRUE, nrow(natural$train))
    }
    # Indicator for having the outcome measured (i.e., no sporadic measurement)
    is_measured <- natural$train$..i..R == 1
    # Indicator for not having been censored
    is_observed <- natural$train$..i..C_1 == 1
    # Combining indicators for the training subset
    train_subset <- at_risk & is_measured & is_observed

    # Create indicators for the subset of rows to use for estimation at this time
    # This uses slightly different logic than the paper to not rely on a variable 's'
    time <- as.numeric(natural$train$time) <= rev(1:task$tau)[t]

    # Remove variables we don't need for estimation
    vars <- setdiff(names(natural$train), c("..i..C_1", "..i..C_1_lag", "..i..wide_id", "..i..N", "..i..D_1", "..i..R"))
    # Don't need the time variable for the last iteration
    if (t == task$tau) {
      vars <- setdiff(vars, "time")
    }

    # Remove variables with no variation
    for(var in setdiff(vars, c("..i..lmtp_id", "..i..Y_1", grep("^..i..A", vars, value = TRUE)))) {
      if(all(duplicated(natural$train[train_subset & time, var])[-1L])) {
        vars <- setdiff(vars, var)
      }
    }

    if (is.null(supplied$oracle_regressions)) {
    fit <- run_ensemble(
      natural$train[train_subset & time, vars],
      "..i..Y_1",
      learners,
      ifelse(t != 1, "continuous", task$outcome_type),
      "..i..lmtp_id",
      learner_split(task, fold, natural$train$..i..lmtp_id[train_subset & time], control$.learners_outcome_folds),
      control$.discrete,
      control$.info,
      if (is.null(supplied$response_engine)) NULL else local({
        depth <- t; engine <- supplied$response_engine
        function(ids, times, valid_ids, valid_times) engine$response(ids, times, valid_ids, valid_times, depth)
      }),
      if (is.null(supplied$response_engine)) NULL else supplied$response_engine$splits
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    # Update the matrices storing all predictions
    mnt <- update_m(mnt, fit, natural$train, t, task$tau, task$survival, tol)
    mst <- update_m(mst, fit, ust, t, task$tau, task$survival, tol)
    mnv <- update_m(mnv, fit, natural$valid, t, task$tau, task$survival, tol)
    msv <- update_m(msv, fit, usv, t, task$tau, task$survival, tol)
    } else {
      oracle <- supplied$oracle_regressions
      for (s in seq_len(task$tau)) {
        mnt[[s]][, seq_len(s)] <- task$scale(oracle$training$natural[[s]])
        mst[[s]][, seq_len(s)] <- task$scale(oracle$training$shifted[[s]])
        mnv[[s]][, seq_len(s)] <- task$scale(oracle$validation$natural[[s]])
        msv[[s]][, seq_len(s)] <- task$scale(oracle$validation$shifted[[s]])
      }
      fits[[t]] <- list(oracle = TRUE)
    }

    # Only When t (l in the paper) = 1 do we apply the sporadic weights
    pseudo <- list()
    for (s in t:task$tau) {
      tau <- ncol(mnt[[s]])

      pseudo_valid[[s]][, t] <- eif(ratios$valid, sporadic_weights$valid, msv[[s]], mnv[[s]], t = (s - t + 1), tau = tau)

      if (s == t) {
        # If the smallest time point in the sequence from t:tau use validation data to calculate EIF for theta at time t
        influence_functions[, t] <- eif(ratios$valid, sporadic_weights$valid, msv[[s]], mnv[[s]], t = 1, tau = tau)
      } else {
        # Otherwise calculate the DR transformation using training values as the pseudo outcome for the next iteration
        dr_transformation <- if (estimator == "sdr") {
          eif(ratios$train, sporadic_weights$train, mst[[s]], mnt[[s]], t = (s - t + 1), tau = tau)
        } else {
          # Initial TMLE nuisance recursion stays entirely within T_j.
          mst[[s]][, s - t + 1]
        }
        pseudo[[length(pseudo) + 1]] <- dr_transformation
      }
    }

    # Add the pseudo outcomes to the training data
    pseudo <- unlist(pseudo)
    natural$train[, "..i..Y_1"] <- c(pseudo, rep(NA_real_, nrow(natural$train) - length(pseudo)))

    # After the first iteration, R becomes 1 for all observations
    # Sub-setting is only based on censoring for iterations > 1
    if (t == 1) {
      natural$train$..i..R <- 1
    }

    pb()
  }

  targeting <- NULL
  if (estimator == "tmle") {
    targeted <- target_curve(task, fold, ratios$valid, sporadic_weights$valid,
                             mnv, msv, tol)
    mnv <- targeted$natural
    msv <- targeted$shifted
    influence_functions <- targeted$influence_functions
    targeting <- targeted$targeting
  }

  list(
    influence_functions = influence_functions,
    shifted = msv,
    natural = mnv,
    pseudo_valid = pseudo_valid,
    fits = fits,
    targeting = targeting,
    lmtp_id = natural$valid$..i..lmtp_id
  )
}

predict_long <- function(fit, newdata, t, tau, survival, tol = NULL) {
  # Create indicators for the subset of rows to use based on time
  time <- as.numeric(newdata$time) <= rev(1:tau)[t]
  # Empty matrix to store predictions
  predictions <- matrix(nrow = nrow(newdata[time, ]), ncol = 1)
  # Indicator for not having been censored at the previous time point
  is_observed <- newdata$..i..C_1_lag == 1

  if (isTRUE(survival)) {
    # Indicator for not experiencing the outcome already
    outcome_free <- newdata$..i..N[time] == 1
    # Indicator for not experiencing competing risk already
    competing_risk_free <- newdata$..i..D_1[time] == 0
  } else {
    outcome_free <- rep(TRUE, nrow(newdata[time, ]))
    competing_risk_free <- rep(TRUE, nrow(newdata[time, ]))
  }
  outcome_free[is.na(outcome_free)] <- FALSE
  competing_risk_free[is.na(competing_risk_free)] <- TRUE

  predict_subset <- is_observed[time] & outcome_free & competing_risk_free
  rows <- which(time)[predict_subset]
  if (length(rows)) predictions[predict_subset, 1] <- predict(fit, newdata[rows, , drop = FALSE], tol)
  predictions[!outcome_free, 1] <- 0
  predictions[!competing_risk_free, 1] <- 1
  predictions[, 1]
}

update_m <- function(m, fit, newdata, t, tau, survival, tol = NULL) {
  pred <- predict_long(fit, newdata, t, tau, survival, tol)
  time <- as.numeric(newdata$time) >= t

  for (l in seq_along(t:tau)) {
    x <- (t:tau)[l]
    j <- as.numeric(newdata$time[time]) == x
    m[[x]][, x - t + 1] <- pred[j]
  }
  m
}
