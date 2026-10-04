cf_sporadic <- function(task, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_sporadic(task, fold, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(weights = recombine(rbind_depth(ans, "weights"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_sporadic <- function(task, fold, learners, control, pb) {
  # Create folds from the "natural" data
  natural <- get_folded_data(task$natural, task$folds, fold)

  # Initialize matrix to store predictions, one column for each time point initially set to 1
  weights <- matrix(data = 1, nrow = nrow(natural$valid), ncol = task$tau)
  # List to store model fits
  fits <- vector("list", length = task$tau)

  for (t in 1:task$tau) {
    # Indicators for censoring
    i <- task$observed(natural$train, t - 1)

    # Figure out the treatment at time t
    if (length(task$vars$A) > 1) {
      # if there are multiple treatments, we need to use the treatment at time t
      A_t <- task$vars$A[[t]]
    } else {
      # if there is only one treatment, we can use the treatment at time 1
      A_t <- task$vars$A[[1]]
    }

    # Create an indicator for sporadic missingness
    # If there is no time-varying outcome, there can't be sporadic outcome measurement
    ..i..R <- sporadic_missing(task, natural$train, t)

    # Adding sporadic measurement indicator to training data
    natural$train$..i..R <- ..i..R

    # Variables needed for estimation
    vars <- c("..i..lmtp_id", task$vars$history("A", t), A_t, "..i..R")

    # Simple test that there is variation in the sporadic measurement indicator
    if (length(unique(..i..R[i])) == 1L) {
      fit <- stats::lm(..i..R ~ 1, data = natural$train[i, vars])
    } else {
      # Fit model for sporadic missingness using the super learner
      fit <- run_ensemble(natural$train[i, vars], "..i..R",
                          learners, "binomial", "..i..lmtp_id",
                          learner_split(task, fold, natural$train$..i..lmtp_id[i], control$.learners_trt_folds),
                          control$.discrete,
                          control$.info)
    }

    # Save fit
    if (inherits(fit, "lmtp_ensemble")) {
      if (control$.return_full_fits) {
        fits[[t]] <- fit
      } else {
        fits[[t]] <- extract_sl_weights(fit)
      }
    }

    # Indicators for censoring in the validation set
    i <- task$observed(natural$valid, t - 1)

    # Create an indicator for sporadic measurement in the validation set
    # Using same logic as for training data
    ..i..R <- sporadic_missing(task, natural$valid, t)

    # Predict on validation set
    pred <- matrix(1, nrow = nrow(natural$valid), ncol = 1)
    pred[i, ] <- pmax(1 - predict(fit, natural$valid[i, ]), 1e-5)

    # Create IPW weights
    weights[, t] <- (!..i..R) / pred

    pb()
  }

  list(weights = weights, fits = fits)
}

sporadic_missing <- function(task, data, t) {
  outcomes <- c(task$vars$N, task$vars$Y)
  if (length(outcomes) == 1L && t < task$tau) return(rep(0, nrow(data)))
  outcome <- outcomes[min(t, length(outcomes))]
  as.numeric(task$observed(data, t) & is.na(data[[outcome]]))
}
