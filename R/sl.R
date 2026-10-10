#' @importFrom nnls nnls
run_ensemble <- function(data, y, learners, outcome_type, id, folds, discrete, info,
                         response_fit = NULL, penalty_folds = NULL) {
  n_ids <- length(unique(data[[id]]))
  if (is.list(folds)) {
    labels <- folds$fold_id
    valid_rows <- lapply(unique(labels), function(label) which(labels == label))
    V <- length(valid_rows)
    if (V < 2L) stop("A learner fit needs at least two nonempty inner CV folds.")
    cv_control <- SuperLearner::SuperLearner.CV.control(V = V, validRows = valid_rows)
  } else {
    if (is.null(folds)) folds <- min(10L, n_ids)
    if (folds > n_ids) {
      cli::cli_abort("The number of Super Learner folds ({folds}) exceeds the number of unique IDs ({n_ids}) in the training data. Reduce the learner fold counts in {.fn lmtp_control}.")
    }
    cv_control <- SuperLearner::SuperLearner.CV.control(V = folds)
  }
  aliases <- c(glm = "SL.glm", mean = "SL.mean", ranger = "SL.ranger", cv_glmnet = "SL.glmnet")
  mapped <- learners %in% names(aliases)
  learners[mapped] <- unname(aliases[learners[mapped]])
  features <- setdiff(names(data), c(id, y))
  learner_env <- new.env(parent = environment(SuperLearner::SuperLearner))
  for (name in c("SL.saturated_l1", "SL.saturated_l1_small", "SL.saturated_l1_large",
                 "SL.saturated_l1_cv", "SL.time_mean", "predict.SL.time_mean",
                 "predict.SL.saturated_l1")) {
    assign(name, get(name, envir = environment(run_ensemble)), envir = learner_env)
  }
  if (!is.null(response_fit)) {
    if (!is.list(folds)) stop("Response refitting requires supplied learner fold labels.")
    times <- if ("time" %in% names(data)) data$time else rep(1L, nrow(data))
    for (valid in valid_rows) {
      train <- setdiff(seq_len(nrow(data)), valid)
      rebuilt <- response_fit(data[[id]][train], times[train], data[[id]][valid], times[valid])
      data[[y]][valid] <- rebuilt$validation
    }
    for (name in unique(learners)) {
      base_learner <- get(name, envir = learner_env)
      assign(name, wrap_response_learner(base_learner, name, response_fit, penalty_folds), learner_env)
    }
  } else if (is.function(penalty_folds) && "SL.saturated_l1_cv" %in% learners) {
    base <- SL.saturated_l1_cv
    wrapped <- local({f <- base; split <- penalty_folds
      function(...) f(..., penalty_folds = split)})
    assign("SL.saturated_l1_cv", wrapped, learner_env)
  }
  boundary_failures <- list()
  for (name in unique(learners)) {
    learner <- get(name, envir = learner_env)
    guarded <- local({f <- learner
      function(...) withCallingHandlers(f(...), error = function(e) {
        if (inherits(e, "cmbridge_penalty_boundary_error"))
          boundary_failures[[length(boundary_failures) + 1L]] <<- e
      })})
    assign(name, guarded, learner_env)
  }
  fit <- SuperLearner::SuperLearner(
    Y = data[[y]], X = data[, features, drop = FALSE],
    family = if (outcome_type == "binomial") binomial() else gaussian(),
    SL.library = learners, id = data[[id]], method = "method.NNLS",
    env = learner_env, cvControl = cv_control,
    verbose = info
  )
  if (length(boundary_failures)) stop(boundary_failures[[1L]])
  if (discrete) {
    fit$coef[] <- 0
    fit$coef[which.min(fit$cvRisk)] <- 1
  }
  class(fit) <- append("lmtp_ensemble", class(fit))
  fit
}

wrap_response_learner <- function(base_learner, name, response_fit, penalty_folds) {
  force(base_learner); force(name); force(response_fit); force(penalty_folds)
  function(Y, X, newX, family, obsWeights, id, ...) {
    times <- if ("time" %in% names(X)) X$time else rep(1L, nrow(X))
    if (name != "SL.saturated_l1_cv") {
      rebuilt <- response_fit(id, times, id[FALSE], times[FALSE])
      return(base_learner(rebuilt$train, X, newX, family, obsWeights, id, ...))
    }
    selected <- .saturated_penalty_selection(Y, X, family, obsWeights, id,
      10^seq(-5, 0, by = .25), penalty_folds, response_fit)
    rebuilt <- response_fit(id, times, id[FALSE], times[FALSE])
    fit <- SL.saturated_l1(rebuilt$train, X, newX, family, obsWeights, id, penalty_scale = selected$scale)
    fit$fit$penalty_cv <- selected$penalty_cv
    fit$fit$penalty_boundary <- selected$boundary_history
    fit$fit$penalty_status <- selected$status
    fit
  }
}

#' @export
predict.lmtp_ensemble <- function(object, newdata, tol = .Machine$double.eps, ...) {
  pred <- NextMethod("predict", newdata = newdata[, object$varNames, drop = FALSE], onlySL = TRUE)
  pred <- pred$pred[, 1]
  if (is.null(tol)) {
    return(pred)
  }
  bound(pred, tol)
}
