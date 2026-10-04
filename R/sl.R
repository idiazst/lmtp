#' @importFrom nnls nnls
run_ensemble <- function(data, y, learners, outcome_type, id, folds, discrete, info) {
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
  fit <- SuperLearner::SuperLearner(
    Y = data[[y]], X = data[, features, drop = FALSE],
    family = if (outcome_type == "binomial") binomial() else gaussian(),
    SL.library = learners, id = data[[id]], method = "method.NNLS",
    env = environment(SuperLearner::SuperLearner), cvControl = cv_control,
    verbose = info
  )
  if (discrete) {
    fit$coef[] <- 0
    fit$coef[which.min(fit$cvRisk)] <- 1
  }
  class(fit) <- append("lmtp_ensemble", class(fit))
  fit
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
