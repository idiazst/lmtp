# Joint categorical cells supply all interactions, including time. The support
# is learned from training predictors only; an unseen cell uses the training
# mean and is recorded in the fitted object.
.saturated_keys <- function(x) {
  if (!ncol(x)) return(rep("constant", nrow(x)))
  columns <- lapply(as.data.frame(x), function(column) {
    values <- unique(column)
    as.character(values)[match(column, values)]
  })
  do.call(paste, c(columns, sep = "|"))
}
.saturated_lambda_index <- function(requested, fitted) {
  index <- vapply(requested, function(value) which.min(abs(fitted - value) / value), 0L)
  if (any(abs(fitted[index] - requested) > 1e-12 * requested))
    stop("Requested penalty was not fitted to numerical precision.")
  index
}
#' Saturated L1 SuperLearner wrapper
#' @param Y,X,newX,family,obsWeights,id Standard SuperLearner wrapper inputs.
#' @param penalty_scale Penalty multiplied by the inverse square root of IDs.
#' @param ... Additional SuperLearner arguments.
#' @return Standard SuperLearner prediction and fit list.
#' @export
SL.saturated_l1 <- function(Y, X, newX, family, obsWeights, id,
                            penalty_scale = 0.01, ...) {
  if (is.null(id)) id <- seq_along(Y)
  key <- .saturated_keys
  levels <- sort(unique(key(X)))
  if (!any(obsWeights > 0)) stop("At least one positive observation weight is required.")
  constant <- length(unique(Y[obsWeights > 0])) == 1L || length(levels) == 1L
  # Repeated finite-support rows have the same likelihood contribution as
  # their weighted cell means and counts. Keep the original ID count for the
  # penalty; aggregation is computational, not a change to the model or CV.
  cell <- match(key(X), levels)
  weights <- as.numeric(rowsum(obsWeights, cell, reorder = TRUE))
  responses <- as.numeric(rowsum(Y * obsWeights, cell, reorder = TRUE)) / weights
  responses[weights == 0] <- 0
  constant <- constant || length(unique(responses[weights > 0])) == 1L
  cell_design <- Matrix::Diagonal(length(levels))
  fit_response <- responses
  fit_weights <- weights
  if (family$family == "binomial") {
    # glmnet interprets a numeric vector as class labels. Its two-column
    # successes/failures input supports aggregated binomial counts.
    fit_response <- cbind((1 - responses) * weights, responses * weights)
    fit_weights <- rep(1, length(weights))
  }
  keep <- weights > 0
  fit_response <- if (is.matrix(fit_response)) fit_response[keep, , drop = FALSE] else fit_response[keep]
  fit <- if (constant) NULL else glmnet::glmnet(
    cbind(cell_design[keep, , drop = FALSE], 0), fit_response,
    weights = fit_weights[keep], family = family$family,
    alpha = 1, standardize = FALSE, lambda = sort(unique(c(.1, .01, .001, penalty_scale)), decreasing = TRUE) / sqrt(length(unique(id))),
    thresh = 1e-9, maxit = 1000000L)
  if (!is.null(fit) && (fit$jerr != 0 || !length(fit$lambda)))
    stop("Saturated L1 regression/classification solver failed to converge.")
  object <- list(object = fit, levels = levels, mean = weighted.mean(Y, obsWeights),
                 penalty = penalty_scale / sqrt(length(unique(id))),
                 family = family$family, predictor_names = names(X))
  # Every row in a cell has the same fitted value. Evaluate the small cell
  # design once, rather than constructing a large repeated prediction matrix.
  object$cell_predictions <- if (is.null(fit)) NULL else as.numeric(predict(fit,
    cbind(Matrix::Diagonal(length(levels)), 0), type = "response")[, .saturated_lambda_index(object$penalty, fit$lambda)])
  class(object) <- "SL.saturated_l1"
  list(pred = predict(object, newX), fit = object)
}

#' @export
predict.SL.saturated_l1 <- function(object, newdata, ...) {
  key <- .saturated_keys(newdata)
  j <- match(key, object$levels)
  out <- rep(object$mean, nrow(newdata))
  keep <- which(!is.na(j))
  if (!is.null(object$object) && length(keep)) {
    values <- object$cell_predictions
    if (is.null(values)) values <- as.numeric(predict(object$object,
      cbind(Matrix::Diagonal(length(object$levels)), 0), s = object$penalty, type = "response"))
    out[keep] <- values[j[keep]]
  }
  out
}

#' @rdname SL.saturated_l1
#' @export
SL.saturated_l1_small <- function(...) SL.saturated_l1(..., penalty_scale = 0.001)
#' @rdname SL.saturated_l1
#' @export
SL.saturated_l1_large <- function(...) SL.saturated_l1(..., penalty_scale = 0.1)

# Prepare the full joint-category design once per penalty training sample.
# All predictors, including current treatment, are used in every grid fit.
.saturated_path_setup <- function(Y, X, family, obsWeights, id) {
  levels <- sort(unique(.saturated_keys(X)))
  cell <- match(.saturated_keys(X), levels)
  weights <- as.numeric(rowsum(obsWeights, cell, reorder = TRUE))
  if (!any(weights > 0)) stop("At least one positive observation weight is required.")
  responses <- as.numeric(rowsum(Y * obsWeights, cell, reorder = TRUE)) / weights
  responses[weights == 0] <- 0
  keep <- weights > 0
  response <- if (family$family == "binomial")
    cbind((1 - responses) * weights, responses * weights) else responses
  list(levels = levels, design = cbind(Matrix::Diagonal(length(levels))[keep, , drop = FALSE], 0),
    response = if (is.matrix(response)) response[keep, , drop = FALSE] else response[keep],
    weights = if (family$family == "binomial") rep(1, sum(keep)) else weights[keep],
    constant = length(unique(Y[obsWeights > 0])) == 1L || length(levels) == 1L ||
      length(unique(responses[keep])) == 1L,
    mean = weighted.mean(Y, obsWeights), n_ids = length(unique(id)),
    family = family$family, predictor_names = names(X))
}

.saturated_path_predict <- function(setup, validation_keys, scales) {
  if (setup$constant) return(matrix(setup$mean, length(validation_keys), length(scales)))
  # Warm starts also cover a downward grid extension; preparation is reused.
  path <- sort(unique(c(.1, .01, .001, scales)), decreasing = TRUE) / sqrt(setup$n_ids)
  fit <- glmnet::glmnet(setup$design, setup$response, weights = setup$weights,
    family = setup$family, alpha = 1, standardize = FALSE, lambda = path,
    thresh = 1e-9, maxit = 1000000L)
  if (fit$jerr != 0 || length(fit$lambda) != length(path))
    stop("Saturated L1 penalty path failed to converge at every requested penalty.")
  values <- as.matrix(predict(fit, cbind(Matrix::Diagonal(length(setup$levels)), 0),
    type = "response"))
  selected <- .saturated_lambda_index(scales / sqrt(setup$n_ids), fit$lambda)
  if (anyNA(selected)) stop("Requested penalty was not fitted exactly.")
  values <- values[, selected, drop = FALSE]
  index <- match(validation_keys, setup$levels)
  answer <- matrix(setup$mean, length(index), length(scales))
  seen <- !is.na(index)
  answer[seen, ] <- values[index[seen], , drop = FALSE]
  answer
}

.saturated_penalty_selection <- function(Y, X, family, obsWeights, id,
    penalty_scales, penalty_folds, response_fit = NULL, max_extensions = 16L) {
  values <- as.character(id)
  codes <- suppressWarnings(abs(as.numeric(values)))
  nonnumeric <- which(!is.finite(codes))
  codes[nonnumeric] <- vapply(values[nonnumeric], function(value)
    strtoi(substr(digest::digest(value, algo = "xxhash32"), 1, 7), 16L), numeric(1))
  for (bit in 0:52) {
    labels <- as.integer(floor(codes / 2^bit) %% 2) + 1L
    if (length(unique(labels)) == 2L) break
  }
  if (is.function(penalty_folds)) labels <- penalty_folds(id)
  if (length(labels) != length(id) || anyNA(labels) || any(labels < 1 | labels != floor(labels)) ||
      length(unique(labels)) < 2L) stop("Penalty cross-validation needs valid person groups.")
  if (any(vapply(split(labels, id), function(group) length(unique(group)) != 1L, TRUE)))
    stop("Repeated person IDs must stay in the same penalty group.")
  times <- if ("time" %in% names(X)) X$time else rep(1L, nrow(X))
  plans <- lapply(sort(unique(labels)), function(label) {
    train <- labels != label; valid <- !train
    rebuilt <- if (is.null(response_fit)) list(train = Y[train], validation = Y[valid]) else
      response_fit(id[train], times[train], id[valid], times[valid])
    list(setup = .saturated_path_setup(rebuilt$train, X[train, , drop = FALSE],
      family, obsWeights[train], id[train]), keys = .saturated_keys(X[valid, , drop = FALSE]),
      response = rebuilt$validation, weights = obsWeights[valid])
  })
  evaluate <- function(scales) {
    losses <- numeric(length(scales))
    for (plan in plans) {
      pred <- .saturated_path_predict(plan$setup, plan$keys, scales)
      error <- if (family$family == "binomial") {
        pred[] <- pmin(1 - 1e-12, pmax(1e-12, pred))
        -(plan$response * log(pred) + (1 - plan$response) * log1p(-pred))
      } else (plan$response - pred)^2
      losses <- losses + colSums(plan$weights * error)
    }
    losses / sum(obsWeights)
  }
  cmbridge::select_penalty_grid(penalty_scales, evaluate, max_extensions = max_extensions)
}

#' Saturated L1 learner with training-only penalty selection
#' @rdname SL.saturated_l1
#' @param penalty_scales Positive initial scales; boundary minima extend the grid.
#' @param penalty_folds Optional function returning shared training-only fold labels from IDs.
#' @param penalty_max_extensions Maximum grid extensions before a boundary failure.
#' @export
SL.saturated_l1_cv <- function(Y, X, newX, family, obsWeights, id,
    penalty_scales = 10^seq(-5, 0, by = .25), penalty_folds = NULL,
    penalty_max_extensions = 16L, ...) {
  if (is.null(id)) id <- seq_along(Y)
  selected <- .saturated_penalty_selection(Y, X, family, obsWeights, id,
    penalty_scales, penalty_folds, max_extensions = penalty_max_extensions)
  result <- SL.saturated_l1(Y, X, newX, family, obsWeights, id,
                            penalty_scale = selected$scale)
  result$fit$penalty_cv <- selected$penalty_cv
  result$fit$penalty_boundary <- selected$boundary_history
  result$fit$penalty_status <- selected$status
  result
}

#' Constant regression within calendar time for misspecification comparisons
#' @param Y,X,newX,family,obsWeights,id Standard SuperLearner wrapper inputs.
#' @param ... Additional SuperLearner arguments.
#' @return Standard SuperLearner prediction and fit list.
#' @export
SL.time_mean <- function(Y, X, newX, family, obsWeights, id, ...) {
  key <- if ("time" %in% names(X)) as.character(X$time) else rep("all", length(Y))
  value <- vapply(split(seq_along(Y), key), function(ii) weighted.mean(Y[ii], obsWeights[ii]), 0)
  object <- list(value = value, mean = weighted.mean(Y, obsWeights), has_time = "time" %in% names(X))
  class(object) <- "SL.time_mean"
  list(pred = predict(object, newX), fit = object)
}

#' @export
predict.SL.time_mean <- function(object, newdata, ...) {
  key <- if (object$has_time) as.character(newdata$time) else rep("all", nrow(newdata))
  result <- unname(object$value[key])
  result[is.na(result)] <- object$mean
  result
}
