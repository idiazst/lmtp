# Joint categorical cells supply all interactions, including time. The support
# is learned from training predictors only; an unseen cell uses the training
# mean and is recorded in the fitted object.
#' Saturated L1 SuperLearner wrapper
#' @param Y,X,newX,family,obsWeights,id Standard SuperLearner wrapper inputs.
#' @param penalty_scale Penalty multiplied by the inverse square root of IDs.
#' @param ... Additional SuperLearner arguments.
#' @return Standard SuperLearner prediction and fit list.
#' @export
SL.saturated_l1 <- function(Y, X, newX, family, obsWeights, id,
                            penalty_scale = 0.01, ...) {
  if (is.null(id)) id <- seq_along(Y)
  key <- function(x) if (!ncol(x)) rep("constant", nrow(x)) else
    do.call(paste, c(as.data.frame(x), sep = "|"))
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
    alpha = 1, standardize = FALSE, lambda = penalty_scale / sqrt(length(unique(id))),
    thresh = 1e-9, maxit = 100000L)
  object <- list(object = fit, levels = levels, mean = weighted.mean(Y, obsWeights),
                 penalty = penalty_scale / sqrt(length(unique(id))),
                 family = family$family)
  class(object) <- "SL.saturated_l1"
  list(pred = predict(object, newX), fit = object)
}

#' @export
predict.SL.saturated_l1 <- function(object, newdata, ...) {
  key <- do.call(paste, c(as.data.frame(newdata), sep = "|"))
  j <- match(key, object$levels)
  out <- rep(object$mean, nrow(newdata))
  keep <- which(!is.na(j))
  if (!is.null(object$object) && length(keep)) {
    design <- Matrix::sparseMatrix(i = seq_along(keep), j = j[keep], x = 1,
                                  dims = c(length(keep), length(object$levels)))
    out[keep] <- as.numeric(predict(object$object, cbind(design, 0),
                                   s = object$penalty, type = "response"))
  }
  out
}

#' @rdname SL.saturated_l1
#' @export
SL.saturated_l1_small <- function(...) SL.saturated_l1(..., penalty_scale = 0.001)
#' @rdname SL.saturated_l1
#' @export
SL.saturated_l1_large <- function(...) SL.saturated_l1(..., penalty_scale = 0.1)
