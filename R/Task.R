#' @importFrom R6 R6Class
LmtpTask <- R6::R6Class(
  "LmtpTask",
  public = list(
    natural = NULL,
    shifted = NULL,
    vars = NULL,
    id = NULL,
    n = NULL,
    tau = NULL,
    outcome_type = NULL,
    survival = NULL,
    folds = NULL,
    learner_folds = NULL,
    weights = NULL,
    initialize = function(data, shifted, A, Y, L, W, C, D, k,
                          id, outcome_type, folds, weights, bounds = NULL,
                          learner_folds = NULL) {
      # Identify tau
      self$tau <- private$tau_is(Y, A, C)
      self$n <- nrow(data)

      # Create Vars object
      self$vars <- LmtpVars$new(W, L, A, C, D, Y, outcome_type, self$tau, k)

      # Additional checks
      assert_numeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
      if (!is.list(folds)) checkmate::assert_count(folds, positive = TRUE)

      # Set outcome types
      self$outcome_type <- ifelse(outcome_type %in% c("binomial", "survival"), "binomial", "continuous")
      self$survival <- outcome_type == "survival"

      # Set outcome bounds
      assert_numeric(bounds, len = 2, unique = TRUE, sorted = TRUE, finite = TRUE, null.ok = TRUE)
      private$bounds <- if (is.null(bounds)) private$y_bounds(data[, self$vars$Y]) else bounds
      observed_y <- unlist(data[, self$vars$Y, drop = FALSE])
      observed_y <- observed_y[!is.na(observed_y)]
      if (any(observed_y < private$bounds[1] | observed_y > private$bounds[2])) {
        stop("Observed outcomes must lie inside `bounds`.")
      }

      # Add cluster IDs
      self$id <- private$add_ids(data, id)

      # Modify data.frames to work with estimators
      self$natural <- private$as_lmtp_data(data)
      self$shifted <- private$as_lmtp_data(shifted)

      assert_lmtp_data(self)

      # Make folds for cross-fitting
      self$folds <- private$make_folds(folds)
      self$learner_folds <- validate_learner_folds(learner_folds, self$folds, self$id)

      # Add or normalize weights
      self$weights <- private$make_weights(weights)
    },

    scale = function(x) {
      (x - private$bounds[1]) / (private$bounds[2] - private$bounds[1])
    },

    rescale = function(x) {
      (x*(private$bounds[2] - private$bounds[1])) + private$bounds[1]
    },

    observed = function(data, t) {
      if (is.null(self$vars$C)) {
        return(rep(TRUE, nrow(data)))
      }

      if (t == 0) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$C[t]]] == 1
    },

    at_risk_D = function(data, t) {
      if (t > self$tau) {
        return(rep(TRUE, nrow(data)))
      }

      # always at risk at first time point
      if (t == 0) {
        return(rep(TRUE, nrow(data)))
      }

      if (is.null(self$vars$D)) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$D[t]]] == 0 & !is.na(data[[self$vars$D[t]]])
    },

    at_risk_N = function(data, t) {
      if (t > self$tau) {
        return(rep(TRUE, nrow(data)))
      }

      if (is.null(self$vars$N) || isFALSE(self$survival)) {
        return(rep(TRUE, nrow(data)))
      }

      # always at risk at first time point
      if (t == 0) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$N[t]]] == 1 & !is.na(data[[self$vars$N[t]]])
    },

    R = function(data, t) {
      if (t > self$tau) {
        return(rep(TRUE, nrow(data)))
      }

      # If not a survival outcome, always at risk
      if (isFALSE(self$survival)) {
        return(rep(TRUE, nrow(data)))
      }

      self$at_risk_D(data, t - 1) & self$at_risk_N(data, t - 1)
    }

  ),
  private = list(
    bounds = NULL,
    tau_is = function(Y, A, C) {
      if (!(length(Y) > 1)) {
        return(max(length(A), length(C)))
      }
      length(Y)
    },

    add_ids = function(data, id) {
      if (is.null(id)) {
        return(1:nrow(data))
      }
      data[[id]]
    },

    make_folds = function(V) {
      id <- self$natural$..i..lmtp_id
      if (is.list(V)) return(validate_folds(V, self$n, id))

      if (length(unique(id)) == self$n & self$outcome_type == "binomial") {
        strata <- self$natural[[last(self$vars$Y)]]
        strata[is.na(strata)] <- 2
        folds <- make_folds(self$natural, V = V, strata_ids = strata)
      } else {
        folds <- make_folds(self$natural, cluster_ids = id, V = V)
      }

      validate_folds(folds, self$n, id)
    },

    y_bounds = function(x) {
      if (self$outcome_type == "binomial") {
        return(c(0, 1))
      }
      bounds <- range(unlist(x), na.rm = TRUE)
      if (any(!is.finite(bounds))) stop("At least one observed outcome is required.")
      if (bounds[1] == bounds[2]) bounds <- bounds + c(-0.5, 0.5)
      bounds
    },

    as_lmtp_data = function(x) {
      data <- data.table::copy(as.data.frame(x))
      data$..i..lmtp_id <- self$id
      data <- fix_censoring_ind(data, self$vars$C)

      if (self$survival) {
        for (y in c(self$vars$N, self$vars$Y)) {
          data.table::set(data, j = y, value = convert_to_surv(data[[y]]))
        }
      }

      data[, self$vars$Y] <- self$scale(data[, self$vars$Y])
      data
    },

    make_weights = function(weights) {
      if (!is.null(weights)) {
        if (is_normalized(weights)) {
          return(weights)
        } else {
          # Normalize weights
          return(weights / mean(weights))
        }
      }

      rep(1, self$n)
    }
  )
)
