# Raw conditional-moment fits are never clipped. The longitudinal estimator
# can request a wider finite interval for its applied bridge predictions.
bridge_predict <- function(model, newdata, control) {
  value <- predict(model, newdata)
  limits <- control$.bridge_prediction_bounds
  if (is.null(limits)) return(value)
  pmax(limits[1L], pmin(limits[2L], value))
}

# Fold-specific construction of sequential-regression responses. All ancestors
# of a learner's training response are fitted within that learner's training IDs.
make_bridge_nested_folds <- function(ids, measurement, labels) {
  ids <- as.character(ids); measurement <- as.matrix(measurement)
  cache <- new.env(parent = emptyenv())
  groups <- length(unique(labels))
  function(current_ids) {
    current <- as.character(current_ids)
    kept <- which(ids %in% current)
    if (length(unique(ids[kept])) < 2L)
      stop("Nested cross-validation requires at least two training clusters.")
    key <- digest::digest(kept)
    if (!exists(key, cache, inherits = FALSE)) {
      if (identical(kept, seq_along(ids))) value <- labels else {
        people <- unique(ids[kept])
        observed <- measurement[kept, , drop = FALSE]
        if (anyDuplicated(ids[kept])) observed <-
          rowsum(observed, ids[kept], reorder = FALSE) > 0
        counts <- colSums(observed)
        incomplete <- counts[counts < length(people)]
        # Larger additional fits keep the economical two-way split. When
        # measured training clusters are scarce, the configured extra group
        # keeps more measured observations in every descendant training fit.
        # This uses only measurement indicators inside the current training set.
        additional_groups <- if (groups >= 3L && length(incomplete) && min(incomplete) <= 16L)
          groups else 2L
        profile <- do.call(paste, c(as.data.frame(observed), sep = "|"))
        assigned <- integer(length(people)); offset <- 0L
        for (stratum in sort(unique(profile))) {
          group <- which(profile == stratum)
          group <- group[order(people[group])]
          assigned[group] <- (seq_along(group) - 1L + offset) %% additional_groups + 1L
          offset <- offset + length(group)
        }
        value <- assigned[match(ids[kept], people)]
      }
      assign(key, value, cache)
    }
    value <- get(key, cache, inherits = FALSE)
    value[match(current, ids[kept])]
  }
}

bridge_response_engine <- function(original, encoded, task, fold, measurement, health,
                                    bridge_library, learners, learners_trt, mtp, control,
                                    estimator, split_function, shared_nuisance = NULL) {
  rows <- task$folds[[fold]]$training_set
  control$.nested_split_function <- split_function
  people <- as.character(task$id)
  beta_cache <- ratio_cache <- prefix_cache <- response_cache <- NULL
  common <- if (is.environment(shared_nuisance$response_nuisance_cache))
    shared_nuisance$response_nuisance_cache else new.env(parent = emptyenv())
  clean_library <- lapply(bridge_library, function(candidate) {
    if (is.list(candidate)) for (name in c("penalty_ids", "penalty_folds")) candidate$control[[name]] <- NULL
    candidate
  })
  signature <- digest::digest(list(task$id[rows], encoded[rows, , drop = FALSE],
    clean_library, learners_trt, control$.cm_kernel, control$.cm_project_psd,
    control$.discrete, control$.ratio_history))
  if (!identical(common$signature, signature) || !identical(common$splits, split_function)) {
    common$beta <- common$ratio <- NULL
    common$signature <- signature; common$splits <- split_function
  }
  if (is.null(common$beta)) common$beta <- new.env(parent = emptyenv())
  if (is.null(common$ratio)) common$ratio <- new.env(parent = emptyenv())
  beta_cache <- common$beta
  ratio_cache <- common$ratio
  prefix_cache <- new.env(parent = emptyenv())
  response_cache <- new.env(parent = emptyenv())
  select_rows <- function(ids) rows[people[rows] %in% as.character(ids)]
  key_for <- function(train, extra = NULL) digest::digest(list(train, extra))
  map_rows <- function(ids, times) {
    out <- integer(length(ids))
    for (time in unique(times)) {
      ii <- which(times == time)
      selected <- select_rows(ids[ii])
      if (length(selected) != length(ii) ||
          !identical(people[selected], as.character(ids[ii])))
        stop("Response row order does not match the supplied training clusters.")
      out[ii] <- selected
    }
    out
  }
  beta_design <- lapply(seq_len(task$tau), function(s) {
    H <- task$vars$history("A", s); A <- task$vars$A[[s]]
    list(B = as.matrix(encoded[, c(H, A), drop = FALSE]),
      V = as.matrix(cbind(encoded[, H, drop = FALSE], original[, health[[s]], drop = FALSE])),
      M = original[[measurement[s]]])
  })
  fit_beta <- function(train, s) {
    key <- key_for(train, s)
    if (!exists(key, beta_cache, inherits = FALSE)) {
      z <- beta_design[[s]]
      if (all(z$M == 1)) model <- NULL else {
        lib <- bridge_library
        for (a in seq_along(lib)) if (!is.null(lib[[a]]$control$penalty_scales)) {
          lib[[a]]$control$penalty_ids <- task$id[train]
          lib[[a]]$control$penalty_folds <- split_function
        }
        # A single candidate has weight one. Its full-data model is unchanged
        # by omitting the otherwise unused ensemble-weight validation fits.
        if (length(lib) == 1L && lib[[1]]$method %in% c("saturated_l1", "sieve_md") &&
            !is.null(lib[[1]]$control$penalty_scales) && identical(control$.cm_kernel, "cell")) {
          candidate <- lib[[1]]; scales <- candidate$control$penalty_scales
          at_penalty <- function(scale, ids) {
            ctrl <- candidate$control
            for (name in c("penalty_scales", "penalty_ids", "penalty_folds",
              "penalty_max_extensions", "penalty_extension_factor", "penalty_extension_points"))
              ctrl[[name]] <- NULL
            if (candidate$method == "sieve_md") ctrl$lambda <- scale else
              ctrl$penalty <- scale / sqrt(length(unique(ids)))
            ctrl
          }
          labels <- split_function(task$id[train])
          paired_count <- sum(vapply(unique(labels), function(label) {
            vi <- train[labels == label]
            scoring <- cmbridge::cell_moment_gram(rep(0, length(vi)), z$B[vi, , drop = FALSE])
            as.integer(attr(scoring, "cell_scoring")$paired_rows)
          }, 0L))
          evaluate <- function(requested) vapply(requested, function(scale) {
            loss <- 0
            for (label in unique(labels)) {
              ti <- train[labels != label]; vi <- train[labels == label]
              if (!any(z$M[ti] == 1)) stop("No measured training observations in a bridge penalty fit.")
              ctrl <- at_penalty(scale, task$id[ti])
              f <- cmbridge::fit_bridge(z$B[ti, , drop = FALSE], z$V[ti, , drop = FALSE], z$M[ti],
                candidate$method, ctrl)
              residual <- rep(-1, length(vi)); observed <- which(z$M[vi] == 1)
              residual[observed] <- predict(f, z$V[vi[observed], , drop = FALSE]) - 1
              # Retain current A as well as every observed-history variable.
              loss <- loss + length(vi) * cmbridge::cell_moment_gram(residual,
                z$B[vi, , drop = FALSE])[1, 1]
            }
            loss / length(train)
          }, 0)
          selected <- cmbridge::select_penalty_grid(scales, evaluate,
            max_extensions = if (is.null(candidate$control$penalty_max_extensions)) 16L else candidate$control$penalty_max_extensions,
            extension_factor = if (is.null(candidate$control$penalty_extension_factor)) 10 else candidate$control$penalty_extension_factor,
            extension_points = if (is.null(candidate$control$penalty_extension_points)) 4L else candidate$control$penalty_extension_points)
          model <- cmbridge::fit_bridge(z$B[train, , drop = FALSE], z$V[train, , drop = FALSE], z$M[train],
            candidate$method, at_penalty(selected$scale, task$id[train]))
          model$penalty_cv <- selected$penalty_cv
          model$penalty_cv$scored_rows <- length(train)
          model$penalty_cv$paired_rows <- paired_count
          model$penalty_cv$informative <- paired_count > 0L
          model$penalty_boundary <- selected$boundary_history
          model$penalty_status <- selected$status
          model$penalty_parameter <- if (candidate$method == "sieve_md") "lambda" else "penalty"
          model$penalty_ids <- task$id[train]
          model$penalty_fold_id <- labels
        } else model <- cmbridge::fit_bridge_ensemble(z$B[train, , drop = FALSE], z$V[train, , drop = FALSE],
          z$M[train], lib, fold_id = split_function(task$id[train]), kernel = control$.cm_kernel,
          project_psd = !identical(control$.cm_project_psd, FALSE))
      }
      assign(key, list(model = model, horizon = s, training_ids = task$id[train]), beta_cache)
    }
    get(key, beta_cache, inherits = FALSE)$model
  }
  terminal <- function(train, valid) {
    out <- matrix(0, length(valid), task$tau)
    for (s in seq_len(task$tau)) {
      observed <- which(beta_design[[s]]$M[valid] == 1)
      if (!length(observed)) next
      model <- fit_beta(train, s)
      beta <- if (is.null(model)) rep(1, length(observed)) else
        bridge_predict(model, beta_design[[s]]$V[valid[observed], , drop = FALSE], control)
      out[observed, s] <- original[[task$vars$Y[s]]][valid[observed]] * beta
    }
    task$scale(out)
  }
  slim_sl <- function(fit) {
    if (is.null(fit)) return(NULL)
    fit[c("Y", "Z", "SL.predict", "library.predict", "validRows", "call")] <- NULL
    fit
  }
  # Only fits on exactly the same training IDs can be reused. Sequential
  # regressions/responses remain estimator-specific, never in this common cache.
  if (!is.null(shared_nuisance$beta_fits)) for (s in seq_len(task$tau))
    if (!is.null(shared_nuisance$beta_fits[[s]])) assign(key_for(rows, s),
      list(model = shared_nuisance$beta_fits[[s]], horizon = s, training_ids = task$id[rows]), beta_cache)
  if (!is.null(shared_nuisance$ratio_fits) && length(shared_nuisance$ratio_fits) == task$tau &&
      all(vapply(shared_nuisance$ratio_fits, function(fit) inherits(fit, "SuperLearner"), logical(1)))) assign(key_for(rows),
    lapply(shared_nuisance$ratio_fits, slim_sl), ratio_cache)
  fit_ratios <- function(train) {
    key <- key_for(train)
    if (!exists(key, ratio_cache, inherits = FALSE)) {
      sub <- task$clone(deep = TRUE)
      sub$folds <- list(list(training_set = train, validation_set = rows))
      sub$learner_folds <- list(split_function(task$id[train]))
      ctrl <- control; ctrl$.return_full_fits <- TRUE
      fits <- estimate_r(sub, 1, learners_trt, mtp, ctrl, function(...) NULL)$fits
      if (any(vapply(fits, function(f) any(f$errorsInCVLibrary) || any(f$errorsInLibrary), FALSE)))
        stop("A density-ratio learner failed inside response construction.")
      assign(key, lapply(fits, slim_sl), ratio_cache)
    }
    get(key, ratio_cache, inherits = FALSE)
  }
  predict_ratios <- function(train, valid) {
    fits <- fit_ratios(train); ans <- matrix(0, length(valid), task$tau)
    for (s in seq_len(task$tau)) {
      A <- if (length(task$vars$A) > 1L) task$vars$A[[s]] else task$vars$A[[1]]
      value <- predict(fits[[s]], task$natural[valid, , drop = FALSE])
      followed <- followed_rule(task$natural[valid, , drop = FALSE], task$shifted[valid, , drop = FALSE], A, mtp)
      if (!mtp) value[followed] <- pmax(value[followed], .5)
      ans[, s] <- value * followed / (1 - pmin(value, .999))
    }
    ans
  }
  response <- prefix <- NULL
  prefix <- function(train, depth) {
    key <- key_for(train, depth)
    if (!exists(key, prefix_cache, inherits = FALSE)) {
      natural <- pivot(task$natural[train, , drop = FALSE], task$vars)
      fits <- vector("list", depth)
      for (d in seq_len(depth)) {
        keep <- natural$time <= task$tau - d + 1L
        vars <- setdiff(names(natural), c("..i..C_1", "..i..C_1_lag", "..i..wide_id", "..i..N", "..i..D_1", "..i..R"))
        if (d == task$tau) vars <- setdiff(vars, "time")
        for (name in setdiff(vars, c("..i..lmtp_id", "..i..Y_1", grep("^..i..A", vars, value = TRUE))))
          if (all(duplicated(natural[keep, name])[-1L])) vars <- setdiff(vars, name)
        rebuild <- local({current_depth <- d
          function(ids, times, valid_ids, valid_times) response(ids, times, valid_ids, valid_times, current_depth)})
        fitted <- run_ensemble(natural[keep, vars, drop = FALSE], "..i..Y_1", learners,
          "continuous", "..i..lmtp_id", list(fold_id = split_function(natural$..i..lmtp_id[keep])),
          control$.discrete, control$.info, rebuild, split_function)
        if (any(fitted$errorsInCVLibrary) || any(fitted$errorsInLibrary))
          stop("A sequential learner failed inside response construction.")
        fits[[d]] <- slim_sl(fitted)
        assign(key_for(train, d), fits[seq_len(d)], prefix_cache)
      }
      assign(key, fits, prefix_cache)
    }
    get(key, prefix_cache, inherits = FALSE)
  }
  evaluate <- function(train, ids, times, depth) {
    if (!length(ids)) return(numeric())
    valid <- map_rows(ids, times)
    if (depth == 1L) {
      eta <- terminal(train, valid)
      return(eta[cbind(seq_along(valid), times)])
    }
    unique_rows <- select_rows(ids)
    index <- match(valid, unique_rows)
    natural <- pivot(task$natural[unique_rows, , drop = FALSE], task$vars)
    shifted <- natural; policy <- pivot(task$shifted[unique_rows, , drop = FALSE], task$vars)
    A <- grep("^..i..A", names(shifted), value = TRUE)
    A <- grep("_lag", A, value = TRUE, invert = TRUE)
    shifted[, A] <- policy[, A]
    qn <- lapply(seq_len(task$tau), function(s) matrix(0, length(unique_rows), s))
    qd <- lapply(seq_len(task$tau), function(s) matrix(0, length(unique_rows), s + 1L))
    eta <- terminal(train, unique_rows)
    for (s in seq_len(task$tau)) qd[[s]][, s + 1L] <- eta[, s]
    fits <- prefix(train, depth - 1L)
    for (d in seq_along(fits)) {
      tol <- if (estimator == "tmle") control$.bound else NULL
      qn <- update_m(qn, fits[[d]], natural, d, task$tau, FALSE, tol)
      qd <- update_m(qd, fits[[d]], shifted, d, task$tau, FALSE, tol)
    }
    ratios <- if (estimator == "sdr") predict_ratios(train, unique_rows) else NULL
    out <- numeric(length(ids))
    for (time in unique(times)) {
      s <- time + depth - 1L; ii <- which(times == time)
      value <- if (estimator == "sdr") eif(ratios, matrix(1, length(unique_rows), task$tau),
        qd[[s]], qn[[s]], t = time + 1L, tau = s) else qd[[s]][, time + 1L]
      out[ii] <- value[index[ii]]
    }
    out
  }
  response <- function(ids, times, valid_ids, valid_times, depth) {
    train <- select_rows(ids)
    key <- digest::digest(list(train, as.character(ids), times, as.character(valid_ids), valid_times, depth))
    if (!exists(key, response_cache, inherits = FALSE)) {
      labels <- split_function(ids); value <- numeric(length(ids))
      for (label in unique(labels)) {
        vi <- which(labels == label); ti <- select_rows(ids[labels != label])
        value[vi] <- evaluate(ti, ids[vi], times[vi], depth)
      }
      ans <- list(train = value, validation = evaluate(train, valid_ids, valid_times, depth))
      if (any(!is.finite(c(ans$train, ans$validation)))) stop("Nonfinite rebuilt regression responses.")
      assign(key, ans, response_cache)
    }
    get(key, response_cache, inherits = FALSE)
  }
  outer_ratios <- function(valid, oof = FALSE) {
    if (!oof) return(predict_ratios(rows, valid))
    labels <- split_function(task$id[rows]); out <- matrix(NA_real_, length(rows), task$tau)
    for (label in unique(labels)) {
      vi <- which(labels == label)
      out[vi, ] <- predict_ratios(rows[labels != label], rows[vi])
    }
    out
  }
  list(response = response, splits = split_function, outer_ratios = outer_ratios)
}
