# Every cache key includes the complete training row set. Descendant responses
# are rebuilt inside each learner CV split, including penalty selection.
fh_engine <- function(z, root, root_labels, learners_y, learners_a, libraries,
                      bounds, control, estimator, nuisance_cache = NULL) {
  if (is.null(nuisance_cache)) nuisance_cache <- new.env(parent = emptyenv())
  signature <- digest::digest(list(data = z$data[root, , drop = FALSE], id = z$id[root],
    rows = root, labels = root_labels, rules = lapply(z$rules, function(rules) lapply(rules, fh_function_key)), nodes = z$epsilon,
    randomizer = fh_function_key(z$randomizer), ratio_draws = z$ratio_draws, mtp = z$mtp,
    trt = z$trt, health = z$health, baseline = z$baseline, k = z$k,
    libraries = libraries, learners_y = learners_y, learners_a = learners_a,
    bounds = bounds, control = control))
  if (!exists(signature, nuisance_cache, inherits = FALSE))
    assign(signature, new.env(parent = emptyenv()), nuisance_cache)
  cache <- get(signature, nuisance_cache, inherits = FALSE)
  split_cache <- new.env(parent = emptyenv())
  key <- function(rows, suffix) paste0(suffix, "_", digest::digest(sort(rows)))
  remember <- function(rows, suffix, fit) {
    # Depth-zero regressions and adjoints using only those regressions are
    # identical for SDR and TMLE. Deeper augmented responses are method-specific.
    if (grepl("^model[1-9]", suffix)) suffix <- paste0(suffix, ":", estimator)
    if (grepl("^adjoint", suffix)) {
      at <- as.integer(sub("adjoint([0-9]+)h.*", "\\1", suffix))
      horizon <- as.integer(sub(".*h([0-9]+)", "\\1", suffix))
      if (horizon - at > 1L) suffix <- paste0(suffix, ":", estimator)
    }
    name <- key(rows, suffix)
    if (!exists(name, cache, inherits = FALSE)) assign(name, fit(), cache)
    get(name, cache, inherits = FALSE)
  }
  split_rows <- function(rows) {
    name <- key(rows, "split")
    if (!exists(name, split_cache, inherits = FALSE)) {
      ids <- as.character(z$id[rows]); people <- sort(unique(ids))
      first <- match(people, ids)
      if (length(people) < 4L) stop("Nested fitting requires at least four training clusters.")
      if (!is.null(root_labels) && identical(rows, root)) assigned <- root_labels else {
        root_groups <- if (is.null(root_labels)) max(c(2L,
          control$.learners_outcome_folds, control$.learners_trt_folds)) else length(unique(root_labels))
        complete <- z$J[rows[first], , drop = FALSE] *
          as.matrix(z$data[rows[first], z$measurement[-1L], drop = FALSE])
        scarce <- min(colSums(complete)) <= 32L
        # Preserve measured training clusters in deep descendants, matching
        # the observed-history estimator's economical 2/3-group rule.
        V <- if (identical(rows, root)) root_groups else if (root_groups >= 3L && scarce) 3L else 2L
        if (V > length(people)) stop("Inner CV folds exceed the number of training clusters.")
        profile <- apply(cbind(z$J[rows[first], , drop = FALSE],
          z$data[rows[first], z$measurement, drop = FALSE]), 1, paste, collapse = "|")
        person_labels <- integer(length(people)); offset <- 0L
        for (stratum in sort(unique(profile))) {
          group <- which(profile == stratum)
          person_labels[group] <- (seq_along(group) - 1L + offset) %% V + 1L
          offset <- offset + length(group)
        }
        assigned <- person_labels[match(ids, people)]
      }
      if (length(unique(assigned)) < 2L) stop("At least two nonempty inner folds are required.")
      assign(name, assigned, split_cache)
    }
    get(name, split_cache, inherits = FALSE)
  }
  sl_splits <- function(ids) {
    rows <- root[match(unique(as.character(ids)), as.character(z$id[root]))]
    labels <- split_rows(rows)
    labels[match(as.character(ids), as.character(z$id[rows]))]
  }
  cm_library <- function(lib, rows, split = sl_splits) lapply(lib, function(candidate) {
    if (is.character(candidate)) candidate <- list(method = candidate)
    if (!is.null(candidate$control$penalty_scales)) {
      candidate$control$penalty_ids <- z$id[rows]
      candidate$control$penalty_folds <- split
    }
    candidate
  })
  beta_fit <- function(train, s) remember(train, paste0("beta", s), function() {
    eligible <- train[z$J[train, s]]
    M <- z$data[[z$measurement[s + 1L]]][eligible]
    if (!length(eligible)) stop("No complete training histories at time ", s, ".")
    if (all(M == 1)) return(list(known = TRUE))
    B <- as.matrix(fh_design(z, eligible, s))
    V <- as.matrix(fh_design(z, eligible, s, bridge = TRUE))
    cmbridge::fit_bridge_ensemble(B, V, M,
      library = cm_library(libraries$bridge, eligible),
      fold_id = split_rows(train)[match(eligible, train)], kernel = control$.cm_kernel,
      project_psd = control$.cm_project_psd)
  })
  beta_values <- function(train, rows, s) {
    out <- rep(1, length(rows)); fit <- beta_fit(train, s)
    observed <- which(z$J[rows, s] & z$data[[z$measurement[s + 1L]]][rows] == 1)
    if (!isTRUE(fit$known) && length(observed)) out[observed] <-
      predict(fit, as.matrix(fh_design(z, rows[observed], s, bridge = TRUE)))
    if (any(!is.finite(out)) || any(out < 1 - 1e-8))
      stop("Full-history bridges must be finite and >=1; use inverse_logit parameterizations.")
    out
  }
  ratio_fit <- function(train, s, branch) remember(train, paste0("ratio", s, branch), function() {
    eligible <- train[z$J[train, s]]
    natural <- fh_design(z, eligible, s)
    if (!control$.ratio_history) natural <- natural[, grepl("^trt[0-9]+$", names(natural)), drop = FALSE]
    natural$..fh_id <- z$id[eligible]; natural$..fh_y <- 0
    if (branch == "d") {
      shifted <- fh_design(z, eligible, s, fh_policy(z, eligible, s, "d"))
      multiplier <- 1L
    } else {
      multiplier <- z$ratio_draws
      shifted <- do.call(rbind, lapply(seq_len(multiplier), function(j) {
        draws <- z$randomizer(length(eligible), s)
        if (length(draws) != length(eligible) || any(!is.finite(draws)))
          stop("randomizer must return n finite draws.")
        fh_design(z, eligible, s, fh_policy(z, eligible, s, "t", draws))
      }))
    }
    if (!control$.ratio_history) shifted <- shifted[, grepl("^trt[0-9]+$", names(shifted)), drop = FALSE]
    shifted$..fh_id <- rep(z$id[eligible], multiplier); shifted$..fh_y <- 1
    stacked <- rbind(natural, shifted)
    labels <- split_rows(train)[match(stacked$..fh_id, z$id[train])]
    fit <- run_ensemble(stacked, "..fh_y", learners_a, "binomial", "..fh_id",
      list(fold_id = labels), control$.discrete, control$.info, penalty_folds = sl_splits)
    p <- predict(fit, natural)
    ratio <- p / (1 - p) / multiplier
    if (branch == "d" && !z$mtp) ratio[rowSums(as.matrix(z$data[eligible, z$trt[[s]], drop = FALSE]) ==
      fh_policy(z, eligible, s, "d")) != length(z$trt[[s]])] <- 0
    cap <- if (control$.trim < 1) as.numeric(quantile(ratio, control$.trim)) else Inf
    list(fit = fit, multiplier = multiplier, cap = cap)
  })
  nuisance <- function(train, rows) {
    beta <- rd <- rt <- matrix(1, length(rows), z$tau)
    fits <- vector("list", z$tau)
    for (s in seq_len(z$tau)) {
      beta[, s] <- beta_values(train, rows, s)
      fits[[s]] <- beta_fit(train, s)
      eligible <- which(z$J[rows, s])
      rd[, s] <- rt[, s] <- 0
      if (length(eligible)) {
        natural <- fh_design(z, rows[eligible], s)
        for (branch in c("d", "t")) {
          if (s == 1L && branch == "t") next
          fit <- ratio_fit(train, s, branch)
          p <- predict(fit$fit, natural)
          value <- pmin(p / (1 - p) / fit$multiplier, fit$cap)
          if (branch == "d" && !z$mtp) value[rowSums(as.matrix(
            z$data[rows[eligible], z$trt[[s]], drop = FALSE]) ==
            fh_policy(z, rows[eligible], s, "d")) != length(z$trt[[s]])] <- 0
          if (branch == "d") rd[eligible, s] <- value else rt[eligible, s] <- value
        }
      }
    }
    omega <- rd
    if (z$tau > 1L) for (s in 2:z$tau) omega[, s] <- omega[, s - 1L] *
      z$data[[z$measurement[s]]][rows] * (rd[, s] + (beta[, s - 1L] - 1) * rt[, s])
    if (any(!is.finite(omega))) stop("Nonfinite full-history propagation weights.")
    list(beta = beta, d = rd, t = rt, omega = omega, fits = fits)
  }
  terminal <- function(train, rows, s) {
    value <- numeric(length(rows))
    ii <- which(z$J[rows, s] & z$data[[z$measurement[s + 1L]]][rows] == 1)
    beta <- beta_values(train, rows, s)
    value[ii] <- beta[ii] * z$data[[z$outcome[s]]][rows[ii]]
    fh_check_bounds(value, bounds)
    value
  }
  predictions <- function(fit, rows, s) {
    out <- list(natural = numeric(length(rows)), d = numeric(length(rows)),
      t = numeric(length(rows)), t_nodes = matrix(0, length(rows), length(z$epsilon[[s]])))
    ii <- which(z$J[rows, s]); active <- rows[ii]
    if (!length(ii)) return(out)
    out$natural[ii] <- predict(fit, fh_design(z, active, s), tol = NULL)
    out$d[ii] <- predict(fit, fh_design(z, active, s, fh_policy(z, active, s, "d")), tol = NULL)
    if (s == 1L) return(out) # R1 is known to be one; no baseline continuation branch.
    for (j in seq_along(z$epsilon[[s]])) out$t_nodes[ii, j] <- predict(fit,
      fh_design(z, active, s, fh_policy(z, active, s, "t", rep(z$epsilon[[s]][j], length(ii)))), tol = NULL)
    out$t <- rowMeans(out$t_nodes)
    out
  }
  response <- function(train, rows, times, depth) {
    value <- numeric(length(rows))
    for (s in unique(times)) {
      ii <- which(times == s); rr <- rows[ii]
      if (depth == 0L) value[ii] <- terminal(train, rr, s) else {
        fit <- model(train, depth - 1L)
        p <- predictions(fit, rr, s + 1L)
        beta <- beta_values(train, rr, s)
        value[ii] <- fh_q(z, rr, s + 1L, p, beta)
        if (estimator == "sdr") {
          future <- response(train, rr, rep(s + 1L, length(rr)), depth - 1L)
          nu <- nuisance(train, rr)
          local_weight <- z$data[[z$measurement[s + 1L]]][rr] *
            (nu$d[, s + 1L] + (beta - 1) * nu$t[, s + 1L])
          value[ii] <- value[ii] + local_weight * (future - p$natural)
        }
      }
    }
    fh_check_bounds(value, bounds)
    value
  }
  model <- function(train, depth) remember(train, paste0("model", depth), function() {
    times <- seq_len(z$tau - depth)
    designs <- lapply(times, function(s) {
      rows <- train[z$J[train, s]]
      x <- fh_design(z, rows, s); x$..fh_y <- 0
      x
    })
    data <- do.call(rbind, designs)
    # Internal row keys, rather than cluster IDs, locate responses for repeated
    # observations. The learner ID remains the user's cluster identifier.
    data$..fh_row <- unlist(lapply(times, function(s) train[z$J[train, s]]))
    # run_ensemble forwards IDs and times to response_fit. A unique row key is
    # safe only when every ID has one row; for clusters use a compound key and
    # preserve their shared folds via the explicit labels below.
    data$..fh_id <- as.character(data$..fh_row)
    labels <- split_rows(train)[match(data$..fh_row, train)]
    data$..fh_row <- NULL
    rebuild <- function(ids, times, valid_ids, valid_times) {
      train_rows <- sort(unique(as.integer(ids)))
      # All rows in a cluster enter or leave together, including rows that do
      # not contribute to this diagonal's complete-state regression.
      scope <- train[z$id[train] %in% z$id[train_rows]]
      list(train = response(scope, as.integer(ids), times, depth),
        validation = response(scope, as.integer(valid_ids), valid_times, depth))
    }
    run_ensemble(data, "..fh_y", learners_y, "continuous", "..fh_id",
      list(fold_id = labels), control$.discrete, control$.info,
      response_fit = rebuild, penalty_folds = function(ids) {
        rr <- unique(as.integer(ids)); ll <- split_rows(rr)
        ll[match(as.integer(ids), rr)]
      })
  })
  loading <- function(train, rows, s, h) {
    nu <- nuisance(train, rows)
    if (s == h) {
      value <- z$data[[z$outcome[h]]][rows]
      value[!is.finite(value)] <- 0
    } else value <- predictions(model(train, h - s - 1L), rows, s + 1L)$t
    nu$omega[, s] * value
  }
  adjoint <- function(train, rows, s, h) {
    fit <- remember(train, paste0("adjoint", s, "h", h), function() {
      eligible <- train[z$J[train, s]]
      M <- z$data[[z$measurement[s + 1L]]][eligible]
      if (all(M == 1)) return(list(known = TRUE))
      # A conditional-moment subset is smaller than the underlying nuisance
      # training sample. Keep a map from its splits to the full person splits;
      # otherwise conditioning on future complete histories would bias earlier
      # bridges and regressions when rebuilding the loading.
      scopes <- new.env(parent = emptyenv())
      scope_key <- function(e) digest::digest(sort(unique(e)))
      assign(scope_key(eligible), train, scopes)
      register <- function(e, scope) {
        labels <- split_rows(scope)
        ee <- labels[match(e, scope)]
        for (label in unique(ee)) assign(scope_key(e[ee != label]),
          scope[labels != label], scopes)
        ee
      }
      labels <- register(eligible, train)
      penalty_split <- function(ids) {
        e <- eligible[as.character(z$id[eligible]) %in% as.character(ids)]
        name <- scope_key(e)
        if (!exists(name, scopes, inherits = FALSE)) stop("Missing full-history loading split.")
        labels <- register(e, get(name, scopes, inherits = FALSE))
        labels[match(as.character(ids), as.character(z$id[e]))]
      }
      phi <- function(ti, vi) {
        e <- eligible[ti]
        name <- scope_key(e)
        if (!exists(name, scopes, inherits = FALSE)) stop("Missing adjoint training scope.")
        scope <- get(name, scopes, inherits = FALSE)
        ll <- split_rows(scope)
        training <- numeric(length(ti))
        for (label in unique(ll)) {
          at <- which(z$id[e] %in% z$id[scope[ll == label]])
          if (length(at)) training[at] <- loading(scope[ll != label], e[at], s, h)
        }
        list(train = training, validation = loading(scope, eligible[vi], s, h))
      }
      cmbridge::fit_adjoint_ensemble(as.matrix(fh_design(z, eligible, s)),
        as.matrix(fh_design(z, eligible, s, bridge = TRUE)), M, phi,
        library = cm_library(libraries$adjoint, eligible, penalty_split),
        fold_id = labels, kernel = control$.cm_kernel, project_psd = control$.cm_project_psd)
    })
    out <- numeric(length(rows)); ii <- which(z$J[rows, s])
    if (!isTRUE(fit$known) && length(ii)) out[ii] <- predict(fit, as.matrix(fh_design(z, rows[ii], s)))
    list(model = fit, value = out)
  }
  list(model = model, nuisance = nuisance, predictions = predictions,
       terminal = terminal, response = response, adjoint = adjoint,
       splits = split_rows, counts = function() table(sub("_.*", "", ls(cache))))
}

# Hash the bindings actually referenced by a callback. Serializing its entire
# enclosing call frame can capture the mutable cache itself, defeating reuse
# and making signatures depend on unrelated validation observations.
fh_function_key <- function(fun, seen = character()) {
  result <- list(formals = formals(fun), body = body(fun))
  globals <- codetools::findGlobals(fun, merge = FALSE)
  names <- unique(c(globals$variables, globals$functions))
  env <- environment(fun)
  values <- list()
  for (name in setdiff(names, seen)) {
    if (!exists(name, env, inherits = TRUE)) next
    value <- get(name, env, inherits = TRUE)
    if (is.function(value)) {
      parent <- environment(value)
      if (is.primitive(value) || identical(parent, baseenv()) || isNamespace(parent)) next
      values[[name]] <- fh_function_key(value, c(seen, name))
    } else values[[name]] <- value
  }
  result$bindings <- values
  result
}
