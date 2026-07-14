theta_lmtp <- function(task, estimates, density_ratios, shift, is_sdr, strata = NULL) {
  if (is_sdr) {
    theta <- fmean(estimates$uncentered_eif, w = task$weights)
  } else {
    theta <- fmean(estimates$predictions[, 1], w = task$weights)
  }

  influence_function <- task$rescale(estimates$uncentered_eif)
  theta <- task$rescale(theta)
  stratified <- stratified_lmtp_estimates(task, estimates, strata, is_sdr)

  out <- list(
    estimator = ifelse(is_sdr, "SDR", "TMLE"),
    estimate = ife(theta, influence_function, task$weights, as.character(task$id)),
    stratified_estimates = stratified$estimates,
    strata_table = stratified$table,
    n = length(task$weights),
    shift = shift,
    outcome_reg = task$rescale(estimates$predictions),
    density_ratios = density_ratios$density_ratios,
    fits_outcome = estimates$learner_outcome_summary,
    fits_treatment = density_ratios$learner_treatment_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}

stratified_lmtp_estimates <- function(task, estimates, strata, is_sdr) {
  if (is.null(strata)) {
    return(list(estimates = NULL, table = NULL))
  }

  strata <- as.data.frame(strata)

  if (nrow(strata) != length(task$weights)) {
    stop("`strata` must have one row per observation.", call. = FALSE)
  }

  if (any(vapply(strata, is.list, logical(1)))) {
    stop("Stratification variables must be atomic vectors.", call. = FALSE)
  }

  point_contribution <- if (is_sdr) {
    estimates$uncentered_eif
  } else {
    estimates$predictions[, 1]
  }

  point_contribution <- task$rescale(point_contribution)
  uncentered_eif <- task$rescale(estimates$uncentered_eif)

  marginal_estimates <- list()
  marginal_rows <- list()
  output_index <- 0L

  for (variable in names(strata)) {
    x <- strata[[variable]]
    missing <- is.na(x)

    # Retain observed levels in their order of first appearance and include
    # missing values as one explicit level.
    first_indices <- which(!missing & !duplicated(x))
    if (any(missing)) {
      first_indices <- c(first_indices, which(missing)[1L])
    }
    first_indices <- sort(first_indices)

    for (first_index in first_indices) {
      value <- x[first_index]

      in_level <- if (is.na(value)) {
        is.na(x)
      } else {
        !is.na(x) & x == value
      }

      probability <- fmean(as.numeric(in_level), w = task$weights)
      if (!is.finite(probability) || probability <= 0) {
        stop(
          "Every observed level of each stratification variable must have positive total weight.",
          call. = FALSE
        )
      }

      theta <- fmean(
        point_contribution[in_level],
        w = task$weights[in_level]
      )

      # This is the EIF for the intervention mean conditional on V_j = v.
      # Other requested stratification variables are marginalized over.
      marginal_eif <- as.numeric(in_level) / probability *
        (uncentered_eif - theta)

      level_label <- if (is.na(value)) "<NA>" else as.character(value)
      label <- paste0(variable, "=", level_label)
      output_index <- output_index + 1L

      marginal_estimates[[output_index]] <- ife(
        theta,
        marginal_eif,
        task$weights,
        as.character(task$id)
      )

      marginal_rows[[output_index]] <- data.frame(
        variable = variable,
        level = level_label,
        probability = probability,
        n = sum(in_level),
        label = label,
        stringsAsFactors = FALSE
      )
    }
  }

  strata_table <- do.call(rbind, marginal_rows)
  strata_table$label <- make.unique(strata_table$label, sep = "__")
  names(marginal_estimates) <- strata_table$label
  rownames(strata_table) <- strata_table$label

  list(estimates = marginal_estimates, table = strata_table)
}

theta_ltmle <- function(task, estimates, propensity_scores, levels, trt_balance, cens_balance) {
  theta <- sapply(estimates, function(x) fmean(x$predictions[, 1], w = task$weights))

  # Rescale estimates
  theta <- sapply(theta, function(x) task$rescale(x))
  influence_functions <- lapply(estimates, function(x) task$rescale(x$uncentered_eif))

  # Create 'ife' objects
  ifes <- lapply(seq_along(theta), function(i) {
    ife::ife(theta[i], influence_functions[[i]], task$weights, as.character(task$id))
  })
  names(ifes) <- levels

  out <- list(
    estimator = "TMLE",
    estimates = ifes,
    outcome_reg = lapply(estimates, function(x) task$rescale(x$predictions)),
    propensity_scores = propensity_scores$propensity_score,
    prob_observed = propensity_scores$prob_observed,
    balance = list(
      treatment = setNames(trt_balance, paste0("time ", seq_len(length(task$vars$A)))),
      censoring = if (!is.null(cens_balance)) setNames(cens_balance, paste0("time ", seq_len(task$time_horizon))) else NULL
    ),
    fits_outcome = lapply(estimates, function(x) x$learner_outcome_summary),
    fits_treatment = propensity_scores$learner_treatment_summary,
    fits_censoring = propensity_scores$learner_cens_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp_ltmle"
  out
}
