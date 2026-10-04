theta_dr <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, augmented = FALSE) {
  ic <- eif(density_ratios, shifted = sequential_regressions$shifted, natural = sequential_regressions$natural)

  if (augmented) {
    theta <- weighted.mean(ic, task$weights)
  } else {
    theta <- weighted.mean(sequential_regressions$shifted[, 1], task$weights)
  }

  ic <- task$rescale(ic)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(augmented, "SDR", "TMLE"),
    folds = task$folds,
    learner_folds = task$learner_folds,
    estimate = ife::ife(theta, ic, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(sequential_regressions$shifted),
    density_ratios = density_ratios,
    fits_m = fits_m,
    fits_r = fits_r,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}

theta_curve <- function(task, influence_functions, sequential_regressions, density_ratios, sporadic_weights, fits_m, fits_r, fits_sporadic, shift,
                        estimator = "sdr", targeting = NULL) {
  ics <- as.list(as.data.frame(influence_functions))

  thetas <- unlist(lapply(ics, \(x) weighted.mean(x, task$weights)))
  if (estimator == "tmle") {
    thetas <- vapply(sequential_regressions, function(x) weighted.mean(x[, 1], task$weights), numeric(1))
  }

  ics <- lapply(ics, \(x) task$rescale(x))
  thetas <- task$rescale(thetas)
  unprojected <- thetas

  if (task$survival) {
    iso_projection <- isotone::gpava(seq_along(thetas), 1 - thetas)
    thetas <- 1 - iso_projection$x[seq_along(thetas)]
  }

  out <- list(
    estimator = paste(toupper(estimator), "curve"),
    folds = task$folds,
    learner_folds = task$learner_folds,
    targeting = targeting,
    unprojected_estimates = unprojected,
    estimates = lapply(
      seq_along(ics),
      \(x) ife::ife(thetas[[x]], ics[[x]], task$weights, as.character(task$id))
    ),
    outcome_reg = lapply(sequential_regressions, \(x) task$rescale(x)),
    density_ratios = density_ratios,
    sporadic_weights = sporadic_weights,
    fits_m = fits_m,
    fits_r = fits_r,
    fits_sporadic = fits_sporadic,
    shift = shift,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp_curve"
  out
}
