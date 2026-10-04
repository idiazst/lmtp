curve_fixture <- function(n = 600, seed = 712) {
  set.seed(seed)
  W <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(.2 * W))
  A2 <- rbinom(n, 1, plogis(-.2 + .3 * W + .2 * A1))
  A3 <- rbinom(n, 1, plogis(.1 + .2 * W + .2 * A2))
  data.frame(W, A1, A2, A3,
             Y1 = rbinom(n, 1, plogis(-1 + .6 * A1 + .2 * W)),
             Y2 = rbinom(n, 1, plogis(-.5 + .5 * A2 + .2 * W)),
             Y3 = rbinom(n, 1, plogis(.1 + .4 * A3 + .2 * W)))
}

curve_control <- function(...) lmtp_control(.learners_trt_folds = 2,
                                           .learners_outcome_folds = 2, ...)

test_that("pooled TMLE solves each held-out score and retains pooling", {
  dat <- curve_fixture()
  dat$weight <- exp(.3 * dat$W)
  fit <- sw(lmtp_curve(dat, paste0("A", 1:3), paste0("Y", 1:3), baseline = "W",
                      weights = dat$weight, shift = static_binary_on, folds = 3,
                      estimator = "tmle", control = curve_control(.return_full_fits = TRUE)))
  expect_identical(fit$estimator, "TMLE curve")
  for (fold in fit$targeting) for (point in fold) expect_lt(max(abs(point$score)), 1e-9)
  for (s in 1:3) {
    expect_equal(weighted.mean(fit$estimates[[s]]@eif, dat$weight), fit$estimates[[s]]@x, tolerance = 1e-9)
    expect_true(all(fit$outcome_reg[[s]][, 1] >= 0 & fit$outcome_reg[[s]][, 1] <= 1))
  }
  # One pooled initial fit per depth, not one full fit per outcome horizon.
  expect_length(fit$fits_m[[1]], 3)
  expect_true("time" %in% fit$fits_m[[1]][[1]]$varNames)
  expect_false("time" %in% fit$fits_m[[1]][[3]]$varNames)
})

test_that("a one-point pooled TMLE agrees with ordinary held-out TMLE", {
  dat <- curve_fixture()
  splits <- make_folds(dat, V = 3)
  args <- list(data = dat, trt = "A1", outcome = "Y1", baseline = "W",
               shift = static_binary_on, folds = splits, control = curve_control())
  pooled <- sw(do.call(lmtp_curve, c(args, list(estimator = "tmle"))))
  ordinary <- sw(do.call(lmtp_tmle, args))
  expect_equal(pooled$estimates[[1]]@x, ordinary$estimate@x, tolerance = 1e-7)
  expect_equal(pooled$estimates[[1]]@eif, ordinary$estimate@eif, tolerance = 1e-7)
})

test_that("both pooled estimators recover known randomized-treatment curves", {
  set.seed(8712)
  n <- 4000
  dat <- as.data.frame(matrix(rbinom(n * 3, 1, .5), n, 3))
  names(dat) <- paste0("A", 1:3)
  for (s in 1:3) dat[[paste0("Y", s)]] <- rbinom(n, 1, plogis(-1.5 + .5 * s + .4 * dat[[s]]))
  truth <- plogis(-1.5 + .5 * (1:3) + .4)
  for (method in c("sdr", "tmle")) {
    fit <- sw(lmtp_curve(dat, paste0("A", 1:3), paste0("Y", 1:3),
                        shift = static_binary_on, estimator = method, folds = 3,
                        control = curve_control(.isotonic_constraint = FALSE)))
    expect_lt(max(abs(vapply(fit$estimates, function(x) x@x, numeric(1)) - truth)), .035)
  }
})

test_that("TMLE handles sporadic measurement including the final outcome", {
  dat <- curve_fixture(1200)
  set.seed(847)
  for (s in 1:3) dat[[paste0("Y", s)]][runif(nrow(dat)) < plogis(-2 + .3 * dat$W)] <- NA_real_
  fit <- sw(lmtp_curve(dat, paste0("A", 1:3), paste0("Y", 1:3), baseline = "W",
                      shift = static_binary_on, estimator = "tmle", folds = 3,
                      control = curve_control()))
  expect_true(all(is.finite(tidy(fit)$estimate)))
  expect_true(all(is.finite(unlist(lapply(fit$estimates, function(x) x@eif)))))
  expect_true(all(fit$sporadic_weights[is.na(dat$Y3), 3] == 0))
  for (fold in fit$targeting) for (point in fold) expect_lt(max(abs(point$score)), 1e-9)
})

test_that("continuous TMLE respects supplied bounds and separated folds", {
  dat <- curve_fixture()
  for (s in 1:3) dat[[paste0("Y", s)]] <- 2 + dat[[paste0("Y", s)]] * 4
  fit <- sw(lmtp_curve(dat, paste0("A", 1:3), paste0("Y", 1:3), baseline = "W",
                      shift = static_binary_on, estimator = "tmle", folds = 3,
                      outcome_type = "continuous", bounds = c(0, 8), control = curve_control()))
  expect_true(all(tidy(fit)$estimate >= 0 & tidy(fit)$estimate <= 8))
  for (point in fit$estimates) expect_equal(mean(point@eif), point@x, tolerance = 1e-9)
  expect_identical(curve_fluctuation(rep(0, 3), rep(.5, 3), rep(1, 3), 1e-5), -Inf)
  expect_identical(curve_fluctuation(rep(1, 3), rep(.5, 3), rep(1, 3), 1e-5), Inf)
})

test_that("outer validation changes cannot change that fold's initial nuisance fits", {
  dat <- curve_fixture(400)
  splits <- make_folds(dat, V = 2)
  inner <- lapply(splits, function(f) rep(1:2, length.out = length(f$training_set)))
  changed <- dat
  valid <- splits[[1]]$validation_set
  changed[valid, paste0("Y", 1:3)] <- 1 - changed[valid, paste0("Y", 1:3)]
  changed$W[valid] <- changed$W[valid] + 2
  args <- list(trt = paste0("A", 1:3), outcome = paste0("Y", 1:3), baseline = "W",
               shift = static_binary_on, folds = splits, learner_folds = inner,
               control = curve_control(.return_full_fits = TRUE, .isotonic_constraint = FALSE))
  coefficients <- function(fit) lapply(fit$fits_m[[1]], function(x) coef(x$fitLibrary[[1]]$object))
  for (method in c("sdr", "tmle")) {
    set.seed(879)
    before <- sw(do.call(lmtp_curve, c(list(data = dat, estimator = method), args)))
    set.seed(879)
    after <- sw(do.call(lmtp_curve, c(list(data = changed, estimator = method), args)))
    expect_equal(coefficients(before), coefficients(after), tolerance = 1e-10)
  }
})

test_that("pooled TMLE handles censored survival curves", {
  dat <- sim_point_surv[1:800, ]
  fit <- sw(lmtp_curve(dat, "trt", paste0("Y.", 1:6), baseline = c("W1", "W2"),
                      cens = paste0("C.", 0:5), shift = static_binary_on,
                      outcome_type = "survival", estimator = "tmle", folds = 3,
                      control = curve_control()))
  expect_true(all(is.finite(tidy(fit)$estimate)))
  expect_true(all(diff(tidy(fit)$estimate) <= 0))
  for (fold in fit$targeting) for (point in fold) expect_lt(max(abs(point$score)), 1e-9)
})

test_that("multivariate curve TMLE uses supplied shifted data", {
  dat <- curve_fixture()
  set.seed(291)
  dat$B1 <- rbinom(nrow(dat), 1, .5)
  dat$B2 <- rbinom(nrow(dat), 1, .5)
  shifted <- dat
  shifted[c("A1", "B1", "A2", "B2")] <- 1
  fit <- sw(lmtp_curve(dat, list(c("A1", "B1"), c("A2", "B2")), c("Y1", "Y2"),
                      baseline = "W", shifted = shifted, mtp = FALSE,
                      estimator = "tmle", folds = 3, control = curve_control()))
  expect_true(all(is.finite(tidy(fit)$estimate)))
  for (fold in fit$targeting) for (point in fold) expect_lt(max(abs(point$score)), 1e-9)
})

test_that("a point-treatment MTP is applied once at every outcome horizon", {
  set.seed(141)
  n <- 1600
  dat <- data.frame(A = runif(n))
  for (s in 1:3) dat[[paste0("Y", s)]] <- 1 + 2 * dat$A + rnorm(n, 0, .1)
  policy <- function(data, trt) .8 * data[[trt]]
  truth <- 1 + 1.6 * mean(dat$A)
  for (method in c("sdr", "tmle")) {
    fit <- sw(lmtp_curve(dat, "A", paste0("Y", 1:3), shift = policy,
                        outcome_type = "continuous", bounds = c(0, 4), folds = 2,
                        estimator = method, control = curve_control()))
    expect_lt(max(abs(tidy(fit)$estimate - truth)), .025)
  }
  # The fixed-horizon path must also avoid shifting point treatment repeatedly.
  dat[c("C1", "C2", "C3")] <- 1
  for (method in c("lmtp_sdr", "lmtp_tmle")) {
    fit <- sw(do.call(method, list(data = dat, trt = "A", outcome = "Y3",
                      cens = c("C1", "C2", "C3"), shift = policy,
                      outcome_type = "continuous", bounds = c(0, 4), folds = 2,
                      control = curve_control())))
    expect_lt(abs(fit$estimate@x - truth), .025)
  }
})
