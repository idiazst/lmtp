response_test_setup <- function(data, folds = NULL, learner_folds = NULL, bridge_library = NULL) {
  data$..bridge_row <- seq_len(nrow(data))
  task <- LmtpTask$new(data, make_shifted(data, "A", NULL, static_binary_on, NULL),
    "A", "Y", NULL, "W", NULL, NULL, Inf, NULL, "continuous", if (is.null(folds)) 3 else folds,
    NULL, c(0, 6), learner_folds)
  prepare_curve_splits(task, lmtp_control(.learners_trt_folds = 2, .learners_outcome_folds = 2))
  train <- task$folds[[1]]$training_set
  split <- make_bridge_nested_folds(task$id[train], data[train, "R", drop = FALSE], task$learner_folds[[1]])
  lib <- list(cv = list(method = "saturated_l1", control = list(
    penalty_scales = 10^seq(-5, 0, by = .25))))
  if (!is.null(bridge_library)) lib <- bridge_library
  engine <- bridge_response_engine(data, data, task, 1, "R", list("Y"), lib,
    c("SL.saturated_l1_cv", "SL.mean"), c("SL.saturated_l1_cv", "SL.mean"), FALSE,
    lmtp_control(.learners_trt_folds = 2, .learners_outcome_folds = 2, .discrete = FALSE), "sdr", split)
  list(task = task, engine = engine, data = data, train = train)
}

test_that("sequential-response sieve CV matches cmbridge and excludes outer validation outcomes", {
  set.seed(948)
  n <- 900
  data <- data.frame(W = rbinom(n, 1, .5), A = rbinom(n, 1, .5))
  data$Y <- rbinom(n, 1, .15 + .7*data$A)
  data$R <- rbinom(n, 1, .3 + .3*data$Y); data$Y[data$R == 0] <- NA_real_
  lib <- list(sieve = list(method = "sieve_md", control = list(target_basis = "cell",
    instrument_basis = "cell", link = "inverse_logit", penalty_scales = 10^seq(-10, 0, by = .5))))
  before <- response_test_setup(data, bridge_library = lib)
  rows <- before$train; task <- before$task
  model <- get("fit_beta", environment(before$engine$response))(rows, 1L)
  H <- task$vars$history("A", 1L)
  B <- as.matrix(task$natural[rows, c(H, "A"), drop = FALSE])
  V <- as.matrix(cbind(task$natural[rows, H, drop = FALSE], data[rows, "Y", drop = FALSE]))
  lib$sieve$control$penalty_ids <- task$id[rows]
  lib$sieve$control$penalty_folds <- before$engine$splits
  reference <- cmbridge::fit_bridge_ensemble(B, V, data$R[rows], lib,
    fold_id = before$engine$splits(task$id[rows]), kernel = "cell")$candidates$sieve
  expect_identical(model$method, "sieve_md")
  # The single-candidate shortcut is always U; cmbridge also records the
  # statistic column introduced with its opt-in V interface. Compare every
  # tuning value while explicitly verifying this extra reference metadata.
  expect_true(all(reference$penalty_cv$statistic == "u"))
  reference_cv <- reference$penalty_cv[, setdiff(names(reference$penalty_cv), "statistic"), drop = FALSE]
  expect_equal(model$penalty_cv, reference_cv, tolerance = 1e-12)
  expect_equal(model$tuning$lambda, reference$tuning$lambda)
  expect_equal(model$penalty_fold_id, before$engine$splits(model$penalty_ids))
  expect_identical(model$penalty_ids, task$id[rows])
  validation <- task$folds[[1]]$validation_set
  data$Y[validation] <- ifelse(data$R[validation] == 1, 1 - data$Y[validation], NA_real_)
  after <- response_test_setup(data, task$folds, task$learner_folds, lib)
  changed <- get("fit_beta", environment(after$engine$response))(rows, 1L)
  expect_equal(model$penalty_cv, changed$penalty_cv)
  expect_equal(model$coefficients, changed$coefficients)
})

test_that("upstream fits for candidate training responses exclude its learner validation outcomes", {
  set.seed(221)
  n <- 800
  data <- data.frame(W = rbinom(n, 1, .5), A = rbinom(n, 1, .5))
  data$Y <- rbinom(n, 1, plogis(-2 + 4 * data$A))
  data$R <- rbinom(n, 1, .4 + .2 * data$Y); data$Y[data$R == 0] <- NA_real_
  before <- response_test_setup(data)
  labels <- before$task$learner_folds[[1]]
  val <- before$train[labels == 1]; train <- before$train[labels != 1]
  one <- before$engine$response(train, rep(1, length(train)), val, rep(1, length(val)), 1)
  data$Y[val] <- ifelse(data$R[val] == 1, 1 - data$Y[val], NA_real_)
  after <- response_test_setup(data, before$task$folds, before$task$learner_folds)
  two <- after$engine$response(train, rep(1, length(train)), val, rep(1, length(val)), 1)
  expect_equal(one$train, two$train, tolerance = 1e-12)
  expect_gt(max(abs(one$validation - two$validation)), .1)
  make_fit <- function(setup) {
    z <- data.frame(W = setup$data$W[setup$train], A = setup$data$A[setup$train],
      id = setup$train, Y = 0)
    rebuild <- function(ids, times, valid_ids, valid_times)
      setup$engine$response(ids, times, valid_ids, valid_times, 1)
    run_ensemble(z, "Y", c("SL.saturated_l1_cv", "SL.mean"), "continuous", "id",
      list(fold_id = labels), FALSE, FALSE, rebuild, setup$engine$splits)
  }
  fit_one <- make_fit(before); fit_two <- make_fit(after)
  expect_false(any(fit_one$errorsInCVLibrary) || any(fit_two$errorsInCVLibrary))
  expect_equal(fit_one$Z[labels == 1, ], fit_two$Z[labels == 1, ], tolerance = 1e-12)
})

test_that("backward responses exclude learner validation outcomes for SDR and TMLE", {
  set.seed(224)
  n <- 1000
  data <- data.frame(W = rbinom(n, 1, .5), A1 = rbinom(n, 1, .5))
  data$Y2 <- rbinom(n, 1, plogis(-1 + 2 * data$A1))
  data$A2 <- rbinom(n, 1, plogis(-.5 + data$Y2))
  data$Y3 <- rbinom(n, 1, plogis(-1 + 2 * data$A2))
  data$Y1 <- 0; data$R2 <- 1; data$R3 <- rbinom(n, 1, .5)
  data$Y3[data$R3 == 0] <- NA_real_
  setup <- function(data, folds = 3, learner_folds = NULL, estimator = "sdr") {
    data$..bridge_row <- seq_len(nrow(data))
    data$C1 <- 0; data$C2 <- data$Y2
    task <- LmtpTask$new(data,
      make_shifted(data, c("A1", "A2"), NULL, static_binary_on, NULL),
      c("A1", "A2"), c("Y2", "Y3"), list("C1", "C2"), "W",
      NULL, NULL, Inf, NULL, "continuous", folds, NULL, c(0, 6), learner_folds)
    ctrl <- lmtp_control(.learners_trt_folds = 2, .learners_outcome_folds = 2, .discrete = FALSE)
    prepare_curve_splits(task, ctrl)
    train <- task$folds[[1]]$training_set
    split <- make_bridge_nested_folds(task$id[train], data[train, c("R2", "R3")], task$learner_folds[[1]])
    lib <- list(cv = list(method = "saturated_l1", control = list(
      penalty_scales = 10^seq(-5, 0, by = .25))))
    engine <- bridge_response_engine(data, data, task, 1, c("R2", "R3"),
      list("Y2", "Y3"), lib, c("SL.saturated_l1_cv", "SL.mean"),
      c("SL.saturated_l1_cv", "SL.mean"), FALSE, ctrl, estimator, split)
    list(task = task, engine = engine, train = train)
  }
  initial <- setup(data)
  labels <- initial$task$learner_folds[[1]]
  val <- initial$train[labels == 1]; train <- initial$train[labels != 1]
  changed <- data
  changed$Y2[val] <- 1 - changed$Y2[val]
  changed$Y3[val] <- ifelse(changed$R3[val] == 1, 1 - changed$Y3[val], NA_real_)
  for (method in c("sdr", "tmle")) {
    before <- setup(data, initial$task$folds, initial$task$learner_folds, method)
    after <- setup(changed, initial$task$folds, initial$task$learner_folds, method)
    one <- before$engine$response(train, rep(1, length(train)), val, rep(1, length(val)), 2)
    two <- after$engine$response(train, rep(1, length(train)), val, rep(1, length(val)), 2)
    expect_equal(one$train, two$train, tolerance = 1e-12)
    expect_true(all(is.finite(c(one$train, one$validation))))
  }
})
test_that("additional bridge validation preserves measured observations in small training sets", {
  ids <- seq_len(90)
  measured <- cbind(R2 = 1, R3 = rep(0:1, 45))
  labels <- rep(1:3, 30)
  split <- make_bridge_nested_folds(ids, measured, labels)
  subset <- 1:24
  assigned <- split(subset)
  expect_setequal(unique(assigned), 1:3)
  expect_true(all(vapply(1:3, function(group) sum(measured[subset[assigned == group], "R3"]) > 0, TRUE)))
  duplicated <- rep(subset, 2)
  expect_equal(split(duplicated), rep(assigned, 2))
  expect_setequal(unique(split(1:60)), 1:2)
})

test_that("SDR and TMLE reuse bridge fits only with identical data and splits", {
  set.seed(2235)
  n <- 240L
  data <- data.frame(W=rbinom(n,1,.5), A=rbinom(n,1,.5))
  data$Y <- rbinom(n,1,.2+.6*data$A)
  data$R <- rbinom(n,1,.35+.3*data$Y); data$Y[data$R==0]<-NA_real_
  lib <- list(sieve=list(method="sieve_md",control=list(target_basis="cell",instrument_basis="cell",link="inverse_logit")))
  setup <- response_test_setup(data,bridge_library=lib)
  cache <- new.env(parent=emptyenv()); split <- setup$engine$splits
  make_engine <- function(dat, estimator) bridge_response_engine(dat,dat,setup$task,1,"R",list("Y"),lib,
    "SL.mean","SL.mean",FALSE,lmtp_control(.discrete=FALSE),estimator,split,
    list(response_nuisance_cache=cache))
  first <- make_engine(setup$data,"sdr")
  beta <- get("fit_beta",environment(first$response))(setup$train,1L)
  expect_length(ls(cache$beta),1L)
  second <- make_engine(setup$data,"tmle")
  reused <- get("fit_beta",environment(second$response))(setup$train,1L)
  expect_identical(beta,reused)
  expect_length(ls(cache$beta),1L)
  # Changing only outer validation data cannot affect the cache signature.
  changed <- setup$data; validation <- setup$task$folds[[1]]$validation_set
  changed$Y[validation] <- ifelse(changed$R[validation]==1,1-changed$Y[validation],NA_real_)
  third <- make_engine(changed,"sdr")
  expect_length(ls(cache$beta),1L)
  expect_identical(get("fit_beta",environment(third$response))(setup$train,1L),beta)
  # A change to training data invalidates that cache before another fit.
  changed$Y[setup$train] <- ifelse(changed$R[setup$train]==1,1-changed$Y[setup$train],NA_real_)
  fourth <- make_engine(changed,"tmle")
  expect_length(ls(cache$beta),0L)
})
