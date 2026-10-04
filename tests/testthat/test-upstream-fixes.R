test_that("shifts depending on past treatments read natural history", {
  dat <- data.frame(A1 = c(0, 1, 2), A2 = c(3, 4, 5))
  shift <- function(data, trt) if (trt == "A1") data$A1 + 10 else data$A2 + data$A1
  shifted <- shift_data(dat, c("A1", "A2"), NULL, shift)
  expect_equal(shifted$A1, dat$A1 + 10)
  expect_equal(shifted$A2, dat$A2 + dat$A1)
})

test_that("multivariate supplied data and policy-following are checked by row", {
  dat <- data.frame(A1 = c(0, 1, 0), B1 = c(1, 0, 0), W = 1:3)
  shifted <- dat
  shifted$A1[1] <- 1
  shifted$B1[2] <- 1
  expect_true(check_shifted_data(dat, shifted, list(c("A1", "B1")), NULL))
  expect_equal(followed_rule(dat, shifted, c("A1", "B1"), FALSE), c(0, 0, 1))
})

test_that("cluster IDs are honored by task-generated folds", {
  set.seed(844)
  dat <- data.frame(id = rep(1:12, each = 2), A = rep(0:1, 12), Y = rep(0:1, 12))
  task <- LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, "id", "binomial", 3, NULL)
  for (fold in task$folds) expect_length(intersect(dat$id[fold$training_set], dat$id[fold$validation_set]), 0)
  expect_error(LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, "id", "binomial",
                           list(list(training_set = 2:24, validation_set = 1),
                                list(training_set = 1, validation_set = 2:24)), NULL), "Cluster IDs")
})

test_that("singleton strata and supplied splits cover all rows without overlap", {
  dat <- data.frame(A = 0:5, Y = rep(0:1, 3))
  folds <- make_folds(dat, V = 3, strata_ids = 1:6)
  expect_equal(sort(unlist(lapply(folds, `[[`, "validation_set"))), 1:6)
  expect_true(all(lengths(lapply(folds, `[[`, "validation_set")) > 0))
  task <- LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, NULL, "binomial", folds, NULL)
  expect_identical(task$folds, folds)
  bad <- folds
  bad[[1]]$training_set <- seq_len(6)
  expect_error(validate_folds(bad, 6, 1:6), "complement")
})

test_that("supplied inner splits group stacked and pooled copies of each subject", {
  dat <- data.frame(id = rep(1:20, 2), A = rep(0:1, 20), Y = rep(0:1, 20))
  outer <- make_folds(dat, V = 2, cluster_ids = dat$id)
  inner <- lapply(outer, function(f) match(dat$id[f$training_set], unique(dat$id[f$training_set])) %% 2 + 1L)
  task <- LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, "id", "binomial", outer, NULL,
                       learner_folds = inner)
  ids <- rep(task$id[outer[[1]]$training_set], 3)
  labels <- learner_split(task, 1, ids, 10)$fold_id
  expect_true(all(vapply(split(labels, ids), function(x) length(unique(x)) == 1L, logical(1))))
  bad <- inner
  bad[[1]][1] <- 9
  expect_error(validate_learner_folds(bad, outer, dat$id), "same inner fold")
  expect_error(run_ensemble(dat, "Y", "glm", "binomial", "id", 30, TRUE, FALSE), "exceeds")
})

test_that("continuous bounds are used by both fixed-horizon estimators", {
  dat <- data.frame(A = rep(0:1, 10), Y = seq(2, 4, length.out = 20))
  task <- LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, NULL, "continuous", 2, NULL, bounds = c(0, 8))
  expect_equal(task$natural$Y, dat$Y / 8)
  expect_error(LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, NULL, NULL, Inf, NULL, "continuous", 2, NULL, bounds = c(0, 3)), "inside")
})

test_that("point treatment supports time-varying censoring", {
  set.seed(584)
  dat <- data.frame(A = rbinom(200, 1, .5), C1 = 1, C2 = 1, C3 = 1,
                    Y = rbinom(200, 1, .5))
  task <- LmtpTask$new(dat, dat, "A", "Y", NULL, NULL, c("C1", "C2", "C3"),
                       NULL, Inf, NULL, "binomial", 2, NULL)
  expect_equal(task$tau, 3)
  fit <- sw(lmtp_tmle(dat, "A", "Y", cens = c("C1", "C2", "C3"),
                      shift = static_binary_on, folds = 2,
                      control = lmtp_control(.learners_trt_folds = 2, .learners_outcome_folds = 2)))
  expect_equal(mean(fit$estimate@eif), fit$estimate@x, tolerance = 1e-7)
})

test_that("pivoting encodes missing observed history without lagging a terminal outcome", {
  dat <- data.frame(A1 = c(0, 1, 0), A2 = c(1, 0, 0), Y1 = c(1, NA, 0), Y2 = c(0, 1, 1), ..i..lmtp_id = 1:3)
  vars <- LmtpVars$new(NULL, NULL, c("A1", "A2"), NULL, NULL, c("Y1", "Y2"), "binomial", 2)
  long <- pivot(dat, vars)
  expect_equal(long$..i..Y_1_lag1[4:6], c(1, 0, 0))
  expect_equal(long$..i..Y_1_lag1_missing[4:6], c(0, 1, 0))
  vars <- LmtpVars$new(NULL, NULL, c("A1", "A2"), NULL, NULL, "Y2", "binomial", 2)
  expect_false(any(grepl("Y.*_lag", names(pivot(dat, vars)))))
})

test_that("EIF calculations keep the observation dimension in small validation folds", {
  ratios <- cbind(c(2, 3), c(4, 5))
  shifted <- cbind(c(.2, .3), c(.1, .2), c(1, 0))
  natural <- cbind(c(.3, .4), c(.5, .4))
  expect_equal(eif(ratios, shifted = shifted[, 1:2], natural = natural[, 1, drop = FALSE], tau = 1), c(-.2, -.3))
  expect_equal(eif(ratios, shifted = shifted, natural = natural, tau = 2), c(3.8, -6.3))
})

test_that("SuperLearner uses supplied grouped CV with an ensemble library", {
  set.seed(637)
  dat <- data.frame(id = rep(1:60, 2), W = rnorm(120))
  dat$Y <- rbinom(120, 1, plogis(dat$W))
  labels <- dat$id %% 3 + 1L
  fit <- sw(run_ensemble(dat, "Y", c("glm", "mean"), "binomial", "id",
                         list(fold_id = labels), FALSE, FALSE))
  for (j in seq_along(fit$validRows)) {
    rows <- fit$validRows[[j]]
    expect_length(unique(labels[rows]), 1)
    expect_length(intersect(dat$id[rows], dat$id[-rows]), 0)
  }
  pred <- predict(fit, dat)
  expect_true(is.numeric(pred) && length(pred) == nrow(dat))
  expect_true(all(is.finite(pred)))
})

test_that("risk-ratio inference uses the log null and positive confidence limits", {
  scores <- c(.1, .2, .5, .8)
  fit <- structure(list(estimate = ife::ife(.4, scores, rep(1, 4), as.character(1:4)),
                        outcome_type = "binomial"), class = "lmtp")
  contrast <- lmtp_contrast(fit, ref = .2, type = "rr")$estimates
  se_log <- sqrt(var(scores) / length(scores)) / .4
  expect_equal(contrast$estimate, 2)
  expect_equal(contrast$p.value, 2 * pnorm(-abs(log(2) / se_log)), tolerance = 1e-10)
  expect_equal(contrast$conf.low, exp(log(2) - qnorm(.975) * se_log), tolerance = 1e-10)
  expect_gt(contrast$conf.low, 0)
})
