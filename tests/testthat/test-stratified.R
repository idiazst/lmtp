context("Stratified LMTP estimates")

test_that("lmtp_tmle exposes a strata argument", {
  expect_true("strata" %in% names(formals(lmtp::lmtp_tmle)))
})

make_mock_task <- function(weights = rep(1, 8)) {
  list(
    weights = weights,
    id = seq_along(weights),
    rescale = function(x) x
  )
}

test_that("multiple stratification variables produce the requested estimates", {
  d <- c(1, 3, 5, 7, 2, 4, 6, 8)
  strata <- data.frame(
    sex = rep(c("F", "M"), each = 4),
    age_group = rep(c("younger", "older"), times = 4)
  )
  task <- make_mock_task()
  estimates <- list(
    predictions = cbind(d),
    uncentered_eif = d
  )

  out <- lmtp:::stratified_lmtp_estimates(task, estimates, strata, is_sdr = FALSE)

  expect_length(out$estimates, 4)
  expect_equal(names(out$estimates), c(
    "sex=F, age_group=younger",
    "sex=F, age_group=older",
    "sex=M, age_group=younger",
    "sex=M, age_group=older"
  ))
  expect_equal(unname(vapply(out$estimates, function(x) x@x, numeric(1))), c(3, 5, 4, 6))
  expect_equal(out$table[["..lmtp_probability.."]], rep(0.25, 4))
  expect_equal(out$table[["..lmtp_n.."]], rep(2L, 4))

  in_first <- strata$sex == "F" & strata$age_group == "younger"
  expected_eif <- as.numeric(in_first) / mean(in_first) * (d - mean(d[in_first]))
  expect_equal(out$estimates[[1]]@eif, expected_eif)
  expect_equal(out$estimates[[1]]@std_error, sqrt(var(expected_eif) / length(d)))
  expect_equal(
    out$estimates[[1]]@conf_int,
    out$estimates[[1]]@x + c(-1, 1) * qnorm(0.975) * out$estimates[[1]]@std_error
  )
})

test_that("stratum probabilities and means respect survey weights", {
  d <- c(1, 3, 5, 7, 2, 4, 6, 8)
  weights <- c(0.5, 1.5, 1, 1, 0.5, 1.5, 1, 1)
  strata <- data.frame(group = rep(c("A", "B"), each = 4))
  task <- make_mock_task(weights)
  estimates <- list(
    predictions = cbind(d),
    uncentered_eif = d
  )

  out <- lmtp:::stratified_lmtp_estimates(task, estimates, strata, is_sdr = FALSE)
  in_a <- strata$group == "A"
  p_a <- weighted.mean(as.numeric(in_a), weights)
  theta_a <- weighted.mean(d[in_a], weights[in_a])
  expected_eif <- as.numeric(in_a) / p_a * (d - theta_a)

  expect_equal(out$estimates[["group=A"]]@x, theta_a)
  expect_equal(out$table["group=A", "..lmtp_probability.."], p_a)
  expect_equal(out$estimates[["group=A"]]@eif, expected_eif)
})

test_that("missing values are retained as an explicit stratum", {
  d <- 1:4
  strata <- data.frame(group = c("A", NA, "A", NA))
  task <- list(weights = rep(1, 4), id = 1:4, rescale = function(x) x)
  estimates <- list(predictions = cbind(d), uncentered_eif = d)

  out <- lmtp:::stratified_lmtp_estimates(task, estimates, strata, is_sdr = FALSE)

  expect_true("group=<NA>" %in% names(out$estimates))
  expect_equal(out$estimates[["group=<NA>"]]@x, mean(c(2, 4)))
})

test_that("tidy includes population and stratified inference", {
  d <- c(1, 3, 5, 7)
  strata <- data.frame(group = rep(c("A", "B"), each = 2))
  task <- list(weights = rep(1, 4), id = 1:4, rescale = function(x) x)
  estimates <- list(predictions = cbind(d), uncentered_eif = d)
  stratified <- lmtp:::stratified_lmtp_estimates(task, estimates, strata, is_sdr = FALSE)

  fit <- list(
    estimator = "TMLE",
    estimate = ife::ife(mean(d), d),
    stratified_estimates = stratified$estimates,
    strata_table = stratified$table,
    shift = "NULL"
  )
  class(fit) <- "lmtp"

  out <- generics::tidy(fit)
  expect_equal(out$stratum, c("Population", "group=A", "group=B"))
  expect_equal(out$probability, c(1, 0.5, 0.5))
  expect_equal(out$n, c(4L, 2L, 2L))
  expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(out)))
})
