#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a(n) lmtp object
#'
#' @param x A `lmtp` object produced by a call to [lmtp::lmtp_tmle()], [lmtp::lmtp_sdr()],
#' [lmtp::lmtp_survival()].
#' @param ... Unused, included for generic consistency only.
#'
#' @examples
#' \donttest{
#' a <- c("A1", "A2")
#' nodes <- list(c("L1"), c("L2"))
#' cens <- c("C1", "C2")
#' y <- "Y"
#' fit <- lmtp_tmle(sim_cens, a, y, time_vary = nodes, cens = cens, shift = NULL, folds = 2)
#' tidy(fit)
#' }
#'
#' @export
tidy.lmtp <- function(x, ...) {
  population <- ife::tidy(x$estimate)
  if (is.null(x$stratified_estimates)) {
    return(population)
  }

  population$stratum <- "Population"
  population$probability <- 1
  population$n <- sum(x$strata_table[["..lmtp_n.."]])

  stratified <- do.call("rbind", lapply(x$stratified_estimates, ife::tidy))
  stratified$stratum <- names(x$stratified_estimates)
  stratified$probability <- x$strata_table[["..lmtp_probability.."]]
  stratified$n <- x$strata_table[["..lmtp_n.."]]

  out <- rbind(population, stratified)
  rownames(out) <- NULL
  out[, c("stratum", "probability", "n", "estimate", "std.error", "conf.low", "conf.high")]
}

#' Tidy a(n) lmtp_survival object
#'
#' @param x A `lmtp_survival` object produced by a call to [lmtp::lmtp_survival()].
#' @param ... Unused, included for generic consistency only.
#'
#' @example inst/examples/lmtp_survival-ex.R
#'
#' @export
tidy.lmtp_survival <- function(x, ...) {
  out <- do.call("rbind", lapply(x, tidy.lmtp))
  out$time <- seq_along(x)
  out[, c(ncol(out), 1:ncol(out) - 1)]
}

#' Tidy a(n) lmtp_ltmle object
#'
#' @param x A `lmtp_ltmle` object produced by a call to [lmtp::ltmle()].
#' @param ... Unused, included for generic consistency only.
#'
#' @example
#'
#' @export
tidy.lmtp_ltmle <- function(x, ...) {
  out <- do.call("rbind", lapply(x$estimates, ife::tidy))
  out$level <- names(x$estimates)
  out[, c(ncol(out), 1:ncol(out) - 1)]
}
