test_that("fold-specific bridge outcomes and correction use the raw scale", {
  set.seed(812)
  n <- 120; dat <- data.frame(W=rbinom(n,1,.5),A=rbinom(n,1,.5),Y=rbinom(n,1,.4),R=rbinom(n,1,.5))
  dat$Y[dat$R==0] <- NA_real_
  supplied <- function(training_data,validation_data,inner_folds,task,fold) {
    make <- function(n,x) matrix(x,n,1)
    bt <- make(nrow(training_data),2);bv <- make(nrow(validation_data),2)
    bt[training_data$R==0] <- NA_real_;bv[validation_data$R==0] <- NA_real_
    list(beta_training=bt,beta_validation=bv,
      lambda_validation=make(nrow(validation_data),.7),ratios_training=make(nrow(training_data),1),
      ratios_validation=make(nrow(validation_data),1),
      oracle_regressions=list(training=list(natural=list(make(nrow(training_data),.4)),shifted=list(make(nrow(training_data),.4))),
        validation=list(natural=list(make(nrow(validation_data),.4)),shifted=list(make(nrow(validation_data),.4)))))
  }
  f <- lmtp_bridge(dat,"A","Y","R",list("Y"),baseline="W",shift=static_binary_on,
      mtp=FALSE,bounds=c(0,6),nuisance_fit=supplied,
      control=lmtp_control(.learners_trt_folds=2,.learners_outcome_folds=2))
  expected <- ifelse(dat$R==1,2*dat$Y,0)-.7*(2*dat$R-1)
  expect_equal(f$estimates[[1]]@x,mean(expected),tolerance=1e-12)
  expect_equal(f$estimates[[1]]@eif,expected,tolerance=1e-12)
  expect_equal(f$estimates[[1]]@std_error,sd(expected)/sqrt(n),tolerance=1e-12)
})

test_that("saturated learner includes joint interactions and unseen-cell fallback", {
  X <- expand.grid(x=0:1,z=0:1)
  X <- X[rep(1:4,each=100),]
  y <- as.numeric(X$x==X$z)
  f <- SL.saturated_l1(y,X,X,gaussian(),rep(1,length(y)),seq_along(y),penalty_scale=1e-6)
  expect_lt(max(abs(f$pred-y)),1e-5)
  expect_equal(predict(f$fit,data.frame(x=2,z=2)),mean(y))
})

test_that("fast saturated keys preserve factor and numeric coding", {
  x <- data.frame(a = rep(c(0, 1 / 6, 1, NA_real_), 100), b = factor(rep(c("off", "on"), 200)))
  expect_identical(.saturated_keys(x), do.call(paste, c(x, sep = "|")))
})

test_that("penalty selection excludes newX and keeps repeated person IDs together", {
  set.seed(440)
  id <- rep(seq_len(200), 2)
  X <- data.frame(x = rep(rbinom(200, 1, .5), 2), z = rep(0:1, each = 200))
  y <- rbinom(400, 1, plogis(-1 + 2 * X$x))
  one <- SL.saturated_l1_cv(y, X, X, binomial(), rep(1, 400), id)
  two <- SL.saturated_l1_cv(y, X, data.frame(x = 2, z = 2), binomial(), rep(1, 400), id)
  expect_equal(one$fit$penalty_cv, two$fit$penalty_cv)
  expect_equal(sum(one$fit$penalty_cv$selected), 1L)
  expect_gt(one$fit$penalty, 0)
  even <- id %% 2 == 0
  expect_silent(SL.saturated_l1_cv(y[even], X[even, ], X, binomial(),
                                 rep(1, sum(even)), id[even]))
  expect_error(SL.saturated_l1_cv(y, X, X, binomial(), rep(1, 400), id,
                                 penalty_scales = 0), "positive")
  # Calendar-time constants preserve the separate endpoint regressions.
  mean_fit <- SL.time_mean(c(1, 3, 10, 20), data.frame(time = c(1, 1, 2, 2)),
    data.frame(time = c(1, 2)), gaussian(), rep(1, 4), 1:4)
  expect_equal(mean_fit$pred, c(2, 15))
})

test_that("cell aggregation preserves the Gaussian and binomial penalized objectives", {
  set.seed(568)
  n <- 800
  X <- data.frame(x=rbinom(n,1,.5),z=rbinom(n,1,.5))
  keys <- do.call(paste,c(X,sep="|"));j <- match(keys,sort(unique(keys)))
  design <- cbind(Matrix::sparseMatrix(i=seq_len(n),j=j,x=1,dims=c(n,4)),0)
  weights <- runif(n,.5,2)
  for(fam in c("gaussian","binomial")) {
    y <- if(fam=="gaussian")rnorm(n)+X$x else rbinom(n,1,plogis(X$x))
    fit <- SL.saturated_l1(y,X,X,get(fam)(),weights,seq_len(n))
    reference <- glmnet::glmnet(design,y,weights=weights,family=fam,alpha=1,
      standardize=FALSE,lambda=.01/sqrt(n),thresh=1e-9,maxit=100000)
    expected <- as.numeric(predict(reference,design,s=.01/sqrt(n),type="response"))
    expect_equal(fit$pred,expected,tolerance=1e-6)
  }
})

test_that("outer validation outcomes cannot alter its fitted bridge nuisances", {
  set.seed(845)
  n <- 600
  dat <- data.frame(W=rbinom(n,1,.5),A=rbinom(n,1,.5))
  dat$Y <- rbinom(n,1,plogis(-2+4*dat$A))
  dat$R <- rbinom(n,1,.4+.2*dat$Y); dat$Y[dat$R==0] <- NA_real_
  folds <- make_folds(dat,V=3)
  args <- list(trt="A",outcome="Y",measurement="R",health=list("Y"),baseline="W",
    shift=static_binary_on,mtp=FALSE,folds=folds,
    control=lmtp_control(.learners_trt_folds=2,.learners_outcome_folds=2))
  set.seed(456); one <- do.call(lmtp_bridge,c(list(data=dat),args))
  for(fold in one$bridge_fits) for(model in c(fold$beta_fits, fold$adjoint_fits)) {
    expect_setequal(names(model$weights), c("sieve_md", "landweber", "pmmr"))
    expect_identical(model$candidates$sieve_md$target_spec$type, "cell")
    for(candidate in model$candidates) {
    expect_false(identical(candidate$method, "saturated_l1"))
    expect_true(is.null(candidate$tuning$lower) || identical(candidate$tuning$lower, -Inf))
    expect_true(is.null(candidate$tuning$upper) || identical(candidate$tuning$upper, Inf))
    }
  }
  for (fold in one$bridge_fits) {
    for (model in c(fold$beta_fits, fold$adjoint_fits)) {
      sieve <- model$candidates$sieve_md
      selected <- which(sieve$penalty_cv$selected)
      expect_length(selected, 1L)
      expect_true(selected > 1 && selected < nrow(sieve$penalty_cv))
      expect_equal(sieve$tuning$lambda, sieve$penalty_cv$scale[selected])
      expect_gt(sieve$tuning$lambda, 0)
      expect_equal(sieve$penalty_fold_id, fold$nested_split_function(sieve$penalty_ids))
    }
    for (model in fold$beta_fits) {
      for (candidate in model$candidates) {
        expect_identical(candidate$tuning$link, "inverse_logit")
        expect_true(all(fitted(candidate)[!is.na(fitted(candidate))] >= 1))
        expect_true(all(is.finite(candidate$link_coefficients)))
      }
      expect_true(all(fitted(model)[!is.na(fitted(model))] >= 1))
    }
    for (model in fold$adjoint_fits) for (candidate in model$candidates)
      expect_identical(candidate$tuning$link, "identity")
  }
  v <- folds[[1]]$validation_set
  dat$Y[v] <- ifelse(dat$R[v]==1,1-dat$Y[v],NA_real_)
  set.seed(456); two <- do.call(lmtp_bridge,c(list(data=dat),args))
  expect_equal(one$bridge_fits[[1]]$beta_fits[[1]]$weights,two$bridge_fits[[1]]$beta_fits[[1]]$weights)
  expect_equal(one$bridge_fits[[1]]$beta_fits[[1]]$candidates$sieve_md$penalty_cv,
               two$bridge_fits[[1]]$beta_fits[[1]]$candidates$sieve_md$penalty_cv)
  for(k in seq_along(one$bridge_fits[[1]]$beta_fits[[1]]$candidates)) {
    expect_equal(one$bridge_fits[[1]]$beta_fits[[1]]$candidates[[k]]$coefficients,
                 two$bridge_fits[[1]]$beta_fits[[1]]$candidates[[k]]$coefficients)
  }
  expect_equal(one$bridge_fits[[1]]$ratios_validation,two$bridge_fits[[1]]$ratios_validation)
  expect_equal(one$bridge_fits[[1]]$adjoint_fits[[1]]$weights,two$bridge_fits[[1]]$adjoint_fits[[1]]$weights)
})
