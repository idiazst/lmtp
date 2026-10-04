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
  for(fold in one$bridge_fits) for(candidate in fold$adjoint_fits[[1]]$candidates) {
    expect_equal(candidate$tuning$lower,-Inf)
    expect_equal(candidate$tuning$upper,Inf)
  }
  v <- folds[[1]]$validation_set
  dat$Y[v] <- ifelse(dat$R[v]==1,1-dat$Y[v],NA_real_)
  set.seed(456); two <- do.call(lmtp_bridge,c(list(data=dat),args))
  expect_equal(one$bridge_fits[[1]]$beta_fits[[1]]$weights,two$bridge_fits[[1]]$beta_fits[[1]]$weights)
  for(k in seq_along(one$bridge_fits[[1]]$beta_fits[[1]]$candidates)) {
    expect_equal(one$bridge_fits[[1]]$beta_fits[[1]]$candidates[[k]]$coefficients,
                 two$bridge_fits[[1]]$beta_fits[[1]]$candidates[[k]]$coefficients)
  }
  expect_equal(one$bridge_fits[[1]]$ratios_validation,two$bridge_fits[[1]]$ratios_validation)
  expect_equal(one$bridge_fits[[1]]$adjoint_fits[[1]]$weights,two$bridge_fits[[1]]$adjoint_fits[[1]]$weights)
})
