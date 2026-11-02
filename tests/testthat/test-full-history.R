fh_test_data <- function(n = 360L) {
  x <- data.frame(C1=rep(0:1,length.out=n),R1=1)
  x$A1 <- rbinom(n,1,.25+.5*x$C1)
  x$C2 <- rbinom(n,1,.1+.25*x$C1+.5*x$A1)
  x$R2 <- rbinom(n,1,.6+.2*x$C2)
  x$A2 <- x$R2*rbinom(n,1,.2+.6*x$C2)
  x$C3 <- rbinom(n,1,.1+.25*x$C2+.5*x$A2)
  x$R3 <- rbinom(n,1,.6+.2*x$C3)
  x
}
fh_test_setup <- function(x, estimator = "sdr", root = seq_len(nrow(x))) {
  alpha <- function(data,trt) rep(0,nrow(data))
  visit <- function(data,trt) rep(1,nrow(data))
  continuation <- function(data,trt,epsilon) as.numeric(epsilon<.35)
  z <- fh_setup(x,c("A1","A2"),c("C2","C3"),c("R1","R2","R3"),
    list("C1","C2","C3"),NULL,1,alpha,visit,continuation,as.list(rep(list(c(.25,.75)),2)),NULL)
  z$randomizer <- function(n,time) rep(c(.2,.4,.8,.1),length.out=n)
  z$ratio_draws<-2L;z$mtp<-FALSE
  lib <- list(sieve=list(method="sieve_md",control=list(target_basis="cell",
    instrument_basis="cell",link="inverse_logit",lambda=.001)))
  control<-lmtp_control(.trim=1,.learners_outcome_folds=2,.learners_trt_folds=2)
  list(z=z, engine=fh_engine(z,root,NULL,"SL.saturated_l1","SL.saturated_l1",
    fh_libraries(lib,NULL),c(-100,100),control,estimator), lib=lib,control=control,
    alpha=alpha,d=visit,t=continuation)
}

test_that("unrecorded health is erased and callbacks cannot see future data", {
  set.seed(123);x<-fh_test_data();a<-fh_test_setup(x)
  x$C2[x$R2==0]<-NA;x$C3[x$R3==0]<-NA;b<-fh_test_setup(x)
  expect_identical(a$z$data,b$z$data)
  expect_true(all(is.na(a$z$data$C2[a$z$data$R2==0])))
  a$z$rules$d[[1]]<-function(data,trt) {
    expect_false(any(c("C2","C3","A2","R2")%in%names(data)))
    data[[trt]]
  }
  expect_equal(as.numeric(fh_policy(a$z,1:20,1,"d")),x$A1[1:20])
  x$A2[x$R2==0]<-1
  expect_error(fh_test_setup(x),"disagrees with alpha")
})

test_that("continuation classifier corrects unequal class sizes", {
  set.seed(124);a<-fh_test_setup(fh_test_data(600))
  rows<-seq_len(a$z$n);nu<-a$engine$nuisance(rows,rows)
  expect_true(all(is.finite(nu$omega)))
  expect_true(all(nu$d[a$z$data$A1==0,1]==0))
  # Its target law puts mass .5 on each treatment. Standard odds without the
  # two-to-one class-size correction would estimate a ratio near two here.
  eligible<-which(a$z$J[,2]);c<-a$z$data$C2[eligible];aa<-a$z$data$A2[eligible]
  truth<-.5/ifelse(aa==1,.2+.6*c,.8-.6*c)
  expect_lt(mean(abs(nu$t[eligible,2]-truth)),.65)
  expect_true(all(nu$omega[a$z$data$R2==0,2]==0))
})

test_that("all nuisance ancestors exclude held-out observations", {
  set.seed(125);x<-fh_test_data(480);root<-1:320;valid<-321:480
  one<-fh_test_setup(x,root=root)
  changed<-x;changed$C2[valid]<-1-changed$C2[valid];changed$C3[valid]<-1-changed$C3[valid]
  two<-fh_test_setup(changed,root=root)
  first<-one$engine$model(root,1L);second<-two$engine$model(root,1L)
  expect_equal(first$coef,second$coef)
  expect_equal(first$library.predict,second$library.predict)
  # An inner validation change must not enter that inner model's training
  # response, even though its validation response changes.
  labels<-one$engine$splits(root);v<-root[labels==1];tr<-root[labels!=1]
  changed<-x;changed$C2[v]<-1-changed$C2[v];changed$C3[v]<-1-changed$C3[v]
  two<-fh_test_setup(changed,root=root)
  for(method in c("sdr","tmle")) {
    one<-fh_test_setup(x,method,root);two<-fh_test_setup(changed,method,root)
    expect_equal(one$engine$response(tr,tr,rep(1,length(tr)),1L),
                 two$engine$response(tr,tr,rep(1,length(tr)),1L),tolerance=1e-10)
  }
})

test_that("full-history SDR and logistic TMLE return complete score components", {
  set.seed(126);x<-fh_test_data();a<-fh_test_setup(x)
  folds<-make_folds(x,V=2)
  for(method in c("sdr","tmle")) {
    fit<-lmtp_full_history(x,c("A1","A2"),c("C2","C3"),c("R1","R2","R3"),
      list("C1","C2","C3"),a$alpha,a$d,a$t,folds=folds,k=1,estimator=method,
      learners_outcome="SL.saturated_l1",learners_trt="SL.saturated_l1",
      bridge_library=a$lib,control=a$control,mtp=FALSE,
      randomizer=a$z$randomizer,ratio_draws=2)
    expect_s3_class(fit,"lmtp_curve");expect_equal(nrow(tidy(fit)),2)
    expect_true(all(is.finite(fit$influence_functions)))
    expect_equal(fit$sequential_estimates+colMeans(fit$bridge_correction),fit$unprojected_estimates)
    expect_equal(dim(fit$bridge_correction),c(nrow(x),2L))
    if(method=="tmle") for(fold in fit$targeting) for(h in fold)
      expect_lt(max(abs(h$score)),1e-7)
  }
})

test_that("fully measured identity policy reduces to empirical means", {
  set.seed(127);x<-fh_test_data(240);x$R2<-x$R3<-1;a<-fh_test_setup(x)
  identity<-function(data,trt) data[,trt,drop=FALSE]
  fit<-lmtp_full_history(x,c("A1","A2"),c("C2","C3"),c("R1","R2","R3"),
    list("C1","C2","C3"),a$alpha,identity,a$t,k=Inf,folds=2,
    learners_outcome="SL.mean",learners_trt="SL.mean",bridge_library=a$lib,
    control=a$control,randomizer=a$z$randomizer,ratio_draws=2)
  expect_equal(fit$unprojected_estimates,unname(colMeans(x[c("C2","C3")])),tolerance=1e-10)
  expect_equal(fit$bridge_correction,matrix(0,nrow(x),2))
  expect_error(lmtp_full_history(x,c("A1","A2"),c("C2","C3"),c("R1","R2","R3"),
    list("C1","C2","C3"),a$alpha,identity,a$t,folds=1),"two outer folds")
})

test_that("SDR and TMLE reuse identical training fits and invalidate changed data", {
  set.seed(128);a<-fh_test_setup(fh_test_data(240));z<-a$z;root<-1:180
  cache<-new.env(parent=emptyenv());libraries<-fh_libraries(a$lib,NULL)
  engine<-function(z,method) fh_engine(z,root,NULL,'SL.mean','SL.mean',libraries,
    c(-100,100),a$control,method,cache)
  one<-engine(z,'sdr');before<-one$nuisance(root,181:240)
  expect_length(ls(cache),1L)
  saved<-one$counts()
  two<-engine(z,'tmle');after<-two$nuisance(root,181:240)
  expect_length(ls(cache),1L)
  expect_identical(before,after)
  expect_equal(saved,two$counts())
  z$data$C3[181:240]<-1-z$data$C3[181:240]
  three<-engine(z,'sdr');expect_length(ls(cache),1L)
  updated<-three$nuisance(root,181:240)
  expect_identical(updated$fits,before$fits)
  expect_gt(max(abs(updated$beta-before$beta)),0)
  z$data$C3[1:180]<-1-z$data$C3[1:180]
  four<-engine(z,'sdr');expect_length(ls(cache),2L)
})

test_that("alpha cannot use current treatment as a circular continuation rule", {
  set.seed(129);a<-fh_test_setup(fh_test_data(80))
  a$z$rules$alpha[[2]]<-function(data,trt) {
    expect_false(trt %in% names(data));rep(0,nrow(data))
  }
  expect_equal(as.numeric(fh_policy(a$z,1:20,2,'alpha')),rep(0,20))
})

test_that("three-time multivariate policies retain cluster splits and finite memory", {
  set.seed(130);n<-160
  data<-data.frame(cluster=rep(seq_len(n/2),each=2))
  for(s in 1:4) {
    data[[paste0('C',s)]]<-rbinom(n,1,.5)
    data[[paste0('X',s)]]<-rnorm(n)
    data[[paste0('R',s)]]<-1
    if(s<4) {
      data[[paste0('A',s)]]<-rbinom(n,1,.5)
      data[[paste0('B',s)]]<-rbinom(n,1,.5)
    }
  }
  trt<-lapply(1:3,function(s) paste0(c('A','B'),s))
  health<-lapply(1:4,function(s) paste0(c('C','X'),s))
  identity<-function(data,trt) data[,trt,drop=FALSE]
  alpha<-function(data,trt) matrix(0,nrow(data),length(trt))
  continuation<-function(data,trt,epsilon) cbind(as.numeric(epsilon<.5),as.numeric(epsilon>=.5))
  control<-lmtp_control(.trim=1,.learners_trt_folds=2,.learners_outcome_folds=2)
  for(k in c(2,Inf)) {
    fit<-lmtp_full_history(data,trt,paste0('C',2:4),paste0('R',1:4),health,
      alpha,identity,continuation,k=k,id='cluster',folds=2,learners_outcome='SL.mean',
      learners_trt='SL.mean',control=control,ratio_draws=1)
    expect_equal(fit$unprojected_estimates,unname(colMeans(data[paste0('C',2:4)])),tolerance=1e-10)
    for(fold in fit$folds) expect_false(any(data$cluster[fold$training_set]%in%data$cluster[fold$validation_set]))
    expect_equal(dim(fit$influence_functions),c(n,3L))
  }
})

test_that("scarce descendant histories retain the failed pilot seed", {
  example<-readRDS(testthat::test_path('fixtures','full-history-scarce-seed.rds'))
  a<-fh_test_setup(example$data)
  control<-lmtp_control(.trim=1,.learners_outcome_folds=3,.learners_trt_folds=3,.discrete=FALSE)
  set.seed(example$seed+2L)
  fit<-lmtp_full_history(example$data,c('A1','A2'),c('C2','C3'),c('R1','R2','R3'),
    list('C1','C2','C3'),a$alpha,a$d,a$t,mtp=FALSE,k=1,folds=example$folds,
    learners_outcome='SL.saturated_l1',learners_trt='SL.saturated_l1',
    bridge_library=a$lib,control=control)
  expect_true(all(is.finite(fit$influence_functions)))
  expect_true(all(is.finite(fit$unprojected_estimates)))
})

test_that("baseline needs no continuation history", {
  set.seed(131);a<-fh_test_setup(fh_test_data(160));z<-a$z
  z$rules$alpha[[1]]<-z$rules$t[[1]]<-function(...) stop('No prescription before baseline.')
  engine<-fh_engine(z,1:120,NULL,'SL.mean','SL.mean',fh_libraries(a$lib,NULL),
    c(-100,100),a$control,'sdr')
  model<-engine$model(1:120,0L)
  expect_true(all(is.finite(engine$predictions(model,121:160,1L)$d)))
})

test_that("alpha is evaluated on the no-visit history even for complete rows", {
  set.seed(132);a<-fh_test_setup(fh_test_data(120));z<-a$z
  z$rules$alpha[[2]]<-function(data,trt) {
    expect_true(all(data$R2==0));expect_true(all(is.na(data$C2)))
    expect_false(trt %in% names(data));rep(0,nrow(data))
  }
  z$rules$t[[2]]<-function(data,trt,epsilon) {
    expect_true(all(data$R2==0));expect_true(all(is.na(data$C2)))
    expect_true(all(data[[trt]]==0));as.numeric(epsilon<.35)
  }
  observed<-which(z$J[,2])
  expect_equal(as.numeric(fh_policy(z,observed,2,'alpha')),rep(0,length(observed)))
  expect_true(all(is.finite(fh_policy(z,observed,2,'t',rep(.2,length(observed))))))
})
