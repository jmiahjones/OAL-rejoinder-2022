# ATE estimation function that should be used directly.
# Minimally, takes Y, A, weights, and
# a method string that then calls the relevant
# estimator.
ATE_est <- function(Y, A, wgt, est_method, Q.hat=NULL, propensity=NULL, ...) {
  if(est_method %in% c("TMLE", "AIPW")){
    stopifnot(
      !is.null(Q.hat), ncol(Q.hat)==2
    )
    if(est_method == "TMLE")
      stopifnot(!is.null(propensity))
  }
  ate <- if(est_method == "TMLE") {
    ATE_est_tmle(Y, A, wgt, Q.hat, propensity)
  } else if(est_method == "IPW") {
    ATE_est_IPW(Y, A, wgt)
  } else if(est_method == "AIPW") {
    ATE_est_AIPW(Y, A, wgt, Q.hat)
  }
  return(ate)
}

scale_01 <- function(Y, lims) {
  (Y - lims[1]) / diff(lims)
}

inv_scale_01 <- function(Y, lims) {
  Y * diff(lims) + lims[1]
}

# Estimate the ATE with TMLE method.
# Currently uses a linear fluctuation
# since we have a parametric linear model
# providing the Q estimates.
ATE_est_tmle <- function(Y, A, wgt, Q.hat, propensity, trunc=0.01){
  QAX <- A * Q.hat[,2] + (1-A) * Q.hat[,1]
  H1 <- A
  H0 <- 1-A

  suppressWarnings(
    tmle_fit <- (lm(Y ~ -1 + offset(QAX) + H0 + H1, weights = wgt))
  )
  tmle_coef <- coef(tmle_fit)
  Q1star <- Q.hat[,2] + tmle_coef[2]
  Q0star <- Q.hat[,1] + tmle_coef[1]
  
  mean(Q1star - Q0star)
}

# Estimate the ATE with a weighted sum.
# Uses the weights in a convex combination.
ATE_est_IPW <- function(Y, A, wgt) {
  wgt1 <- A * wgt
  wgt0 <- (1-A) * wgt
  res <- sum(wgt1 * Y) / sum(wgt1) - (sum(wgt0 * Y) / sum(wgt0))
  return(res)
}

# Estimate the ATE with AIPW estimator.
# Plug-in the weights and the outcome regression
# into the efficient influence function.
# Note that we are not necessarily using the
# efficient influence function if AIPW is used
# without the IPW weights.
ATE_est_AIPW <- function(Y, A, wgt, Q.hat) {
  QAX <- A * Q.hat[,2] + (1-A) * Q.hat[,1]
  eps.hat <- Y - QAX
  influences <- (-(1-A) + A)*wgt*eps.hat + Q.hat[,2] - Q.hat[,1]
  mean(influences)
}

# Helper function: takes a propensity and truncation
# level and returns the truncated propensities.
# When called without a trunc argument, this simply
# returns the propensity.
trunc_propensity <- function(propensity, trunc = 0.0) {
  stopifnot(trunc >= 0.0)
  if (trunc > 0.0) {
    propensity <- pmax(propensity, trunc)
    propensity <- pmin(propensity, 1 - trunc)
  }
  return(propensity)
}

# Weight creation function that should be used directly.
# Dispactches to the truncated weight and overlap weight
# creation functions.
create_weights <- function(propensity, A, type = c("trunc", "overlap"), ...) {
  if(length(type) > 1) {
    type = type[1]
  }

  weights <- switch(
    tolower(substr(type, 1, 1)),
    t = create_trunc_weights(
      propensity=propensity, A=A, ...
    ),
    o = create_overlap_weights(
      propensity=propensity, A=A, ...
    )
  )
  return(weights)
}

# Helper function: creates IPW at the desired truncation level.
create_trunc_weights <- function(propensity, A, trunc = 0.0, ...) {
  propensity <- trunc_propensity(propensity, trunc)
  weights <- rep(0, length(A))
  weights[A == 1] <- 1 / propensity[A == 1]
  weights[A == 0] <- 1 / (1 - propensity[A == 0])
  return(weights)
}

# Helper function: creates the overlap weights.
create_overlap_weights <- function(propensity, A, ...) {
  weights <- abs( A - propensity )
  return(weights)
}

cross_fitting <- function(
  Data, target, folds, learner
) {
  require(mlr3)
  require(mlr3learners)

  stopifnot(
    is.character(target) && length(target) == 1,
    is.character(learner) && length(learner) == 1,
    is.numeric(folds) && length(folds) == 1
  )

  # create Q outcome regression matrix for TMLE
  Y_task <- mlr3::as_task_regr(Data, target = target)
  # using parametric model for simplicity
  lm_lrn <- mlr3::lrn(learner)
  resampling <- mlr3::rsmp("cv", folds = folds)
  rr <- mlr3::resample(Y_task, lm_lrn, resampling, store_models = TRUE)

  Q <- matrix(data = NA, nrow = n, ncol = 2)
  colnames(Q) <- c("Q0", "Q1")
  for (i in seq.int(folds)) {
    test_idxs <- rr$resampling$test_set(i)
    newdata <- Data[test_idxs, ]
    newdata$A <- 0
    Q[test_idxs, 1] <- rr$learners[[i]]$
      predict_newdata(newdata, task = Y_task)$response
    newdata$A <- 1
    Q[test_idxs, 2] <- rr$learners[[i]]$
      predict_newdata(newdata, task = Y_task)$response
  }
  return(Q)
}