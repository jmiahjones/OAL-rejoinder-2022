require(glmnet)

ATE_est <- function(Y, wgt, A) {
  num_ATE <- Y * wgt
  res <- ((sum(num_ATE[A == 1]) / sum(wgt[A == 1])) - (sum(num_ATE[A == 0]) / sum(wgt[A == 0])))
  return(res)
}

trunc_propensity <- function(propensity, trunc = 0.0) {
  stopifnot(trunc >= 0.0)
  if (trunc > 0.0) {
    propensity <- pmax(propensity, trunc)
    propensity <- pmin(propensity, 1 - trunc)
  }
  return(propensity)
}
create_trunc_weights <- function(propensity, A, trunc = 0.0) {
  propensity <- trunc_propensity(propensity, trunc)
  weights <- rep(0, length(A))
  weights[A == 1] <- 1 / propensity[A == 1]
  weights[A == 0] <- 1 / (1 - propensity[A == 0])
  return(weights)
}
wAMD_function <- function(X, A, wgt, beta) {
  diff_vec <- sapply(seq.int(ncol(X)), function(jj) {
    this_var <- X[, jj] * wgt
    abs((sum(this_var[A == 1]) / sum(wgt * A)) - (sum(this_var[A == 0]) / sum(wgt * (1 - A))))
  })
  wdiff_vec <- diff_vec * abs(beta)
  wAMD <- sum(wdiff_vec)
  ret <- list(diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD)
  return(ret)
}

oal_fitting <- function(
  X, A, Y, betaXY, gamma, 
  lambda1 = NULL, lambda2 = 0, trunc = 0.0,
  exact = F
) {
  stopifnot(length(lambda2) == 1, lambda2 >= 0)
  if (!is.null(lambda1)) {
    stopifnot(length(lambda1) == 1)
  }
  if (lambda2 == 0) {
    X1 <- X
    A1 <- A
  } else {
    X1 <- rbind(X, sqrt(lambda2) * diag(ncol(X)))
    A1 <- c(A, rep(0, ncol(X)))
  }
  pen <- abs(betaXY)^(-gamma)
  
  call_list <- list(
    x = X1, y = A1, family = binomial(), intercept = T,
    standardize = F, penalty.factor = pen
    # , type.logistic = "modified.Newton"
  )
  logit_oal <- do.call(glmnet::glmnet, call_list)

  # NOTE: lambda1 is on a different scale than the likelihood function
  #       ln in glmnet. Because we added artificial samples, we need
  #       to correct for this by making sure we multiply by (n+p)/n.
  d <- dim(X)
  lam2_adj <- ifelse(lambda2 == 0, 1, d[1] / sum(d))
  s <- mean(pen) / lam2_adj * lambda1

  # concatenate the call_list with the exact argument
  coef_call_list <- c(list(logit_oal), call_list, list(exact = exact, s = s))
  coefs_all <- do.call(coef, coef_call_list)
  coefs_all <- (1 + lambda2) * coefs_all
  
  predict_call_list <- c(
    coef_call_list,
    list(type = "response", newx = X)
  )
  est_propens <- do.call(predict, predict_call_list)

  wgt <- create_trunc_weights(est_propens, A = A, trunc = trunc)

  # estimate weighted absolute mean different over all covariates using this lambda to generate weights
  wAMD <- wAMD_function(
    X = X, A = A,
    wgt = wgt, beta = betaXY
  )$wAMD

  min_coefs <- as.numeric(coefs_all)
  min_ate <- ATE_est(Y = Y, wgt = wgt, A = A)

  ret <- list(ate = min_ate, wAMD = wAMD, coefs = min_coefs)
  return(ret)
}

oal_fitting <- function(
  X, A, Y, betaXY, gamma, 
  lambda1 = NULL, lambda2 = 0, trunc = 0.0,
  exact = F
) {
  require(lqa)
  stopifnot(length(lambda2) == 1, lambda2 >= 0)
  if (!is.null(lambda1)) {
    stopifnot(length(lambda1) == 1)
  }
  if (lambda2 == 0) {
    X1 <- X
    A1 <- A
  } else {
    X1 <- rbind(X, sqrt(lambda2) * diag(ncol(X)))
    A1 <- c(A, rep(0, ncol(X)))
  }
  pen <- abs(betaXY)^(-gamma)
  
  oal_pen = adaptive.lasso(lambda = lambda1, al.weights = pen)
  # run outcome-adaptive lasso model with appropriate penalty
  logit_oal = lqa.default(
    x = X1, y = A1, penalty = oal_pen, family = binomial(logit),
    intercept = T, standardize = F
  )
  
  coefs_all <- coef(logit_oal)
  coefs_all <- (1 + lambda2) * coefs_all
  
  est_propens <- plogis(cbind(1, X) %*% coefs_all)

  wgt <- create_trunc_weights(est_propens, A = A, trunc = trunc)

  # estimate weighted absolute mean different over all covariates using this lambda to generate weights
  wAMD <- wAMD_function(
    X = X, A = A,
    wgt = wgt, beta = betaXY
  )$wAMD

  coefs_all <- as.numeric(coefs_all)
  min_ate <- ATE_est(Y = Y, wgt = wgt, A = A)

  ret <- list(ate = min_ate, wAMD = wAMD, coefs = coefs_all)
  return(ret)
}

grid_search_oal_fit <- function(
  gamma_vals, lambda_vec,
  X, A, Y, betaXY,
  trunc_vals = 0, lambda2_vals = 0,
  exact = F
) {
  stopifnot(length(lambda_vec) == length(gamma_vals))
  p <- ncol(X)
  min_wamd <- Inf
  tmp_coefs <- rep(NA, 1 + p)
  tmp_ate <- NA
  tmp_lam2 <- NA
  tmp_trunc <- NA
  tmp_lam1 <- NA
  tmp_gam <- NA
    

  for (lam2_idx in seq_along(lambda2_vals)) {
    lambda2 <- lambda2_vals[lam2_idx]
    for (gam_idx in seq_along(gamma_vals)) {
      gamma <- gamma_vals[gam_idx]
      # gamma <- 1
      lambda1 <- lambda_vec[gam_idx]

      for (trunc_idx in seq_along(trunc_vals)) {
        trunc <- trunc_vals[trunc_idx]

        this_oal <- try(
          oal_fitting(X, A, Y, betaXY, gamma,
            lambda1 = lambda1, lambda2 = lambda2,
            trunc = trunc, exact = exact
          ),
          silent = T
        )

        if (
          inherits(this_oal, "try-error") || 
          is.null(this_oal) ||
          is.na(this_oal$ate)
        ) {
          message(sprintf(
            "Failed: lambda2=%.2f, lambda1=%.2f, gamma=%.2f",
            lambda2, lambda1, gamma
          ))
          next
        }

        if (this_oal$wAMD < min_wamd) {
          min_wamd <- this_oal$wAMD
          tmp_coefs <- this_oal$coefs
          tmp_ate <- this_oal$ate
          tmp_lam2 <- lambda2
          tmp_trunc <- trunc
          tmp_lam1 <- lambda1
          tmp_gam <- gamma
        }
      }
    }
  }

  return(list(
    wamd = min_wamd, coefs = tmp_coefs,
    ate = tmp_ate, lambda1 = tmp_lam1,
    gamma = tmp_gam, lambda2 = tmp_lam2,
    trunc = tmp_trunc
  ))
}