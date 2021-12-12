require(glmnet)
require(mlr3)

ATE_est <- function(Y, A, wgt, est_method, weight_type, ...) {
  if(est_method == "TMLE") {
    # tmle code
    ATE_est_tmle(Y, A, wgt, propensity)
  } else if(est_method == "IPW") {
    ATE_est_IPW(Y, wgt, A)
  } else if(est_method == "AIPW") {
    ATE_est_AIPW(Y, wgt, A)
  }
}

ATE_est_tml <- function(Y, A, wgt, propensity){
  H1 <- A
  H0 <- 1-A
  # scale Y 0 to 1
  Ymin <- min(Y)
  Y <- Y - Ymin
  Ymax <- max(Y)
  Y <- Y / Ymax
  
  suppressWarnings(
    tmle_coef <- coef(glm(Y ~ -1 + offset(Q.hat) + H1 + H0, weights = wgt))
  )
  Q1star <- plogis(Q.hat + H1*tmle_coef[1])
  Q0star <- plogis(Q.hat + H0*tmle_coef[2])

  # unscale
  Y = Ymax * Y + Ymin
  Q1 = Ymax * Q1star + Ymin
  Q0 = Ymax * Q0star + Ymin

  mean(Q1 - Q0)
}

ATE_est_IPW <- function(Y, wgt, A) {
  num_ATE <- Y * wgt
  res <- ((sum(num_ATE[A == 1]) / sum(wgt[A == 1])) - (sum(num_ATE[A == 0]) / sum(wgt[A == 0])))
  return(res)
}

ATE_est_AIPW <- function(Y, X, A, wgt, Q.hat) {
  eps.hat <- Y - Q.hat
  influences <- wgt*eps.hat + Q.hat
  mean(influences)
}

trunc_propensity <- function(propensity, trunc = 0.0) {
  stopifnot(trunc >= 0.0)
  if (trunc > 0.0) {
    propensity <- pmax(propensity, trunc)
    propensity <- pmin(propensity, 1 - trunc)
  }
  return(propensity)
}

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

create_trunc_weights <- function(propensity, A, trunc = 0.0, ...) {
  propensity <- trunc_propensity(propensity, trunc)
  weights <- rep(0, length(A))
  weights[A == 1] <- 1 / propensity[A == 1]
  weights[A == 0] <- 1 / (1 - propensity[A == 0])
  return(weights)
}

create_overlap_weights <- function(propensity, A, ...) {
  weights <- abs( A - propensity )
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
  #       Additionally, to make it easier to compare glmnet and lqa,
  #       let's scale by 1/n internally in this function rather than
  #       externally in the simulate_oal function.
  d <- dim(X)
  lam2_adj <- ifelse(lambda2 == 0, 1, d[1] / sum(d))
  lam1_adj <- 1 / d[1]
  s <- mean(pen) / lam2_adj * lam1_adj * lambda1

  # concatenate the call_list with the exact argument
  coef_call_list <- c(list(logit_oal), call_list, list(exact = exact, s = s))
  coefs_all <- do.call(coef, coef_call_list)
  coefs_all <- (1 + lambda2) * coefs_all
  
  predict_call_list <- c(
    coef_call_list,
    list(type = "response", newx = X)
  )
  est_propens <- do.call(predict, predict_call_list)
  min_coefs <- as.numeric(coefs_all)

  ret <- list(propensity = est_propens, coefs = min_coefs)
  return(ret)
}

grid_search_oal_fit <- function(gamma_vals, lambda_vec,
                                X, A, Y, betaXY,
                                trunc_vals = 0, lambda2_vals = 0,
                                exact = F, weight_types = "trunc") {
  stopifnot(length(lambda_vec) == length(gamma_vals))
  p <- ncol(X)

  num_weights <- length(weight_types)
  out <- lapply(1:num_weights, function(x) list(wamd = Inf))

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
            is.null(this_oal)
        ) {
          message(sprintf(
            "Failed: lambda2=%.2f, lambda1=%.2f, gamma=%.2f",
            lambda2, lambda1, gamma
          ))
          next
        }

        for (wgt_idx in seq.int(num_weights)) {
          wgt_type <- weight_types[wgt_idx]
          wgt <- create_weights(this_oal$propensity, A,
            type = wgt_type,
            trunc = trunc
          )

          # estimate weighted absolute mean difference over all covariates using this lambda to generate weights
          wAMD <- wAMD_function(
            X = X, A = A,
            wgt = wgt, beta = betaXY
          )$wAMD

          if (wAMD < out[[wgt_idx]]$wamd) {
            out[[wgt_idx]] <- list(
              wgt_type = wgt_type,
              wamd = wAMD, coefs = this_oal$coefs,
              wgt = wgt, lambda1 = lambda1,
              gamma = gamma, lambda2 = lambda2,
              trunc = trunc
            )
          }
        }
      }
    }
  }
  
  return(out)
}
