require(glmnet)

# Helper function that translates a named weight type
# into a numeric truncation value.
get_trunc_val <- function(weight_type) {
  stopifnot(is.character(weight_type) || 
    is.factor(weight_type))
  trunc <- switch(
    weight_type,
    "trunc@.05" = 0.05,
    "trunc@.01" = 0.01,
    "overlap" = 0.0,
    "notrunc" = 0.0
  )
  if(is.null(trunc))
    stop(paste0(
      "Error in get_trunc_val: truncation level for",
      "weight type ", weight_type, " is undefined."
    ))
  return(trunc)
}


# Helper function that translates a named weight type
# into a character weight class. The weight class
# is used directly by the create_weights function.
get_weight_class <- function(weight_type) {
  stopifnot(is.character(weight_type) || 
    is.factor(weight_type))
  trunc <- switch(
    weight_type,
    "trunc@.05" = "trunc",
    "trunc@.01" = "trunc",
    "overlap" = "overlap",
    "notrunc" = "trunc"
  )
  if(is.null(trunc))
    stop(paste0(
      "Error in get_trunc_val: truncation level for",
      "weight type ", weight_type, " is undefined."
    ))
  return(trunc)
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
  lambda1 = NULL, lambda2 = 0,
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
                                lambda2_vals = 0,
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

      this_oal <- try(
        oal_fitting(X, A, Y, betaXY, gamma,
          lambda1 = lambda1, lambda2 = lambda2,
          exact = exact
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
        # translates the string wgt_type into a truncation level
        # and class--e.g., trunc@.05 => class trunc, trunc=.05,
        # while overlap => class overlap, trunc=0.0 (unused)
        trunc <- get_trunc_val(wgt_type)
        wgt_class <- get_weight_class(wgt_type)
        wgt <- create_weights(this_oal$propensity, A,
          type = wgt_class,
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
            trunc = trunc, wgt_class = wgt_class,
            propensity = this_oal$propensity
          )
        }
      }
    }
  }

  return(out)
}
