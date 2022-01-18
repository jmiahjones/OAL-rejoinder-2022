
source("./R/oal_funs.R")
source("./R/ate_funs.R")

# test trunc_propensity
props <- c(0.01, 0.1, 0.5, 0.9, 0.99)
stopifnot(
  all(trunc_propensity(props, trunc = 0) == props),
  all(trunc_propensity(props, trunc = 0.05) == c(.05, props[2:4], .95)),
  all(trunc_propensity(props, trunc = 0.2) == c(.2, .2, props[3], .8, .8))
)


# test weight creation
stopifnot(
  all(create_trunc_weights(c(.5, .5), A = c(0, 1), trunc = 0) == c(2, 2)),
  all(create_trunc_weights(c(.2, .8), A = c(0, 1), trunc = 0) == c(1, 1) / .8)
)
create_weights(c(.2, .8), A = c(0, 1), type = "o", trunc = 0)


# test ATE estimation
set.seed(2021)
A <- c(1, 0)
Y <- c(3.1, 1.9)
propensity <- rep(.5, 2)
wgt <- create_trunc_weights(propensity, A, trunc = 0)
ATE_est(Y, A, wgt, est_method = "IPW") == 1.2

# should error -- need to pass Q.hat to ATE with AIPW
error_est <- tryCatch(
  ATE_est(Y, A, wgt, est_method = "AIPW"),
  error = function(e) {
    return(NULL)
  }
)
stopifnot(is.null(error_est))

set.seed(2021)
n <- 10000
propensity <- rep(.5, n)
A <- rbinom(n, 1, propensity)
Y <- A + rnorm(n)
wgt <- create_trunc_weights(propensity, A, trunc = 0)
stopifnot(
  abs(
    ATE_est(Y - A, A, wgt, est_method = "AIPW", Q.hat = matrix(0, nrow = n, ncol = 2))
  ) < 1e-2,
  abs(
    ATE_est(Y, A, wgt, est_method = "AIPW", Q.hat = cbind(rep(0, n), rep(1, n))) - 1
  ) < 1e-2,
  abs(
    ATE_est(Y - A, A, wgt, est_method = "IPW")
  ) < 1e-2,
  abs(
    ATE_est(Y, A, wgt, est_method = "IPW") - 1
  ) < 1e-2,
  abs(
    ATE_est(Y - A, A, wgt,
      est_method = "TMLE",
      Q.hat = matrix(0, nrow = n, ncol = 2),
      propensity = propensity
    )
  ) < 1e-2,
  abs(
    ATE_est(Y, A, wgt,
      est_method = "TMLE",
      Q.hat = cbind(rep(0, n), rep(1, n)),
      propensity = propensity
    ) - 1
  ) < 1e-2
)
