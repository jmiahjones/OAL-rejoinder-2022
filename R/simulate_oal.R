# modified by Jeremiah Jones, 2021 to investigate positivity violations
options(future.rng.onMisuse = "ignore")

library(MASS)
# library(glmnet)
library(lqa)
library(dplyr)
library(foreach)
library(doFuture)
library(purrr)
library(tidyr)
library(mlr3)
library(mlr3learners)
library(data.table)
lgr::get_logger("mlr3")$set_threshold("warn")
plan(multicore, workers = parallel::detectCores() - 1)
# plan(sequential, split = TRUE)
registerDoFuture()


source("./R/oal_funs.R")
source("./R/ate_funs.R")

# n=200; p=100; num_simulations=100; rho=0; sig_x=.6; scenario=1
# sel_method="OAL"; est_method="IPW"; verbose=1; pos_viol_cut=.05
simulate_oal <- function(n, p, num_simulations = 100L,
                         rho, sig_x, scenario,
                         sel_method="OAL",
                         est_method = c("IPW", "AIPW", "TMLE"),
                         weight_types = c("trunc@.05", "trunc@.01", "overlap", "notrunc"),
                         verbose = 1,
                         pos_viol_cut = .05) {
  use_ridge = (sel_method == "GOALn")
  lambda2_vals <- if (use_ridge) {
    c(0, 
      # 10^c(-2,-1.5,-1,-.75,-.5,-.25,0,.25,.5,1)
      10^-2:1
    )
  } else {
    0
  }
  message(paste("Arguments: ", n, p, num_simulations,
                         rho, sig_x, scenario,
                         sel_method, sep = " "))
  stopifnot(all(est_method %in% c("IPW", "AIPW", "TMLE")))
  

  # set information for simulating coviariates
  mean_x <- 0
  # sig_x = 1/4
  # sig_x = 1
  # rho = 0.0
  # sample size
  # n = 200
  # total number of predictors
  # p = 100
  pC <- pI <- pP <- 2
  pS <- p - (pC + pI + pP)
  var.list <- c(paste("Xc", 1:pC, sep = ""), paste("Xp", 1:pP, sep = ""), paste("Xi", 1:pI, sep = ""), paste("Xs", 1:pS, sep = ""))
  var.list.plus.int <- c("(Intercept)", var.list)

  # Set strength of relationship between covariates and outcome (beta_v)
  # Set strength of relationship between covariates and treatment (alpha_v)
  switch(scenario,
    "1" = {
      beta_v <- c(0.6, 0.6, 0.6, 0.6, 0, 0, rep(0, p - 6)) #* 2
      alpha_v <- c(1.0, 1.0, 0, 0, 1, 1, rep(0, p - 6))
    },
    "2" = {
      beta_v <- c(0.6, 0.6, 0.6, 0.6, 0, 0, rep(0, p - 6)) #* 2
      alpha_v <- c(0.4, 0.4, 0, 0, 1, 1, rep(0, p - 6))
    }
  )

  names(beta_v) <- names(alpha_v) <- var.list
  ### set true average treatment effect
  bA <- 0

  # set vector of possible lambda's to try
  lambda_vec <- c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) <- as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor <- 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals <- 2 * (gamma_convergence_factor - lambda_vec + 1)
  names(gamma_vals) <- names(lambda_vec)
  # gamma_vals <- c(0.5, 1, 2, 3)

  ### define some functions for generating data, ATE estimates, and the wAMD,
  expit <- plogis
  fam <- binomial()

  result <- foreach(
    sim_idx = seq.int(num_simulations),
    .combine = dplyr::bind_rows,
    .errorhandling = "remove"
  ) %dopar% {
    
    set.seed(2021 + sim_idx)
    ### simulate data
    Sigma_x <- matrix(rho * sig_x^2, nrow = p, ncol = p)
    diag(Sigma_x) <- sig_x^2
    # Sigma_x[1,2] <- Sigma_x[2,1] <- 0.99*sig_x^2 # make the two confounders nearly-collinear
    Mean_x <- rep(mean_x, p)
    Data <- as.data.frame(mvrnorm(n = n, mu = Mean_x, Sigma = Sigma_x, empirical = FALSE))
    # Data <- scale(Data, T, T)
    names(Data) <- var.list

    gA_x <- as.numeric(as.matrix(Data[, var.list]) %*% alpha_v)
    pA <- expit(gA_x)
    summary(pA)
    prop_positivity <- mean(pA < pos_viol_cut)
    Data$A <- as.numeric(runif(n = length(pA)) < pA) # simulate A
    gY_xA <- as.numeric(as.matrix(Data[, var.list]) %*% beta_v)
    Data$Y <- gY_xA + rnorm(n = n)
    Data$Y <- Data$Y + Data$A * bA

    X <- as.matrix(Data[, 1:p])
    X <- scale(X, T, T)

    # estimate outcome model
    y.form <- formula(paste("Y~A+", paste(var.list, collapse = "+")))
    lm.Y <- lm(y.form, data = Data)
    betaXY <- coef(lm.Y)[var.list]
    unpen <- betaXY

    ######################################################################################
    #####  Run outcome adaptive lasso for each lambda value
    ######################################################################################
    grid_min <- grid_search_oal_fit(
      gamma_vals, 
      lambda_vec = n^(lambda_vec),
      X = X, A = Data$A, Y = Data$Y,
      betaXY = betaXY,
      lambda2_vals = lambda2_vals,
      weight_types = weight_types
    )
    # save coefficients, wamd metric, and hyperparameters
    coefs <- lapply(grid_min, function(x) x$coefs) %>% do.call(c, .)
    sel_idx <- rep(seq.int(1+p), length(grid_min))
    stopifnot(length(coefs) == length(sel_idx))
    
    wamds <- purrr::map_dbl(grid_min, ~.$wamd)
    lambda2_chosens <- purrr::map_dbl(grid_min, ~.$lambda2)
    lambda1_chosens <- purrr::map_dbl(grid_min, ~.$lambda1)
    gamma_chosens <- purrr::map_dbl(grid_min, ~.$gamma)
    truncs <- purrr::map_dbl(grid_min, ~.$trunc)
    wgt_classes <- purrr::map_chr(grid_min, ~.$wgt_class)
    
    # validation: ensure that the weight types are in 
    # the exact order in which they were requested
    stopifnot(all(purrr::map_chr(grid_min, ~.$wgt_type)
      == weight_types))
    positivity_tol <- min(min(pA), min(1 - pA))
    
    Q <- cross_fitting(Data, target = "Y", folds = 5L,
      learner = "regr.lm")
        # estimate the ATE using different methods, for each of the created weights
    ate_result <- tibble(
      weight_type = weight_types,
      wgt = purrr::map(grid_min, ~.$wgt),
      propensity = purrr::map(grid_min, ~.$propensity)
    )
    ate_result <- tidyr::expand_grid(est_method = est_method, ate_result)
    ate_result <- ate_result %>%
      mutate(
        ates = purrr::pmap_dbl(., ATE_est, Y = Data$Y, A = Data$A, 
          Q.hat = Q
        )
      )
    sel_result <- tibble(
      coefs = coefs,
      sel_idx = sel_idx
    ) %>%
    tidyr::nest(sel_result = c(coefs, sel_idx))
    
    # compress ate results down to length(weight_types) to
    # match the length of grid_min
    ate_result <- tidyr::nest(ate_result, 
      data = c(est_method, wgt:ates)
    )
    sim_result <- tibble(
      simulation = sim_idx,
      method_results = list(tibble(
        ate_result,
        lambda2_chosens = lambda2_chosens,
        lambda1_chosens = lambda1_chosens,
        gamma_chosens = gamma_chosens,
        trunc = truncs,
        wgt_class = wgt_classes,
        wamds = wamds,
        sel_result
      )),
      outcome_coefs = list(tibble(unpen = unpen)),
      prop_positivity = prop_positivity,
      positivity_tol = min(min(pA), min(1 - pA))
    )

    if (
      (verbose == 3) | 
      (verbose == 2 & sim_idx %% 10 == 1) |
      (verbose == 1 & sim_idx %% 100 == 1)
    ) {
      print(sprintf("Completed simulation %i.", sim_idx))
    }

    return(sim_result)
  }

  return(result)
}

# demo <- simulate_oal(200, 20, 4, rho = 0.75, sig_x = 1, scenario = 1, sel_method = "OAL", verbose=2)
