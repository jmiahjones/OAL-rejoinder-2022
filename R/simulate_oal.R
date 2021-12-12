# modified by Jeremiah Jones, 2021 to investigate positivity violations
options(future.rng.onMisuse = "ignore")

library(MASS)
# library(glmnet)
library(lqa)
library(dplyr)
library(foreach)
library(doFuture)
library(purrr)
library(mlr3)
library(mlr3learners)
library(mlr3pipelines)
library(data.table)
lgr::get_logger("mlr3")$set_threshold("warn")
plan(multicore, workers = 4)
registerDoFuture()


source("./R/oal_funs.R")

n=200; p=100; num_simulations=100; rho=0; sig_x=.6; scenario=1
sel_method="OAL"; est_method="IPW"; verbose=1; pos_viol_cut=.05
simulate_oal <- function(n, p, num_simulations = 100L,
                         rho, sig_x, scenario,
                         sel_method="OAL",
                         est_method="IPW",
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

  stopifnot(all(est_method %in% c("IPW", "AIPW", "Overlap", "Trunc", "TMLE")))
  num_ests <- length(est_method)
  

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

  sim_res <- foreach(sim_idx = seq.int(num_simulations)) %dopar% {
    
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
    # TODO: create weight_types
    weight_types = c("overlap", "trunc")
    grid_min <- grid_search_oal_fit(
      gamma_vals, 
      lambda_vec = n^(lambda_vec),
      X = X, A = Data$A, Y = Data$Y,
      betaXY = betaXY,
      lambda2_vals = lambda2_vals,
      weight_types = weight_types
    )
    grid_min_exact <- grid_min

    wamds <- purrr::map_dbl(grid_min_exact, ~.$wamd)
    lambda2_chosens <- purrr::map_dbl(grid_min_exact, ~.$lambda2)
    lambda1_chosens <- purrr::map_dbl(grid_min_exact, ~.$lambda1)
    positivity_tol <- min(min(pA), min(1 - pA))

    # save coefficients
    coefs <- grid_min_exact[[1]]$coefs
    propensity <- plogis(coefs[1] + X %*% coefs[-1])
    

    Y_task <- as_task_regr(Data, target = "Y")
    lm_lrn <- lrn("regr.lm")
    resampling <- rsmp("cv", folds=5L)
    rr = resample(Y_task, lm_lrn, resampling, store_models = TRUE)
    
    # Q <- matrix(data = NA, nrow = n, ncol = 2)
    # colnames(Q) <- c("Q0", "Q1")
    # for(i in 1:5) {
    #   test_idxs <- rr$resampling$test_set(i)
    #   newdata <- Data[test_idxs,]; newdata$A <- 0
    #   Q[test_idxs, 1] <- rr$learners[[i]]$
    #     predict_newdata(newdata)
    #   newdata$A <- 1
    #   Q[test_idxs, 2] <- rr$learners[[i]]$
    #     predict_newdata(newdata)
    # }
    
    pred_df <- as.data.table(rr$prediction())
    setorder(pred_df, "row_ids")
    Q.hat <- pred_df$response

    # estimate the ATE using different methods, for each of the created weights
    # TODO: Make ATE_est take an est_method and weight_type argument first
    sim_result <- tibble(
      est_method = est_method, weight_type = weight_types,
      wgt = purrr::map(grid_min_exact, ~.$wgt),
      propensity = purrr::map(grid_min_exact, ~.$propensity)
    )
    sim_result <- sim_result %>%
      mutate(
        ates = purrr::pmap(., ATE_est, Y = Data$Y, A = Data$A, 
          Q.hat = Q.hat
        )
      )

    sim_result <- tibble(
      sim_result, 
      lambda2_chosens = lambda2_chosens,
      lambda1_chosens = lambda1_chosens,
      coefs = coefs,
      wamds = wamds,
      unpen = unpen,
      prop_positivity = prop_positivity,
      positivity_tol = min(min(pA), min(1 - pA))
      
      # TODO: Put in other necessary columns for packing,
      # then modify the unpacking to use tibbles rather than lists
    )
    

    if (verbose == 2 | (verbose == 1 & sim_idx %% 10 == 1)) {
      print(sprintf("Completed simulation %i.", sim_idx))
    }

    # pack the results to send back
    list(
      ates = ates,
      lambda2_chosens = lambda2_chosens,
      lambda1_chosens = lambda1_chosens,
      coefs = coefs,
      wamds = wamds,
      unpen = unpen,
      prop_positivity = prop_positivity,
      positivity_tol = min(min(pA), min(1 - pA))
    )
  }
  
  ates <- vector("numeric", num_simulations)
  lambda1_chosens <- vector("numeric", num_simulations)
  lambda2_chosens <- vector("numeric", num_simulations)
  coefs <- matrix(NA, nrow = num_simulations, ncol = 1 + p)
  wamds <- vector("numeric", num_simulations)
  unpen <- matrix(NA, nrow = num_simulations, ncol = p)
  prop_positivity <- rep(NA, num_simulations)
  colnames(coefs) <- var.list.plus.int

  # unpack the reduced results
  positivity_tol <- sapply(sim_res, function(x) x$positivity_tol)
  positivity_tol <- min(positivity_tol)
  for(sim_idx in seq.int(num_simulations)){
    ates[sim_idx] = sim_res[[sim_idx]]$ates
    coefs[sim_idx, ] = sim_res[[sim_idx]]$coefs
    wamds[sim_idx] = sim_res[[sim_idx]]$wamds
    lambda2_chosens[sim_idx] = sim_res[[sim_idx]]$lambda2_chosens
    lambda1_chosens[sim_idx] = sim_res[[sim_idx]]$lambda1_chosens
    unpen[sim_idx, ] = sim_res[[sim_idx]]$unpen
    prop_positivity[sim_idx] = sim_res[[sim_idx]]$prop_positivity 
  }
  # clean the baggage
  rm(sim_res)


  is_selected <- 1 * (abs(coefs) > 1e-6)

  result <- tibble(
    positivity_tol = positivity_tol,
    prop_positivity = mean(prop_positivity),
    method = ifelse(use_ridge, "GOAL", "OAL"),
    metrics = list(tibble(
      wamd = wamds, ate = ates,
      lambda1 = lambda1_chosens,
      lambda2 = lambda2_chosens
    )),
    sel_perf = list(tibble(
      sel_idx = 1:p,
      prop_sel = colMeans(is_selected)[-1],
      coefs = colMeans(coefs[, -1])
    ))
  )

  return(result)
}

# foo<-simulate_oal(200, 100, 100, 0.5, 1, 1, T, 2, F)
# bar<-simulate_oal(200, 100, 100, 0.5, 1, 1, F, 2, T)
# foo$metrics[[1]]
# bar$metrics[[1]]

# mean(foo$metrics[[1]]$ate)
# mean(bar$metrics[[1]]$ate)
