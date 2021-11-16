# modified by Jeremiah Jones, 2021 to investigate positivity violations
library(MASS)
library(glmnet)
library(dplyr)

source("./R/oal_funs.R")

n <- 200
p <- 20
num_simulations <- 5L
rho <- 0.00
sig_x <- 1
scenario <- 1L
use_ridge <- F
verbose <- 1L

lambda2_vals <- if (use_ridge) {
  c(0, 10^(-2:2))
} else {
  0
}


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
    beta_v <- c(0.6, 0.6, 0.6, 0.6, 0, 0, rep(0, p - 6)) * 2
    alpha_v <- c(1.0, 1.0, 0, 0, 1, 1, rep(0, p - 6))
  },
  "2" = {
    beta_v <- c(0.6, 0.6, 0.6, 0.6, 0, 0, rep(0, p - 6)) * 2
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

# truncation grid
trunc_vals <- seq(0, 1, by=.2)

### define some functions for generating data, ATE estimates, and the wAMD,
expit <- plogis
fam <- binomial()


# ates <- sapply(1L:num_simulations, function(sim_idx){
ates <- vector("numeric", num_simulations)
lambda2_chosens <- vector("numeric", num_simulations)
coefs <- matrix(NA, nrow = num_simulations, ncol = 1 + p)
wamds <- vector("numeric", num_simulations)
unpen <- matrix(NA, nrow = num_simulations, ncol = p)
colnames(coefs) <- var.list.plus.int
trunc_chosens <- vector("numeric", num_simulations) #TODO
for (sim_idx in seq.int(num_simulations)) {
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
  unpen[sim_idx, ] <- betaXY

  coeff_XA <- (matrix(NA, nrow = 1 + p, ncol = 100))
  rownames(coeff_XA) <- var.list.plus.int

  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value
  ######################################################################################

  grid_min <- grid_search_oal_fit(
    gamma_vals, lambda_vec = n^(lambda_vec - 1), 
    X = X, A = Data$A, Y = Data$Y,
    betaXY = betaXY
  )
  
  grid_min_exact <- grid_search_oal_fit(
    gamma_vals = grid_min$gamma, lambda_vec = grid_min$lambda1, 
    X = X, A = Data$A, Y = Data$Y,
    betaXY = betaXY, exact = T,
    lambda2_vals = grid_min$lambda2,
    trunc_vals = grid_min$trunc
  )

  
  wamds[sim_idx] <- grid_min_exact$wamd
  # save coefficients
  coefs[sim_idx, ] <- grid_min_exact$coefs

  # print out ATE corresponding to smallest wAMD value
  ates[sim_idx] <- grid_min_exact$ate

  lambda2_chosens[sim_idx] <- grid_min_exact$lambda2

  if (verbose == 2 | (verbose == 1 & sim_idx %% 10 == 1)) {
    message(sprintf("Completed simulation %i.", sim_idx))
  }
}

# mean(ates)
# colMeans(coefs[,1:8])

is_selected <- 1 * (abs(coefs) > 1e-6)

result <- tibble(
  positivity_tol = min(min(pA), min(1 - pA)),
  method = ifelse(use_ridge, "GOAL", "OAL"),
  metrics = list(tibble(wamd = wamds, ate = ates, lambda2 = lambda2_chosens)),
  sel_perf = list(tibble(
    sel_idx = 1:p,
    prop_sel = colMeans(is_selected)[-1],
    coefs = colMeans(coefs[, -1])
  ))
)