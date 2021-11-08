# modified by Jeremiah Jones, 2021 to investigate positivity violations
library(MASS)
library(glmnet)
library(dplyr)


simulate_oal <- function(
  n, p, num_simulations=100L, 
  rho, sig_x, scenario, use_ridge=F,
  verbose=2
) {
  
  lambda2_vals <- if(use_ridge){
    c(0, 10^(-2:1))
  } else {
    0
  }
  
  
  # set information for simulating coviariates
  mean_x = 0 
  # sig_x = 1/4
  # sig_x = 1
  # rho = 0.0
  # sample size
  # n = 200
  # total number of predictors
  # p = 100
  pC = pI = pP = 2
  pS = p - (pC+pI+pP)
  var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
  var.list.plus.int <- c("(Intercept)", var.list)
  
  # Set strength of relationship between covariates and outcome (beta_v)
  # Set strength of relationship between covariates and treatment (alpha_v)
  switch (scenario,
          "1" = {
            beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )*2
            alpha_v = c( 1.0, 1.0,   0,   0, 1, 1,  rep(0,p-6) )
          },
          "2" = {
            beta_v =  c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )*2
            alpha_v = c( 0.4, 0.4,   0,   0, 1, 1,  rep(0,p-6) )
          }
  )
  
  names(beta_v) = names(alpha_v) = var.list
  ### set true average treatment effect
  bA = 0
  
  # set vector of possible lambda's to try
  lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  # gamma_vals <- c(0.5, 1, 2, 3)
  
  ### define some functions for generating data, ATE estimates, and the wAMD,
  expit = plogis
  fam <- binomial()
  ATE_est = function(fY,fw,fA){
    t_ATE = fY*fw
    tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
    return(tt_ATE) 
  }
  create_weights = function(fp,fA){
    fw = (fp)^(-1)
    fw[fA==0] = (1 - fp[fA==0])^(-1)
    return(fw)
  }
  wAMD_function = function(X,A,wgt,beta){
    diff_vec <- sapply(seq.int(ncol(X)), function(jj){
      this_var <- X[,jj] * wgt
      abs((sum( this_var[A==1]) / sum(wgt*A)) - (sum(this_var[A==0]) / sum(wgt*(1-A))))
    })
    wdiff_vec = diff_vec * abs(beta) 
    wAMD = sum(wdiff_vec)
    ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
    return(ret) 
  }
  
  oal_fitting <- function(X, A, Y, betaXY, gamma, lambda1=NULL, lambda2=0) {
    
    stopifnot(length(lambda2) == 1, lambda2 >= 0)
    if(!is.null(lambda1)){
      stopifnot(length(lambda1) == 1)
    }
    if(lambda2 == 0){
      X1 = X
      A1 = A
    } else {
      X1 = rbind(X, sqrt(lambda2)*diag(ncol(X)))
      A1 = c(A, rep(0, ncol(X)))
    }
    pen <- abs(betaXY)^(-gamma)
    logit_oal <- glmnet(X1, A1, family="binomial", intercept=F,
                        standardize = F,
                        penalty.factor=pen)
    
    
    
    
    if(is.null(lambda1)){
      coefs_all <- coef(logit_oal)
      # propensity and coefficient matrix: n x nlambda
      est_propens = predict(logit_oal, type="response", newx=X)
      
      
      wgt_mat <- apply(est_propens, 2, create_weights, fA=Data$A)
      
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      this_wAMD_vec <- apply(wgt_mat, 2, function(wgt){ 
        wAMD_function(X=X, A=A,
                      wgt=wgt, beta=betaXY)$wAMD
      })
      
      this_min_idx = which.min(this_wAMD_vec)
      wAMD = this_wAMD_vec[this_min_idx]
      min_wgts = wgt_mat[,this_min_idx]
      
      min_coefs <- coefs_all[,this_min_idx]
      min_ate = ATE_est(fY=Y,fw=min_wgts,fA=A)
    } else {
      
      coefs_all <- coef(logit_oal, s=mean(pen)*lambda1)
      est_propens = predict(logit_oal, type="response", newx=X, s=mean(pen)*lambda1)
      
      wgt <- create_weights(est_propens, fA=Data$A)
      
      # estimate weighted absolute mean different over all covariates using this lambda to generate weights
      wAMD = wAMD_function(X=X, A=A,
                           wgt=wgt, beta=betaXY)$wAMD
      
      min_coefs <- as.numeric(coefs_all)
      min_ate = ATE_est(fY=Y,fw=wgt,fA=A)
      
    }
    
    ret <- list(ate=min_ate, wAMD=wAMD, coefs=min_coefs)
    return(ret)
  }
  
  
  # ates <- sapply(1L:num_simulations, function(sim_idx){
  ates <- vector("numeric", num_simulations)
  coefs <- matrix(NA, nrow=num_simulations, ncol=1+p)
  wamds <- vector("numeric", num_simulations)
  unpen <- matrix(NA, nrow=num_simulations, ncol=p)
  colnames(coefs) <- var.list.plus.int
  for( sim_idx in seq.int(num_simulations)) {
    
    set.seed(2021 + sim_idx)
    ### simulate data
    Sigma_x = matrix(rho*sig_x^2,nrow=p,ncol=p) 
    diag(Sigma_x) = sig_x^2
    # Sigma_x[1,2] <- Sigma_x[2,1] <- 0.99*sig_x^2 # make the two confounders nearly-collinear
    Mean_x = rep(mean_x,p)
    Data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
    # Data <- scale(Data, T, T)
    names(Data) = var.list
    
    gA_x = as.numeric(as.matrix(Data[,var.list]) %*% alpha_v)
    pA = expit( gA_x )
    summary(pA)
    Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A 
    gY_xA = as.numeric(as.matrix(Data[,var.list]) %*% beta_v)
    Data$Y = gY_xA + rnorm(n=n)
    Data$Y = Data$Y + Data$A*bA
    
    X <- as.matrix(Data[,1:p])
    X <- scale(X, T, T)
    
    # estimate outcome model
    y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
    lm.Y = lm(y.form,data=Data)
    betaXY = coef(lm.Y)[var.list] 
    unpen[sim_idx,] <- betaXY
    
    coeff_XA = (matrix(NA,nrow=1+p,ncol=100))
    rownames(coeff_XA) = var.list.plus.int
    
    ######################################################################################
    #####  Run outcome adaptive lasso for each lambda value 
    ######################################################################################
    
    # alpha_vals <- c(.2,.5,.75,1)
    # alpha_vals <- c(1)
    
    min_wamd <- Inf
    tmp_coefs <- rep(NA, 1+p)
    tmp_ate <- NA
    
    for( lam2_idx in seq_along(lambda2_vals)){
      lambda2 = lambda2_vals[lam2_idx]
      for( gam_idx in seq_along(gamma_vals)){
        gamma = gamma_vals[gam_idx]
        lambda1 = n^(lambda_vec[gam_idx])
        
        this_oal <- oal_fitting(X, Data$A, Data$Y, betaXY, gamma, 
                                lambda1=lambda1*n, lambda2=lambda2)
        
        if(is.null(this_oal))
          next
        
        if(is.na(this_oal$ate))
          next
        
        if(this_oal$wAMD < min_wamd){
          min_wamd <- this_oal$wAMD
          tmp_coefs <- this_oal$coefs
          tmp_ate <- this_oal$ate
        }
      }
    }
    
    wamds[sim_idx] <- min_wamd
    # save coefficients
    coefs[sim_idx,] <- tmp_coefs
    
    # print out ATE corresponding to smallest wAMD value
    ates[sim_idx] <- tmp_ate
    
    if(verbose==2 | (verbose==1 & sim_idx %% 10 == 1))
      message(sprintf("Completed simulation %i.", sim_idx))
  }
  
  # mean(ates)
  # colMeans(coefs[,1:8])
  
  is_selected <- 1*(abs(coefs) > 1e-6)
  
  result=tibble(
    positivity_tol=min( min(pA), min(1-pA) ),
    method=ifelse(use_ridge, "GOAL", "OAL"),
    metrics=list(tibble(wamd=wamds, ate=ates)),
    sel_perf=list(tibble(
      sel_idx=1:p,
      prop_sel=colMeans(is_selected)[-1]
    ))
  )
  
  return(result)
}
