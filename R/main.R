# main.R
# creates loops, etc.
library(MASS)
library(glmnet)

num_simulations = 1000L

scenarios <- expand.grid(
  rho=c(0,0.2,0.5,0.75),
  size_label=c("a", "b"),
  scenario=c(1,2)
)
num_scenarios <- nrow(scenarios)

for(scenario_idx in 1:num_scenarios) {
  
  rho <- scenarios$rho[scenario_idx]
  scenario <- scenarios$scenario[scenario_idx]
  size_label <- scenarios$size_label[scenario_idx]
  
  ate_save_file <- sprintf(
    "./results/rho-%.2f-size-%s-scen-%i.qs",
    rho, size_label, scenario
  )
  
  n <- switch (size_label,
               "a" = 200,
               "b" = 500
  )
  
  p <- switch (size_label,
               "a" = 100,
               "b" = 200
  )
  
  # Set strength of relationship between covariates and outcome
  beta_v <- c( 0.6, 0.6, 0.6, 0.6, 0, 0, rep(0,p-6) )
  
  # Set strength of relationship between covariates and treatment
  alpha_v <- switch(scenario,
                    "1" = c( 1.0, 1.0,   0,   0, 1, 1,  rep(0,p-6) ),
                    "2" = c( 0.4, 0.4,   0,   0, 1, 1,  rep(0,p-6) )
  )
  
  ### set true average treatment effect
  bA = 0
  
  
  start <- Sys.time()
  
  # set information for simulating coviariates
  mean_x = 0 
  sig_x = 1/8
  pC = pI = pP = 2
  pS = p - (pC+pI+pP)
  var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
  names(beta_v) = names(alpha_v) = var.list
  
  
  Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)) 
  diag(Sigma_x) = sig_x^2
  
  # set vector of possible lambda's to try
  lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  source(file="./R/run-simulations.R", echo=T)
  
  qs::qsave(ates, file=ate_save_file)
  
  stop <- Sys.time()
  
  difft <- difftime(stop, start, units="mins")
  
  message(sprintf("Completed scenario %i of %i in %.2f mins.", 
                  scenario_idx, num_scenarios,
                  difft))
  
}