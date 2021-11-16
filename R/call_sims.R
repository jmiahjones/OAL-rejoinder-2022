# call_sims.R
library(tibble)
library(tidyr)
library(purrr)
library(qs)

rho <- c(0, 0.25, 0.5, 0.75)

base_params <- tribble(
  ~n,   ~p,   ~num_simulations,
  200,  100,  100,
  500,  200,  100
)

vary_params <- tidyr::expand_grid(rho=c(0, 0.25, 0.5, 0.75), 
  sig_x=seq(.2, 1, by=.2), 
  scenario=1, 
  use_ridge=c(F,T))

params <- tidyr::expand_grid(base_params, vary_params)
stopifnot(
  nrow(params) == (nrow(base_params) * nrow(vary_params)),
  ncol(params) == (ncol(base_params) + ncol(vary_params))
)

# params <- tribble(
#   ~n, ~p, ~num_simulations, ~rho, ~sig_x, ~scenario, ~use_ridge,
#   200, 100, 100, 0.00, 1/4, 1, F,
#   200, 100, 100, 0.00, 1/4, 1, T,
#   200, 100, 100, 0.50, 1/4, 1, F,
#   200, 100, 100, 0.50, 1/4, 1, T,
#   200, 100, 100, 0.75, 1/4, 1, F,
#   200, 100, 100, 0.75, 1/4, 1, T,
  
#   200, 100, 100, 0.00,   1, 1, F,
#   200, 100, 100, 0.75,   1, 1, F,
#   200, 100, 100, 0.00,   1, 1, T,
#   200, 100, 100, 0.75,   1, 1, T,
  
#   500, 200, 100, 0.00, 1/4, 1, F,
#   500, 200, 100, 0.00, 1/4, 1, T,
#   500, 200, 100, 0.50, 1/4, 1, F,
#   500, 200, 100, 0.50, 1/4, 1, T,
#   500, 200, 100, 0.75, 1/4, 1, F,
#   500, 200, 100, 0.75, 1/4, 1, T
  
# )

source("./R/simulate_oal.R")

# foo = simulate_oal(200, 100, 2, 0.00, 1/4, 1, F, verbose=2)

results <- params %>% 
  mutate(
    data=pmap(., simulate_oal)
  )

qsave(results, file="./results/full_results.qs")
