# call_sims.R
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(qs)

base_params <- tribble(
  ~n,   ~p,   ~num_simulations,
  200,  100,  100,
  500,  200,  100
)

vary_params <- tidyr::expand_grid(
  rho = 0.75,
  sig_x = 1,
  scenario = 1, 
  sel_method = factor(1:2, levels = 1:2, 
    labels = c("OAL", "GOALn"))
)

params <- tidyr::expand_grid(base_params, vary_params)

source("./R/simulate_oal.R")

results <- params %>% 
  mutate(
    data = pmap(., simulate_oal)
  )

qsave(results, file="./results/extensive_results.qs")
