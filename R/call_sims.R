# call_sims.R
library(tibble)
library(tidyr)
library(dplyr)
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
  method=factor(1:4, levels=1:4, 
    labels=c(
      "OAL", "OAL+overlap", "GOALn", "GLM"
    )
  )
)

params <- tidyr::expand_grid(base_params, vary_params)
# params <- params %>%
#   mutate(use_ridge = (method == "GOALn"),
#          use_overlap = (method == "OAL+overlap"),
#          use_glm = (method == "GLM"))


source("./R/simulate_oal.R")

results <- params %>% 
  mutate(
    data=pmap(., simulate_oal)
  )

qsave(results, file="./results/full_results.qs")
