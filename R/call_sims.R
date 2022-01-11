# call_sims.R
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(qs)

base_params <- tribble(
  ~n,   ~p,   ~num_simulations,
  200,  100,  1000,
  500,  200,  1000
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

# Extensive results: compare multiple estimators in one pathological setting
# results <- params %>% 
#   mutate(
#     data = pmap(., simulate_oal)
#   )

# qsave(results, file = "./results/2_extensive_results.qs")

# Full results: look at multiple scenarios with only the vanilla IPW
vary_params <- tidyr::expand_grid(
  rho = c(0, 0.25, 0.5, 0.75),
  sig_x = seq(.2, 1, by=.2),
  scenario = 1, 
  sel_method = factor(1:2, levels = 1:2, 
    labels = c("OAL", "GOALn"))
)

params <- tidyr::expand_grid(base_params, vary_params)
stopifnot(
  nrow(params) == (nrow(base_params) * nrow(vary_params)),
  ncol(params) == (ncol(base_params) + ncol(vary_params))
)

# call simulate_oal with est_method = "IPW" and weight_types = "notrunc"
full_results <- params %>% 
  mutate(
    data = pmap(., simulate_oal,
      est_method = "IPW", weight_types = "notrunc"
    )
  )

qsave(results, file = "./results/2_full_results.qs")