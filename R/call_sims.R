# call_sims.R
library(tibble)
library(purrr)
library(qs)
# 
# params <- tribble(
#   ~n, ~p, ~num_simulations, ~rho, ~sig_x, ~scenario, ~use_ridge,
#   200, 100, 100, 0.00, 1/4, 1, F,
#   200, 100, 100, 0.00, 1/4, 1, T,
#   200, 100, 100, 0.50, 1/4, 1, F,
#   200, 100, 100, 0.50, 1/4, 1, T,
#   200, 100, 100, 0.75, 1/4, 1, F,
#   200, 100, 100, 0.75, 1/4, 1, T,
#   
#   200, 100, 100, 0.00,   1, 1, F,
#   200, 100, 100, 0.75,   1, 1, F,
#   
#   500, 200, 100, 0.00, 1/4, 1, F,
#   500, 200, 100, 0.00, 1/4, 1, T,
#   500, 200, 100, 0.50, 1/4, 1, F,
#   500, 200, 100, 0.50, 1/4, 1, T,
#   500, 200, 100, 0.75, 1/4, 1, F,
#   500, 200, 100, 0.75, 1/4, 1, T
#   
# )
# 
# source("./R/simulate_oal.R")
# 
# results <- params %>% 
#   mutate(
#     data=pmap(., simulate_oal)
#   )
# 
# qsave(results, file="./results/full_results.qs")


## Add additional sims

params <- tribble(
  ~n, ~p, ~num_simulations, ~rho, ~sig_x, ~scenario, ~use_ridge,
  
  200, 100, 100, 0.00,   1, 1, T,
  200, 100, 100, 0.75,   1, 1, T
  
)

source("./R/simulate_oal.R")

results <- params %>% 
  mutate(
    data=pmap(., simulate_oal)
  )

foo <- simulate_oal(200,100,100,0.00,1,1,T)

qsave(results, file="./results/add_results.qs")

if(file.exists("./results/add_results.qs")){
  if(file.exists("./results/full_results.qs")){
    full <- qread("./results/full_results.qs")
    join <- rbind(full, results)
    qsave(join, "./results/join_results.qs")
  }
}
