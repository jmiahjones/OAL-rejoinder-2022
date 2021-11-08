# visualize.R
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(qs)

results <- qread("./results/full_results.qs")

plots <-
  results %>% 
  unnest(cols=c(data)) %>% 
  mutate(
    plot=map(sel_perf, ~ ggplot(., aes(x = sel_idx, y = prop_sel)) + geom_line()),
    ate_bias=map_dbl(metrics, ~mean(.$ate)),
    wamd_med=map_dbl(metrics, ~median(.$wamd)),
    title_str=sprintf("%s selection", method),
    sub_str=sprintf("sd(X): %.2f, p/n: %i/%i, rho: %.2f", sig_x, p, n, rho),
    foot_str=sprintf("ATE Bias: %.2f, wAMD Median: %.2f", ate_bias, wamd_med),
    plot=pmap(list(plot, title_str, sub_str, foot_str), ~..1 +
                labs(title=..2, subtitle=..3,
                     x="Covariate", y="Proportion of Time Selected",
                     caption=..4) )
  )

plots %>% 
  select(n:method, ate_bias, wamd_med) %>% View
# d_plots$plot

plots %>% 
  # filter(rho==0.75) %>% View
  mutate(np=factor(sprintf("%i/%i", n, p)),
         rho=factor(rho)) %>% 
  ggplot(aes(x=positivity_tol, y=ate_bias, color=np)) + 
  geom_point(aes(shape=rho), size=3, alpha=0.7) + #scale_shape(solid=F) +
  facet_grid(method~.)


# plot(x=seq.int(p), y=colMeans(is_selected)[-1], type="l", 
#      main=sprintf("%s selection: p/n=%i/%i, rho=%.2f", ifelse(use_ridge, "GOAL", "OAL"), p, n, rho),
#      sub=sprintf("ATE Bias: %.2f", mean(ates)),
#      xlab="Covariate", ylab="Proportion of Time Selected")

# summary(wamds)
# summary(pA)