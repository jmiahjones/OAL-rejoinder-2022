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
    # lam2_med=map_dbl(metrics, ~quantile(.$lambda2, probs=.5)),
    lam2_med=map_dbl(metrics, ~median(.$lambda2)),
    title_str=sprintf("%s selection", method),
    sub_str=sprintf("sd(X): %.2f, p/n: %i/%i, rho: %.2f", sig_x, p, n, rho),
    foot_str=sprintf("ATE Bias: %.2f, wAMD Median: %.2f", ate_bias, wamd_med),
    plot=pmap(list(plot, title_str, sub_str, foot_str), ~..1 +
                labs(title=..2, subtitle=..3,
                     x="Covariate", y="Proportion of Time Selected",
                     caption=..4) )
  )



# plots %>% 
#   select(n:method, ate_bias:lam2_med) %>% View


# png("./plots/bias-by-pos.png")
# plots %>% 
#   # filter(rho==0.75) %>% View
#   mutate(np=factor(sprintf("%i/%i", n, p)),
#          rho=factor(rho)) %>% 
#   ggplot(aes(x=positivity_tol, y=ate_bias, color=np)) + 
#   geom_point(aes(shape=rho), size=3, alpha=0.7) + #scale_shape(solid=F) +
#   facet_grid(method~.) + 
#   labs(
#     title = "ATE Bias vs. Positiviy Violations",
#     y = "ATE Bias",
#     x = expression(paste("Min of ", pi, " and ", 1 - pi))
#   )
# dev.off()

# png("./plots/pos-by-sd.png")
# plots %>% 
#   # filter(rho==0.75) %>% View
#   mutate(np=factor(sprintf("%i/%i", n, p)),
#          rho=factor(rho),
#          sig_x=factor(sig_x)) %>% 
#   ggplot(aes(x=sig_x, y=positivity_tol, color=rho)) + 
#   geom_point(aes(shape=rho), size=3, alpha=0.5) + #scale_shape(solid=F) +
#   facet_grid(method~.) + 
#   labs(
#     title = "Positiviy Violations by X Std. Dev.",
#     x = "sd(X)",
#     y = expression(paste("Min of ", pi, " and ", 1 - pi))
#   )
# dev.off()

png("./plots/ate_bias.png", width = 700, height = 1000)
plots %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(sprintf("\u03c1 = %.2f", rho))
  ) %>%
  ggplot(aes(x = sig_x, y = ate_bias)) +
  geom_point(aes(shape = method, group = method),
    size = 3,
    position = position_dodge(width = 0.1)
  ) +
  # geom_line(aes(linetype=method), position="dodge") +
  scale_x_continuous(breaks = seq(.2, 1, by = .2)) +
  facet_grid(rho ~ np) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(
    title = "ATE Bias vs. sd(X)",
    y = "ATE Bias",
    x = "sd(X)", shape = "Method"
  )
dev.off()


png("./plots/positivity_violations.png", width = 700, height = 700)
plots %>%
  filter(use_ridge == F) %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    # rho = factor(sprintf("\u03c1 = %.2f", rho))
    rho = factor(rho)
  ) %>%
  ggplot(aes(x = sig_x, y = prop_positivity)) +
  geom_line(aes(linetype = rho)) +
  facet_grid(~np) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  labs(
    title = "Positivity Violations vs. sd(X)",
    y = "Pr{ min(\u03c0, 1-\u03c0) < .05}",
    x = "sd(X)", linetype = "\u03c1"
  )
dev.off()


# png("./plots/tmp.png")
# plots %>% 
#   mutate(np=factor(sprintf("n/p = %i/%i", n, p)),
#          rho=factor(sprintf("\u03c1 = %.2f", rho))) %>% 
#   ggplot(aes(x=sig_x, y=ate_bias)) + 
#   # geom_point(aes(shape=method), size=3, alpha=0.7) +
#   geom_line(aes(linetype=method)) +
#   facet_grid(rho~np) + 
#   labs(
#     title = "ATE Bias vs. sd(X)",
#     y = "ATE Bias",
#     x = expression(paste("Min of ", pi, " and ", 1 - pi))
#   )
# dev.off()
