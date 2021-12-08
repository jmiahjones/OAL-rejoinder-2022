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
    method=factor(if_else(
      use_ridge, "GOALn",
      if_else(use_overlap, "OAL+overlap", "OAL")
    )),
    plot=map(sel_perf, ~ ggplot(., aes(x = sel_idx, y = prop_sel)) + geom_line()),
    ate_bias=map_dbl(metrics, ~abs(mean(.$ate))),
    ate_mse=map_dbl(metrics, ~mean(.$ate^2)),
    wamd_med=map_dbl(metrics, ~median(.$wamd)),
    lam2_med=map_dbl(metrics, ~median(.$lambda2)),
    title_str=sprintf("%s selection", method),
    sub_str=sprintf("sd(X): %.2f, p/n: %i/%i, rho: %.2f", sig_x, p, n, rho),
    foot_str=sprintf("ATE Bias: %.2f, wAMD Median: %.2f", ate_bias, wamd_med),
    plot=pmap(list(plot, title_str, sub_str, foot_str), ~..1 +
                labs(title=..2, subtitle=..3,
                     x="Covariate", y="Proportion of Time Selected",
                     caption=..4) +
                theme_bw() +
                theme(text = element_text(size = 20))
    )
  )


png("./plots/replication.png")
plots %>%
  filter(sig_x == 1) %>% 
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(sprintf("%.2f", rho))
  ) %>%
  ggplot(aes(x = rho, y = ate_bias)) +
  geom_point(aes(shape = method, group = method, color = method),
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  facet_grid(cols=vars(np)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  labs(
    title = "ATE Bias vs. \u03c1",
    caption = "Scenario 1 with sd(X)=1",
    y = "ATE Bias",
    x = "\u03c1", shape = "Method", color = "Method"
  )
dev.off()

png("./plots/ate_bias.png", width = 700, height = 1000)
plots %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(sprintf("\u03c1 = %.2f", rho))
  ) %>%
  ggplot(aes(x = sig_x, y = ate_bias)) +
  geom_point(aes(shape = method, group = method, color = method),
    size = 3,
    position = position_dodge(width = 0.1)
  ) +
  scale_x_continuous(breaks = seq(.2, 1, by = .2)) +
  facet_grid(rho ~ np) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  labs(
    title = "ATE Bias vs. sd(X)",
    y = "ATE Bias",
    x = "sd(X)", shape = "Method", color = "Method"
  )
dev.off()

png("./plots/ate_mse.png", width = 700, height = 1000)
plots %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(sprintf("\u03c1 = %.2f", rho))
  ) %>%
  ggplot(aes(x = sig_x, y = ate_mse)) +
  geom_point(aes(shape = method, group = method, color = method),
    size = 3,
    position = position_dodge(width = 0.1)
  ) +
  scale_x_continuous(breaks = seq(.2, 1, by = .2)) +
  facet_grid(rho ~ np) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  labs(
    title = "ATE MSE vs. sd(X)",
    y = "ATE MSE",
    x = "sd(X)", shape = "Method", color = "Method"
  )
dev.off()


png("./plots/positivity_violations.png", width = 700, height = 700)
plots %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
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

selection_plots <- results %>% 
  unnest(cols = c(data)) %>% 
  mutate(
    method = factor(if_else(
      use_ridge, "GOALn",
      if_else(use_overlap, "OAL+overlap", "OAL")
    ))
  ) %>%
  unnest(cols = sel_perf)


png("./plots/selection-rho-0.png", width = 700, height = 700)
selection_plots %>% 
  filter(
    method != "OAL+overlap",
    abs(sig_x - 0.4) < 1e-3, rho==0.0, n==200
  ) %>% 
  ggplot(aes(x = sel_idx, y = prop_sel)) + 
  geom_line(aes(
    color = method, linetype = method,
    group = method
  )) +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = "bottom") +
  labs(
    title = "OAL vs GOALn Selection",
    subtitle = "\u03c1: 0.0",
    caption = "sd(X): 0.4, p/n: 100/200",
    x="Covariate", y="Proportion of Time Selected",
    color = "Method", linetype = "Method"
  )
dev.off()


png("./plots/selection-rho-075.png", width = 700, height = 700)
selection_plots %>% 
  filter(
    method != "OAL+overlap",
    abs(sig_x - 0.4) < 1e-3, rho==0.75, n==200
  ) %>% 
  ggplot(aes(x = sel_idx, y = prop_sel)) + 
  geom_line(aes(
    color = method, linetype = method,
    group = method
  )) +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = "bottom") +
  labs(
    title = "OAL vs GOALn Selection",
    subtitle = "\u03c1: 0.75",
    caption = "sd(X): 0.4, p/n: 100/200",
    x="Covariate", y="Proportion of Time Selected",
    color = "Method", linetype = "Method"
  )
dev.off()
