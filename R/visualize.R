# visualize.R
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(qs)

sig_x_str <- "\u03c3\u2093"
GT_ZERO_CONST <- 1e-5

results <- qread("./results/extensive_results.qs")

validation <- results %>% 
  unnest(data) %>%
  unnest(method_results) %>%
  unnest(data) %>%
  group_by(across(c(n:sel_method, est_method, weight_type))) %>%
  summarize(count = n(), .groups = "drop") %>%
  filter(count != num_simulations) %>%
  nrow() > 1

if(validation)
  stop("Error in grouping! Double-check.")

ate_perf <- results %>% 
  unnest(data) %>%
  unnest(method_results) %>%
  unnest(data) %>%
  group_by(across(c(n:sel_method, est_method, weight_type))) %>%
  summarize(
    ate_bias = abs(mean(ates)),
    ate_rmse = sqrt(mean(ates^2)),
    .groups = "drop"
  ) %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    weight_type = factor(
      weight_type, levels = c("notrunc", "trunc@.01", "trunc@.05", "overlap"), 
      labels = c("Inverse Propensity", "Truncated @ .01", "Truncated @ .05", "Overlap")
    ),
    est_method = factor(est_method, levels = c("IPW", "AIPW", "TMLE"))
  )

# Truncation performs poorly and clutters the plot. Remove.
ate_perf <- ate_perf %>%
  filter(weight_type %in% c("Inverse Propensity", "Overlap"))

png("./plots/extended_ate_bias.png", width = 700, height = 800)
ate_perf %>% 
  ggplot(
    aes(x = est_method, y = ate_bias, 
      color = weight_type, shape = weight_type
    )
  ) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(np), cols = vars(sel_method)) +
  theme_bw() +
  theme(text = element_text(size = 26),
        # legend.position = "bottom",
        legend.position = "none"
        # ,
        # legend.justification = "center",
        # legend.direction = "vertical"
        ) +
  labs(
    title = "ATE Bias vs. Estimator",
    # subtitle = "Propensities from OAL and GOALn",
    # caption = "Scenario 1 with sd(X)=1 and \u03c1=0.75",
    y = "ATE Bias",
    x = "Estimator", shape = "Weighting Method", 
    color = "Weighting Method"
  )
dev.off()


png("./plots/extended_ate_rmse.png", width = 700, height = 800)
ate_perf %>% 
  ggplot(
    aes(x = est_method, y = ate_rmse, 
      color = weight_type, shape = weight_type
    )
  ) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(np), cols = vars(sel_method)) +
  theme_bw() +
  theme(text = element_text(size = 26),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "vertical") +
  labs(
    title = "ATE RMSE vs. Estimator",
    # subtitle = "Propensities from OAL and GOALn",
    # caption = "Scenario 1 with sd(X)=1 and \u03c1=0.75",
    y = "ATE RMSE",
    x = "Estimator", shape = "Weighting Method", 
    color = "Weighting Method"
  )
dev.off()


############################
# Positivity vs Bias
############################
results <- qread("./results/full_results.qs")

png("./plots/positivity_violations.png", width = 700, height = 700)
results %>%
  filter(sel_method == "OAL") %>%
  unnest(data) %>%
  group_by(n,p,num_simulations,rho,sig_x,scenario) %>%
  summarize(prop_positivity = mean(prop_positivity), .groups = "drop") %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(rho)
  ) %>%
  ggplot(aes(x = sig_x, y = prop_positivity)) +
  geom_line(aes(linetype = rho)) +
  facet_grid(np ~.) +
  theme_bw() +
  theme(text = element_text(size = 26),
        legend.position = "bottom",
        legend.justification = "center") +
  labs(
    title = paste0("Positivity Violations vs. ", sig_x_str),
    y = "Pr{ min(\u03c0, 1-\u03c0) < .05}",
    x = sig_x_str, linetype = "\u03c1"
  )
dev.off()

ate_perf <- results %>% 
  unnest(data) %>%
  unnest(method_results) %>%
  unnest(data) %>%
  group_by(across(c(n:sel_method, est_method, weight_type))) %>%
  summarize(
    ate_bias = abs(mean(ates)),
    ate_rmse = sqrt(mean(ates^2)),
    .groups = "drop"
  ) %>%
  filter(weight_type == "notrunc", est_method == "IPW") %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p))
  )

png("./plots/ate_bias.png", width = 700, height = 700)
ate_perf %>%
  filter(rho == 0.75) %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(sprintf("\u03c1 = %.2f", rho))
  ) %>%
  ggplot(aes(x = sig_x, y = ate_bias)) +
  geom_point(aes(shape = sel_method, group = sel_method, color = sel_method),
    size = 3,
    position = position_dodge(width = 0.1)
  ) +
  scale_x_continuous(breaks = seq(.2, 1, by = .2)) +
  # facet_grid(rho ~ np) +
  facet_grid(np ~.) +
  theme_bw() +
  theme(text = element_text(size = 26),
        legend.position = "bottom",
        legend.justification = "center") +
  labs(
    title = paste0("Bias of the Standard IPW Estimator vs. ", sig_x_str),
    subtitle = "\u03c1 = 0.75",
    y = "Bias",
    x = sig_x_str, shape = "Method", color = "Method"
  )
dev.off()



#########################################################################
#### Selection Frequency
#########################################################################
png("./plots/selection.png", width = 700, height = 700)
results %>%
  filter(n == 200, rho == 0.75, sig_x == 1) %>%
  mutate(
    np = factor(sprintf("n/p = %i/%i", n, p)),
    rho = factor(rho)
  ) %>%
  unnest(data) %>%
  unnest(method_results) %>%
  filter(weight_type == "notrunc") %>%
  unnest(sel_result) %>%
  mutate(prop_sel = 1 * (abs(coefs) > GT_ZERO_CONST)) %>%
  select(np, scenario, sel_method, sel_idx, prop_sel) %>%
  group_by(np, scenario, sel_method, sel_idx) %>%
  summarize(prop_sel = mean(prop_sel)) %>%
  ggplot(
    aes(
      x = sel_idx, y = prop_sel,
      group = sel_method, linetype = sel_method, color = sel_method
    )
  ) +
  geom_vline(xintercept = 6L, size = 0.2) +
  geom_line(size=1) +
  scale_linetype_manual(values = 2:3) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = "bottom"
  ) +
  labs(
    title = "Covariate Selection Frequency",
    # subtitle = "n/p=200/100, rho=0.75, sd(X)=1",
    x = "Covariate", y = "Proportion of Time Selected",
    color = "Method", linetype = "Method"
  )
dev.off()  
