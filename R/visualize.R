# visualize.R
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(qs)

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
  mutate(np = factor(sprintf("n/p = %i/%i", n, p)),
    weight_type = factor(
      weight_type, levels = c("notrunc", "trunc@.01", "trunc@.05", "overlap"), 
      labels = c("No Truncation", "Truncated @ .01", "Truncated @ .05", "Overlap")
    ),
    est_method = factor(est_method, levels = c("IPW", "AIPW", "TMLE"))
  )


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
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "vertical") +
  labs(
    title = "ATE Bias vs. Estimator",
    subtitle = "Propensities from OAL and GOALn",
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
    subtitle = "Propensities from OAL and GOALn",
    # caption = "Scenario 1 with sd(X)=1 and \u03c1=0.75",
    y = "ATE RMSE",
    x = "Estimator", shape = "Weighting Method", 
    color = "Weighting Method"
  )
dev.off()
