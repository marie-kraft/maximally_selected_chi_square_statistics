###### Create tables and plots to display results

# load libraries 
library(tidyr)
library(dplyr)
library(ggplot2)

# load results
results_all <- readRDS("./results_all.rds")
# or source("./results.R")

order_of_methods <- c("basic" ,"bonferroni", "k_test", "our", "conditional", "unconditional", "asymptotic")

### First part with equal probs: #####################################
results_equal <- results_all %>%
  filter(is.na(power)) %>%
  filter(method != "our_greater") %>%
  mutate(method = factor(method, levels = order_of_methods))

results_equal$LowerCI <- results_equal$alpha_error - 1.96 * results_equal$monte_carlo_se
results_equal$UpperCI <- results_equal$alpha_error + 1.96 * results_equal$monte_carlo_se

# Format the alpha_error and monte_carlo_se into a single string
results_equal <- results_equal %>%
  mutate(error_combined = paste0(round(alpha_error, digits = 3),
                                 " (",
                                 round(monte_carlo_se, digits = 3),
                                 ")"))

## For n=50 :
results_equal_50 <- results_equal %>%
  filter(nobs == "50, 50, 50, 50")

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_equal_50 %>%
  select(method, error_combined, probability) %>%
  pivot_wider(names_from = probability, values_from = error_combined)# Reorder the columns with method as the first column

result_df <- result_df %>%
  mutate(method = factor(method, levels = order_of_methods)) %>%
  arrange(method)
View(result_df)

# Preparation for the plots
results_equal_50$Method <- factor(results_equal_50$method,
                                  levels = rev(order_of_methods),
                                  labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test", 
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

prob_labels <- c(
  "0.01, 0.01, 0.01, 0.01" = "p = 0.01",
  "0.11, 0.11, 0.11, 0.11" = "p = 0.11",
  "0.21, 0.21, 0.21, 0.21" = "p = 0.21",
  "0.31, 0.31, 0.31, 0.31" = "p = 0.31",
  "0.41, 0.41, 0.41, 0.41" = "p = 0.41",
  "0.51, 0.51, 0.51, 0.51" = "p = 0.51"
)

# Plot:

ggplot(results_equal_50, aes(x = Method, y = alpha_error)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_wrap(~ probability, ncol = 2, labeller = labeller(probability = prob_labels)) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated type-I-error") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

# ratio 11 x 5.5 inches

## For n = 20:
results_equal_20 <- results_equal %>%
  filter(nobs == "20, 20, 20, 20")

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_equal_20 %>%
  select(method, error_combined, probability) %>%
  pivot_wider(names_from = probability, values_from = error_combined)# Reorder the columns with method as the first column

result_df <- result_df %>%
  arrange(method)
View(result_df)

# for plots
results_equal_20$Method <- factor(results_equal_20$method,
                                  levels = rev(order_of_methods),
                                  labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test",
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

prob_labels <- c(
  "0.01, 0.01, 0.01, 0.01" = "p = 0.01",
  "0.11, 0.11, 0.11, 0.11" = "p = 0.11",
  "0.21, 0.21, 0.21, 0.21" = "p = 0.21",
  "0.31, 0.31, 0.31, 0.31" = "p = 0.31",
  "0.41, 0.41, 0.41, 0.41" = "p = 0.41",
  "0.51, 0.51, 0.51, 0.51" = "p = 0.51"
)

ggplot(results_equal_20, aes(x = Method, y = alpha_error)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_wrap(~ probability, ncol = 2, labeller = labeller(probability = prob_labels)) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated type-I-error") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

## For n = 10:
results_equal_10 <- results_equal %>%
  filter(nobs == "10, 10, 10, 10")

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_equal_10 %>%
  select(method, error_combined, probability) %>%
  pivot_wider(names_from = probability, values_from = error_combined)# Reorder the columns with method as the first column

result_df <- result_df %>%
  arrange(method)
View(result_df)

# for plots
results_equal_10$Method <- factor(results_equal_10$method,
                                  levels = rev(order_of_methods),
                                  labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test", 
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

prob_labels <- c(
  "0.01, 0.01, 0.01, 0.01" = "p = 0.01",
  "0.11, 0.11, 0.11, 0.11" = "p = 0.11",
  "0.21, 0.21, 0.21, 0.21" = "p = 0.21",
  "0.31, 0.31, 0.31, 0.31" = "p = 0.31",
  "0.41, 0.41, 0.41, 0.41" = "p = 0.41",
  "0.51, 0.51, 0.51, 0.51" = "p = 0.51"
)

ggplot(results_equal_10, aes(x = Method, y = alpha_error)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_wrap(~ probability, ncol = 2, labeller = labeller(probability = prob_labels)) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated type-I-error") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))


### Second part with equal group size and different probs ##############
results_equal_diff <- results_all %>%
  filter(is.na(alpha_error)) %>%
  filter(method != "our_greater") %>%
  mutate(method = factor(method, levels = order_of_methods))

results_equal_diff <- results_equal_diff[1:63, ]

results_equal_diff$LowerCI <- results_equal_diff$power - 1.96 * results_equal_diff$monte_carlo_se
results_equal_diff$UpperCI <- results_equal_diff$power + 1.96 * results_equal_diff$monte_carlo_se

# Format the power and monte_carlo_se into a single string
results_equal_diff <- results_equal_diff %>%
  mutate(power_combined = paste0(round(power, digits = 3),
                                 " (",
                                 round(monte_carlo_se, digits = 3),
                                 ")"))

results_equal_diff <- results_equal_diff %>%
  mutate(probabilities = ifelse(probability == "0.2, 0.2, 0.2, 0.15",
                                "small",
                                ifelse(probability == "0.3, 0.3, 0.25, 0.15",
                                       "medium",
                                       ifelse(probability == "0.3, 0.3, 0.15, 0.35",
                                              "big", probability))))

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_equal_diff %>%
  select(method, power_combined, probabilities, nobs) %>%
  pivot_wider(names_from = probabilities, values_from = power_combined) %>%
  pivot_wider(names_from = nobs, values_from = c(small, medium, big))
  
  
result_df <- result_df %>%
  mutate(method = factor(method, levels = order_of_methods)) %>%
  arrange(method)
View(result_df)

# for plots
results_equal_diff$Method <- factor(results_equal_diff$method,
                                    levels = rev(order_of_methods),
                                    labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test",  # Correct subscript
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

nobs_order <- c("40, 40, 40, 40", "20, 20, 20, 20", "10, 10, 10, 10")
results_equal_diff$nobs <- factor(results_equal_diff$nobs, levels = nobs_order, labels = c(
  "n = 40", "n = 20", "n = 10"
))

probs_order <- c("0.2, 0.2, 0.2, 0.15", "0.3, 0.3, 0.25, 0.15", "0.3, 0.3, 0.15, 0.35")
results_equal_diff$probability <- factor(results_equal_diff$probability,
                                         levels = probs_order,
                                         labels = c(
                                           "p = (0.2, 0.2, 0.2, 0.15)",
                                           "p = (0.3, 0.3, 0.25, 0.15)",
                                           "p = (0.3, 0.3, 0.15, 0.35)"
                                         ))

ggplot(results_equal_diff, aes(x = Method, y = power)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_grid(nobs ~ probability) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated statistical power") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))
# 13x6

### Third part with different group size and different probs ###########
results_diff_diff <- results_all %>%
  filter(is.na(alpha_error)) %>%
  filter(method != "our_greater") %>%
  mutate(method = factor(method, levels = order_of_methods))

results_diff_diff <- results_diff_diff[64:77, ]

results_diff_diff$LowerCI <- results_diff_diff$power - 1.96 * results_diff_diff$monte_carlo_se
results_diff_diff$UpperCI <- results_diff_diff$power + 1.96 * results_diff_diff$monte_carlo_se

# Format the power and monte_carlo_se into a single string
results_diff_diff <- results_diff_diff %>%
  mutate(power_combined = paste0(round(power, digits = 3),
                                 " (",
                                 round(monte_carlo_se, digits = 3),
                                 ")"))

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_diff_diff %>%
  select(method, power_combined, probability, nobs) %>%
  pivot_wider(names_from = nobs, values_from = c(power_combined))


result_df <- result_df %>%
  mutate(method = factor(method, levels = order_of_methods)) %>%
  arrange(method)
View(result_df)

# for plots
results_diff_diff$Method <- factor(results_diff_diff$method,
                                   levels = rev(order_of_methods),
                                   labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test",  # Correct subscript
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

nobs_order <- c("40, 20, 20, 20", "20, 10, 10, 10")
results_diff_diff$nobs <- factor(results_diff_diff$nobs,
                                 levels = nobs_order,
                                 labels = c(
  "n = (40, 20, 20, 20)", "n = (20, 10, 10, 10)"
))


ggplot(results_diff_diff, aes(x = Method, y = power)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_grid(~ nobs) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated statistical power") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))
# 11 x 4

### Fourth part with different group size and different probs ############
results_size <- results_all %>%
  filter(is.na(alpha_error)) %>%
  filter(method != "our_greater") %>%
  mutate(method = factor(method, levels = order_of_methods))

results_size <- results_size[78:98, ]

results_size$LowerCI <- results_size$power - 1.96 * results_size$monte_carlo_se
results_size$UpperCI <- results_size$power + 1.96 * results_size$monte_carlo_se

# Format the power and monte_carlo_se into a single string
results_size <- results_size %>%
  mutate(power_combined = paste0(round(power, digits = 3),
                                 " (",
                                 round(monte_carlo_se, digits = 3),
                                 ")"))

# Reshape the data so that each probability vector is a column, and methods are rows
result_df <- results_size %>%
  select(method, power_combined, probability, nobs) %>%
  pivot_wider(names_from = probability, values_from = c(power_combined))


result_df <- result_df %>%
  mutate(method = factor(method, levels = order_of_methods)) %>%
  arrange(method)
View(result_df)

# for plots
results_size$Method <- factor(results_size$method,
                              levels = rev(order_of_methods),
                              labels = rev(c(
  "Naive~chi^2-test",
  "Bonferroni-adjusted~chi^2-test",
  "K-group~chi^2-test",
  "Flexible~exact~chi[max]^2-test",  # Correct subscript
  "Conditional~exact~test",
  "Unconditional~exact~test",
  "Asymptotic~Simultaneous~CI"
)))

probs_order <- c("0.3, 0.3, 0.15", "0.3, 0.3, 0.25", "0.2, 0.2, 0.15, 0.25, 0.2, 0.18")
results_size$probability <- factor(results_size$probability, levels = probs_order, labels = c(
  "p = (0.3, 0.3, 0.15)", "p = (0.3, 0.3, 0.25)", "p = (0.2, 0.2, 0.15, 0.25, 0.2, 0.18)"
))

start <- proc.time()
ggplot(results_size, aes(x = Method, y = power)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_segment(aes(x = Method, xend = Method, y = LowerCI, yend = UpperCI),
               color = "gray50", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 4) +
  facet_grid(~ probability) +
  coord_flip() +
  labs(x = "Method",
       y = "Estimated statistical power") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(text = element_text(size = 12))
end <- proc.time()
end - start

# 13x4


