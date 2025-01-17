#### Create the results from the lists

# load libraries
library(tidyverse)
library(rsimsum)
library(dplyr)

# Label method variables

method_names <- c("our",
                  "basic",
                  "bonferroni",
                  "asymptotic",
                  "k_test",
                  "unconditional",
                  "conditional")


# Function to process each entry
get_results <- function(list1, list2, method_names, equal_probs = TRUE) {
    # initialize first dataframe
  results1 <- data.frame(
    repetition = integer(),
    probability = character(),
    nobs = character(),
    method = character(),
    p_value = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
    )
  
  
  # get repetition, probability and group size measurements for each scenario
  results1 <- do.call(rbind, lapply(list1, function(entry) {
    repetition <- entry$repetition
    probability <- paste(entry$probabilities, collapse = ", ")
    nobs <- paste(entry$nobs, collapse = ", ")
  
    # get the calculates values for each approach
    data.frame(
      repetition = repetition,
      probability = paste(probability, sep = " ,"),
      nobs = paste(nobs, sep = " ,"),
      method = method_names[-7],
      p_value = sapply(method_names[-7], function(m) entry[[paste0(m, "_p_value")]]),
      decision = sapply(method_names[-7], function (m) entry[[paste0(m, "_decision")]]),
     stringsAsFactors = FALSE
    )
    }
    ))
  # rejection indicator
  results1$reject_flag <- ifelse(results1$decision == "reject", 1, 0)

  # if second list (conditional):
  if (!missing(list2)) {
    results2 <- data.frame(
      repetition = integer(),
      probability = character(),
      nobs = character(),
      method = character(),
      p_value = numeric(),
      decision = character(),
      stringsAsFactors = FALSE
      )
  
    results2 <- do.call(rbind, lapply(list2, function(entry) {
      repetition <- entry$repetition
      probability <- paste(entry$probabilities, collapse = ", ")
      nobs <- paste(entry$nobs, collapse = ", ")
  
      data.frame(
        repetition = repetition,
        probability = paste(probability, sep = " ,"),
        nobs = paste(nobs, sep = " ,"),
        method = method_names[7],
        p_value = sapply(method_names[7], function(m) entry[[paste0(m, "_p_value")]]),
        decision = sapply(method_names[7], function (m) entry[[paste0(m, "_decision")]]),
        stringsAsFactors = FALSE
        )
      }
      ))
    
    results2$reject_flag <- ifelse(results2$decision == "reject", 1, 0)

    results <- rbind(results1, results2)

    } else {
      results <- results1
      }

  # get monte carlo SE's for all scenarios
  monte_carlo_results <- aggregate(
    reject_flag ~ probability + method + nobs,
    data = results,
    FUN = function(x) sqrt(((sum(x, na.rm = TRUE) / length(na.omit(x))) * (1 - (sum(x, na.rm = TRUE) / length(na.omit(x))))) / length(na.omit(x)))
    )

  colnames(monte_carlo_results)[colnames(monte_carlo_results) == "reject_flag"] <- "monte_carlo_se"
  
  # get number for seeing if NA values occurred
  full_counts <- aggregate(
    reject_flag ~ probability + method +nobs,
    data = results,
    FUN = function(x) length(na.omit(x))
  )
  colnames(full_counts)[colnames(full_counts) == "reject_flag"] <- "full_counts"

  if (equal_probs) { # equal probabilities are true, calculate alpha_error
    alpha_error_results <- aggregate(
      reject_flag ~ probability + method + nobs,
      data = results,
      FUN = function(x) sum(x, na.rm = TRUE) / length(na.omit(x))
      )
    
    colnames(alpha_error_results)[colnames(alpha_error_results) == "reject_flag"] <- "alpha_error"
    
    estimates <- full_join(alpha_error_results, monte_carlo_results, by = c("probability", "method", "nobs"))
    estimates <- estimates <- full_join(estimates, full_counts, by = c("probability", "method", "nobs"))

    } else { # if different probabilities, calculate power
      power_results <- aggregate(
        reject_flag ~ probability + method + nobs,
        data = results,
        FUN = function(x) (1 - (sum(x, na.rm = TRUE) / length(na.omit(x))))
        )
      
      colnames(power_results)[colnames(power_results) == "reject_flag"] <- "power"
      
      estimates <- full_join(power_results,
                             monte_carlo_results,
                             full_counts,
                             by = c("probability", "method", "nobs"))
      
      estimates <- estimates <- full_join(estimates, full_counts, by = c("probability", "method", "nobs"))
      }
  
  return(estimates)
  
  }


### load the lists, get the results and delete the lists

list_names <- data.frame(normal_lists = c("list_3_no",
                                          "list_3_small",
                                          "list_3_small_equal",
                                          "list_3_big_small",
                                          "list_3_medium_small",
                                          "list_3_small_small",
                                          "list_3_big_big",
                                          "list_3_medium_big",
                                          "list_3_small_big",
                                          "list_3_big_diff",
                                          "list_3_medium_diff",
                                          "list_3_small_diff",
                                          "list_3_big_size",
                                          "list_3_small_size",
                                          "list_2_small_big",
                                          "list_2_small_small",
                                          "list_5_small_small"),
                         conditional_lists = c("conditional_list_3_no",
                                               "conditional_list_3_small",
                                               "conditional_list_3_small_equal",
                                               "conditional_list_3_big_small",
                                               "conditional_list_3_medium_small",
                                               "conditional_list_3_small_small",
                                               "conditional_list_3_big_big",
                                               "conditional_list_3_medium_big",
                                               "conditional_list_3_small_big",
                                               "conditional_list_3_big_diff",
                                               "conditional_list_3_medium_diff",
                                               "conditional_list_3_small_diff",
                                               "conditional_list_3_big_size",
                                               "conditional_list_3_small_size",
                                               "conditional_list_2_small_big",
                                               "conditional_list_2_small_small",
                                               "conditional_list_5_small_small"),
                         h0_true = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                     FALSE, FALSE, FALSE, FALSE))

results_all <- data.frame(probability = character(),
                          method = character(),
                          nobs = character(),
                          power = numeric(),
                          alpha_error = numeric(),
                          monte_carlo_se = numeric())

for (i in 1:(nrow(list_names))) {
  normal_list <- readRDS(paste0("./Simulation_results/", list_names[i, 1], ".rds"))
  conditional_list <- readRDS(paste0("./Simulation_results/", list_names[i, 2], ".rds"))
  
  results_temp <- get_results(normal_list, conditional_list, method_names = method_names,
                              equal_probs = list_names[i, 3])
  
  if ("power" %in% colnames(results_temp)) {
    results_temp$alpha_error <- NA
  } else if ("alpha_error" %in% colnames(results_temp)) {
    results_temp$power <- NA
  }
  
  results_all <- rbind(results_all, results_temp)
  
  rm(normal_list)
  rm(conditional_list)
}

# save results
saveRDS(results_all, file = "./results_all.rds")



