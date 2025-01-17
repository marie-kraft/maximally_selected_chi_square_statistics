### Maximally selected chi square statistics, a exact flexible approach

# load libraries
library(combinat)
library(checkmate)

#' Maximally Selected Chi-Square Test
#'
#' Performs a maximally selected Chi-square test on a 2xK contingency table to detect significant group differences.  
#' The function allows custom design matrices and supports one-sided and two-sided alternatives.  
#' 
#' @param data A 2xK numeric matrix with counts for two categories across K groups.  
#' @param level Significance level for hypothesis testing (default = 0.05).  
#' @param trace Logical; if TRUE, prints progress and results (default = TRUE).  
#' @param designmatrix Character or matrix defining group comparisons (default = "groupwise").  
#' @param alternative Specifies the alternative hypothesis: "two.sided", "less", or "greater".  
#' @return A named numeric vector with the computed p-value and standard Chi-square p-value.  
#' @examples  
#' example_dataset <- matrix(c(11, 4, 2, 13, 1, 14, 6, 9), nrow = 2)  
#' calculate_exact_uncond_p_value(example_dataset) 
#' @export  

max.chisq.test <- function(data, level = 0.05, trace = TRUE, designmatrix = "groupwise", alternative = "two.sided") {
  # Check if input is correct
  assert_true(nrow(data) == 2 & (ncol(data) >= 2))
  assertNumeric(data, lower = 0, upper = 200, any.missing = FALSE)
  assert_double(level, lower = 0, upper = 1)
  assert_flag(trace)
  if (is.character(designmatrix)) {
    assert_character(designmatrix, fixed = "groupwise")
  } else {
    assert_true(ncol(designmatrix) == ncol(data))
    assert_true(nrow(designmatrix) >= 1)
    assertNumeric(designmatrix, lower = 0, upper = 2, any.missing = FALSE)
  }
  assert_true(alternative %in% c("two.sided", "less", "equal"))
  
  # Input Setup
  nevents <- sum(data[1, ])
  N_vec <- apply(data, MARGIN = 2, FUN = sum)
  k <- length(N_vec)
  
  # If there are no events or no people without events, return p-value 1 
  if (nevents == 0) {
    return(NA)
  } else if (sum(data[2, ]) == 0) {
    return(NA)
  }
  
  # Define Design Matrix if Needed
  if (all(designmatrix == "groupwise")) {
    designmatrix <- matrix(data = NA, nrow = k - 1, ncol = k)
    designmatrix[, 1] <- 1
    for (i in 2:k) {
      designmatrix[i - 1, i] <- 2
    }
  }
  
  # Change to better version
  test_design <- apply(designmatrix, MARGIN = 1, FUN = function(x) {
    c(which(x == 1), which(x == 2))
  })
  
  # Critical Values and p-Values from Chi-Square Test
  if (k == 2) {
    critical_value <- chisq.test(data, correct = FALSE)$statistic
    normal_p_value <- chisq.test(data, correct = FALSE)$p.value
  } else {
    chi <- numeric()
    normal_p_values <- numeric()
    for (i in 1:ncol(test_design)) {
      test_object <- chisq.test(data[, test_design[, i]], correct = FALSE)
      chi <- c(chi, test_object$statistic)
      normal_p_values <- c(normal_p_values, test_object$p.value)
    }
    critical_value <- max(chi, na.rm = TRUE)
    normal_p_value <- normal_p_values[which.max(chi)]
  }
  
  if (is.na(critical_value)) {
    if (trace) {
      cat("There is No event or no observation without an event. Impossible to compute Chi^2 value")
    }
    return(NA)
  }
  
  # Define Helper Function to Create Paths and Compute Max Chi-Square
  create_max_and_num_path <- function(x, N_vec, nevents, k, crit_values, test_design, direction) {
    if (any(x > N_vec) | sum(x) != nevents)
      return(0)
    else {
      test_mat <- N_vec - x # create a counter matrix with y = 0 entries
      x_mat <- rbind(x, test_mat)
      
      # initialize chi values
      chi_list <- list()
      
      # create transformed chi value
      for (i in 1:ncol(test_design)) {
        test <- x_mat[, test_design[, i]]
        chi_list[[i]] <- abs(test[1, 1] * test[2, 2] - test[1, 2] * test[2, 1]) / sqrt(sum(test[1, ]) * sum(test[2, ]))
      }
      # get max as critical value
      actual_chi_max <- max(unlist(chi_list), na.rm = TRUE)
      
      # compute path possibilities if satisfied
      if (all(N_vec == N_vec[1])) {
        if (!(length(actual_chi_max) < 1) && (!is.na(crit_values)) && (!is.na(actual_chi_max)) && get(direction)(actual_chi_max, crit_values) && !all(is.na(unlist(chi_list)))) {
          return(prod(choose(N_vec, x[1:k])))
        } else {
          return(0)
        }
      } else {
        if (!is.na(actual_chi_max) && get(direction)(actual_chi_max, crit_values[[which.max(unlist(chi_list))]])) {
          return(prod(choose(N_vec, x[1:k])))
        } else {
          return(0)
        }
      }
      
    }
  }
  
  # calculate transformed observed chi
  if (all(N_vec == N_vec[1])) {
    crit_value <- sqrt((critical_value * N_vec[1] * N_vec[1]) / (N_vec[1] + N_vec[1]))
  } else {
    crit_value <- list()
    for (i in 1:ncol(test_design)) {
      crit_value[[i]] <- sqrt((critical_value * N_vec[test_design[1, i]] * N_vec[test_design[2, i]]) / (N_vec[test_design[1, i]] + N_vec[test_design[2, i]]))
    }
  }
  
  # create all possible divisions
  a_mat <- xsimplex(k, nevents)
  
  # get p_values by summing over all path possibilities, where equation is
  # satisfied, dividing by all possible divisions
  if (alternative == "greater") {
    test.direction <- "<="
    p_value <- sum(apply(a_mat, MARGIN = 2, FUN = create_max_and_num_path,
                         N_vec = N_vec, nevents = nevents, k = k,
                         crit_values = crit_value, test_design = test_design, direction = test.direction)) / choose(sum(N_vec), nevents)
  } else if (alternative == "less") {
    test.direction <- ">="
    p_value <- sum(apply(a_mat, MARGIN = 2, FUN = create_max_and_num_path,
                         N_vec = N_vec, nevents = nevents, k = k,
                         crit_values = crit_value, test_design = test_design, direction = test.direction)) / choose(sum(N_vec), nevents)
  } else if (alternative == "two.sided") {
    test.direction <- "<="
    p_value_less <- sum(apply(a_mat, MARGIN = 2, FUN = create_max_and_num_path,
                              N_vec = N_vec, nevents = nevents, k = k,
                              crit_values = crit_value, test_design = test_design, direction = test.direction)) / choose(sum(N_vec), nevents)
    test.direction <- ">="
    p_value_greater <- sum(apply(a_mat, MARGIN = 2, FUN = create_max_and_num_path,
                                 N_vec = N_vec, nevents = nevents, k = k,
                                 crit_values = crit_value, test_design = test_design, direction = test.direction)) / choose(sum(N_vec), nevents)
    p_value <- min(2 * p_value_less, 2 * p_value_greater, 
                   1)
  }
  
  # decide on the null hypothesis
  test_entscheidung <- p_value < level
  
  # Construct Test Decision and Summary
  if (test_entscheidung) {
    test_decision <- "H_0 can be rejected at level "
  } else {
    test_decision <- "H_0 cannot be rejected at level "
  }
  
  summary_message <- paste(
    "Maximally selected Chi^2 Test: \n",
    "Observed Chi^2:", critical_value, "\n",
    "P-value: ", round(p_value, digits = 5), "\n",
    "Test decision: ", test_decision, level, "\n"
  )
  
  if (trace) {
    cat(summary_message, sep = "\n")
  }
  
  # Return Values
  values <- c(p_value, normal_p_value)
  names(values) <- c("Our p_value", "casual p_value")
  
  return(values)
}

# Use case examples:
#example_dataset <- matrix(data = c(9, 6, 11, 4), nrow = 2, byrow = FALSE)
#max.chisq.test(example_dataset, level = 0.05)

#example_dataset <- matrix(data = c(11, 4, 2, 13, 1, 14, 6, 9),nrow = 2, byrow = FALSE)
#result <- max_chisq_test(example_dataset, level = 0.05, alternative = "two.sided")



