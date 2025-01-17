## Code to produce the p value for the exact unconditional test

#' Exact Unconditional P-Value Calculation
#'
#' Computes the exact unconditional p-value for comparing a control group to multiple treatment groups  
#' using a test statistic based on unpooled variance estimation. This method accounts for small sample sizes  
#' and avoids reliance on asymptotic approximations.  
#'
#' @param data A 2xK numeric matrix with event counts for the control (first column) and treatment groups.  
#' @return A numeric value representing the exact unconditional p-value.  
#' @details The function calculates the maximum standardized difference in event rates across groups,  
#' evaluates cumulative binomial probabilities over a grid of nuisance parameters, and determines the  
#' tail probability for hypothesis testing.  
#' @examples  
#' example_dataset <- matrix(c(11, 4, 2, 13, 1, 14, 6, 9), nrow = 2)  
#' calculate_exact_uncond_p_value(example_dataset)  
#' @export  


calculate_exact_uncond_p_value <- function(data) {
  
  N_vec <- apply(data, MARGIN = 2, FUN = sum) # Vector with group sizes
  n0 <- N_vec[1] # Sample size for the control group
  n1 <- N_vec[2] # Sample size for each treatment group when those are equally big
  k <- length(N_vec) - 1 # Number of treatment groups
  nevents <- data[1, ] # Vector of events
  
  teststat <- function(N_vec, nevents) {
    nc <- N_vec[1] # ich weiß, dass ich mich hier doppel, sorry. Code ist tzd schnell
    pc <- nevents[1]/nc # observed probability for events \hat{\pi_0}
    ntreat <- N_vec[-1] 
    ptreat <- nevents[-1]/ntreat # observed probability for the groups \hat{\pi_i}
    tstat <- numeric(length = length(ptreat)) # store Test statistics in vector
    for (e in 1:length(ptreat)) {
      # Statistic with unpooled variance estimation as in 2. Methods formula 2
      # the + 1e-04 is to make sure that we don't divide by 0
      tstat[e] <- (ptreat[e] - pc)/sqrt(pc * (1 - pc)/nc + 
                                          ptreat[e] * (1 - ptreat[e])/ntreat[e] + 1e-04)
    }
    return(max(tstat)) # return the T^{UP}
  }
  
  t0 <- teststat(N_vec, nevents) # Observed value of the max statistic 

  # Grid for pi values (nuisance parameter)
  # As described in Section 3.1
  pi_values <- seq(0, 0.99, by = 0.01)

  # Function to compute cumulative binomial probability
  # For the first part of Equation 2.3
  binom_cdf <- function(n, p, x) {
    return(pbinom(x, size = n, prob = p))
  }
  
  # Getting the limit index for each i=0,...,n0
  # Initialize the limit index data, is_smaller_as_t0 is for checking in debugging
  limit_index <- data.frame(index_number = numeric(),
                            limit = numeric(),
                            is_smaller_as_t0 = logical())
  
  # For each i=0,...,n0
  for (i in 0:n0) {
    j <- 0
    is_smaller <- TRUE
    
    # Initialize the latest limit index
    latest <- data.frame(index_number = 0,
                         limit = 0,
                         is_smaller_as_t0 = TRUE)
    
    # Loop as long as S(j,i)<t0
    while (is_smaller & (j <= n1)) {
    p_i <- i/n0
    p_j <- j/n1
    
    # Compute S(j,i)
    statistic <- (p_j - p_i)/sqrt((p_i * (1-p_i)) / n0 +
                                    (p_j * (1-p_j)) /n1 + 1e-04)
    
    # Store information
    latest[1, ] <- list(i, j, is_smaller)
     
    # Prepare for the next step, check if still S(j,i)<t0
    is_smaller <- statistic < t0
    j <- j + 1
    }
    
    # For each i=0,...,n0 there is one limit
    limit_index <- rbind(limit_index, latest)
  }
  
  # Define L as minimum index with l(j,t0)=n1 for each j >= L
  L <- min(limit_index[which(limit_index$limit == n1), "index_number"], n0 + 1)

  
  # Function to calculate tail probability for each pi as Equation 2.3
  calculate_tail_probability <- function(pi) {
    #Check for n0 == 0 to avoid NaNs
    if (n0 == 0) {
      warning("n0 is 0; setting L to 0 and adjusting calculations.")
      L <- 0 
    }
    
    # Compute cumulative probability B(n0, pi, L - 1)
    B_n0_pi_L_minus_1 <- binom_cdf(n0, pi, L - 1)
  
    # Initialize sum_product
      sum_product <- 0
      # Outer Sum
    for (i0 in 0:L) {
      b_n0_pi_i0 <- dbinom(x = i0, size = n0, prob = pi)
      inner_sum <- 0
      # Inner Sum
      if (i0 == n0 + 1) {
        inner_sum <- inner_sum + 0
      } else {
      for (i1 in 0:limit_index[which(limit_index$index_number == i0), "limit"]) {
        inner_sum <- inner_sum + dbinom(x = i1, size = n1, prob = pi)
      }
      }
      sum_product <- sum_product + b_n0_pi_i0 * (inner_sum ^ k)
    }
    
  
    # Return tail probability
    tail_probability <- B_n0_pi_L_minus_1 - sum_product
  
    return(max(0, tail_probability))
  }

  # Compute the tail probability for each pi
  tail_probs <- sapply(pi_values, calculate_tail_probability)

  # The p-value is the maximum of the tail probabilities
  p_value <- max(tail_probs)
  return(p_value)

}

# use case example
#example_dataset <- matrix(data = c(9, 6, 11, 4), nrow = 2, byrow = FALSE)
#calculate_exact_uncond_p_value(example_dataset)
