###### Using the ADEMP structure for simulations:

# load libraries
library(checkmate)
library(binMto)

# load implemented functions
source("./full function.R")
source("./unconditional_exact_test.R")

# If trying to reproduce results, then load:
# initial_states <- readRDS("./initial_states.rds")


### Function to generate one contingency table and apply the methods on it ####

# recommended is commenting the conditional out and performing
# the calculations additionally see add_conditional.R

onerep <- function(repetition, nobs = c(100, 100), probs = c(0.5, 0.5)) {
  # check if input is correct
  assert_true(length(nobs) == length(probs) & length(nobs) >= 2)
  assertNumeric(nobs, lower = 1, any.missing = FALSE)
  assert_double(probs, lower = 0, upper = 1)
  
  # initialize
  event_vector <- numeric()
  treatment_vector <- character()
  
  # generate data
  for (i in 1:length(probs)) {
    assign(paste0("trt", i), rbinom(n = nobs[i], prob = probs[i], size = 1))
    event_vector <- c(event_vector, get(paste0("trt", i)))
    treatment_vector <- c(treatment_vector,
                          base::rep(x = paste0("trt", i), times = nobs[i]))
  }
  
  # transform into contingency tables
  df <- data.frame(event = event_vector,
                   treatment = treatment_vector)
  
  contingency_table <- as.matrix(table(df$event, df$treatment))
  if (nrow(contingency_table) < 2) {
    if (rownames(contingency_table) == "1") {
      contingency_table <- rbind(contingency_table, rep(0, length(nobs)))
      } else if (rownames(contingency_table) == "0") {
        contingency_table <- rbind(rep(0, length(nobs)), contingency_table)
        }
    } else if (nrow(contingency_table) == 2) {
      contingency_table <- contingency_table[c(2, 1), ] # flip rows such that events are on the upper side
      }
  
  # Apply methods on contingency table
  # function for flexible exact approach
  test_object <- max.chisq.test(contingency_table, trace = FALSE)
  our_decision <- ifelse(test_object[1] <= 0.05, "reject", "cant reject")
  # naive chi-square test
  basic_decision <- ifelse(test_object[2] <= 0.05, "reject", "cant reject")
  # bonferroni adjustment
  bonferroni_p_value <- test_object[2] * (length(nobs) - 1)
  bonferroni_decision <- ifelse(test_object[2] <= 0.05 / (length(nobs) - 1), "reject", "cant reject")
  # unconditional
  unconditional <- calculate_exact_uncond_p_value(contingency_table) 
  unconditional_decision <- ifelse(unconditional <= 0.05, "reject", "cant reject")
  # asymptotic
  asymptotic_ci <- binMto(x = contingency_table[1, ], n = nobs)
  asymptotic_decision <- ifelse(any((asymptotic_ci$conf.int[, 2] <= 0 &
                                       0 <= asymptotic_ci$conf.int[, 3]) == FALSE), "reject", "cant reject")
  # conditional
  conditional <- ec.mto(n = nobs, x = contingency_table[1, ], alternative = "two.sided")
  conditional_decision <- ifelse(conditional <= 0.05, "reject", "cant reject")
  
  
  output <- list(
    repetition = repetition, 
    probabilities = probs,
    nobs = nobs,
    obs_data = contingency_table,
    our_p_value = test_object[1],
    basic_p_value = test_object[2],
    conditional_p_value = conditional,
    unconditional_p_value = unconditional,
    asymptotic_ci = asymptotic_ci,
    our_decision = our_decision,
    basic_decision = basic_decision,
    bonferroni_decision = bonferroni_decision,
    conditional_decision = conditional_decision,
    unconditional_decision = unconditional_decision,
    asymptotic_decision = asymptotic_decision
  )
  
  return(output)
}


# How does this work with different parameters:
# onerep(nobs = c(100, 100, 100), probs = c(0.5, 0.5, 0.5), repetition = 1)
# this works very well

##### Simulation #####################################################

# prepare to run nsim repetitions
set.seed(111224) # or .Random.seed <- initial_states[1, ]
nsim <- 10000

# prepare a grid to evaluate the differences
list_probs <- seq(from = 0.01, to = 0.51, by = 0.1)
nobs <- c(50, 50, 50, 50)

# create a list to store the information
list_3_no <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (j in 1:length(list_probs)) {
  for (r in 1:nsim) {
    states[(j - 1) * nsim + r, ] <- .Random.seed
    list_3_no[[(j - 1) * nsim + r]] <- onerep(repetition = r, nobs = nobs, probs = rep(list_probs[j], times = length(nobs)))
  }
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_no, file = "list_3_no.rds")
saveRDS(states, file = "states_3_no.rds")



##### Same problem with smaller group size 
# .Random.seed <- initial_states[2, ]

# prepare a grid to evaluate the differences
list_probs <- seq(from = 0.01, to = 0.51, by = 0.1)
nobs <- c(20, 20, 20, 20)

# create a list to store the information
list_3_small <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (j in 1:length(list_probs)) {
  for (r in 1:nsim) {
    states[(j - 1) * nsim + r, ] <- .Random.seed
    list_3_small[[(j - 1) * nsim + r]] <- onerep(repetition = r, nobs = nobs, probs = rep(list_probs[j], times = length(nobs)))
  }
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small, file = "list_3_small.rds")
saveRDS(states, file = "states_3_small.rds")


##### Same problem with smaller sample size 
# .Random.seed <- initial_states[3, ]

# prepare a grid to evaluate the differences
list_probs <- seq(from = 0.01, to = 0.51, by = 0.1)
nobs <- c(10, 10, 10, 10)

# create a list to store the information
list_3_small_equal <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (j in 1:length(list_probs)) {
  for (r in 1:nsim) {
    states[(j - 1) * nsim + r, ] <- .Random.seed
    list_3_small_equal[[(j - 1) * nsim + r]] <- onerep(repetition = r, nobs = nobs, probs = rep(list_probs[j], times = length(nobs)))
  }
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small_equal, file = "list_3_small_equal.rds")
saveRDS(states, file = "states_3_small_equal.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[4, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.2, 0.2, 0.2, 0.15)
nobs <- c(40, 40, 40, 40)

# create a list to store the information
list_3_big_small <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    list_3_big_small[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
  }
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_big_small, file = "list_3_big_small.rds")
saveRDS(states, file = "states_3_big_small.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[5, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.2, 0.2, 0.2, 0.15)
nobs <- c(20, 20, 20, 20)

# create a list to store the information
list_3_medium_small <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_medium_small[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_medium_small, file = "list_3_medium_small.rds")
saveRDS(states, file = "states_3_medium_small.rds")

##### Equal group size, different probs

# prepare a grid to evaluate the differences
list_probs <- c(0.2, 0.2, 0.2, 0.15)
nobs <- c(10, 10, 10, 10)

# create a list to store the information
list_3_small_small <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_small_small[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small_small, file = "list_3_small_small.rds")
saveRDS(states, file = "states_3_medium_small.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[7, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.25, 0.15)
nobs <- c(40, 40, 40, 40)

# create a list to store the information
list_3_big_big <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_big_big[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_big_big, file = "list_3_big_big.rds")
saveRDS(states, file = "states_3_big_big.rds")

##### Equal group size, different probs
# .Random.seed <- initial_states[8, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.25, 0.15)
nobs <- c(20, 20, 20, 20)

# create a list to store the information
list_3_medium_big <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_medium_big[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_medium_big, file = "list_3_medium_big.rds")
saveRDS(states, file = "states_3_medium_big.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[9, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.25, 0.15)
nobs <- c(10, 10, 10, 10)

# create a list to store the information
list_3_small_big <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_small_big[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small_big, file = "list_3_small_big.rds")
saveRDS(states, file = "states_3_small_big.rds")



##### Equal group size, different probs
# .Random.seed <- initial_states[10, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15, 0.35)
nobs <- c(40, 40, 40, 40)

# create a list to store the information
list_3_big_diff <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_big_diff[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_big_diff, file = "list_3_big_diff.rds")
saveRDS(states, file = "states_3_big_diff.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[11, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15, 0.35)
nobs <- c(20, 20, 20, 20)

# create a list to store the information
list_3_medium_diff <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_medium_diff[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_medium_diff, file = "list_3_medium_diff.rds")
saveRDS(states, file = "states_3_medium_diff.rds")


##### Equal group size, different probs
# .Random.seed <- initial_states[12, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15, 0.35)
nobs <- c(10, 10, 10, 10)

# create a list to store the information
list_3_small_diff <- list(list())
states <- matrix(ncol = 626, nrow = nsim)

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_small_diff[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small_diff, file = "list_3_small_diff.rds")
saveRDS(states, file = "states_3_small_diff.rds")




##### Different group sizes, different probabilities
# .Random.seed <- initial_states[13, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15, 0.35)
nobs <- c(40, 20, 20, 20)

# create a list to store the information
list_3_big_size <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    list_3_big_size[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
  }
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_big_size, file = "list_3_big_size.rds")
saveRDS(states, file = "states_3_big_size.rds")


##### Different group sizes, different probabilities
# .Random.seed <- initial_states[14, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15, 0.35)
nobs <- c(20, 10, 10, 10)

# create a list to store the information
list_3_small_size <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_3_small_size[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_3_small_size, file = "list_3_small_size.rds")
saveRDS(states, file = "states_3_small_size.rds")



##### Equal group sizes, different probabilities, k = 3
# .Random.seed <- initial_states[15, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.15)
nobs <- c(20, 20, 20)

# create a list to store the information
list_2_small_big <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_2_small_big[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_2_small_big, file = "list_2_small_big.rds")
saveRDS(states, file = "states_2_small_big.rds")

##### Equal group sizes, different probabilities, k = 3
# .Random.seed <- initial_states[16, ]

# prepare a grid to evaluate the differences
list_probs <- c(0.3, 0.3, 0.25)
nobs <- c(20, 20, 20)

# create a list to store the information
list_2_small_small <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_2_small_small[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_2_small_small, file = "list_2_small_small.rds")
saveRDS(states, file = "states_2_small_small.rds")


##### Equal group sizes, different probabilities, k = 6
# .Random.seed <- initial_states[17, ]

## small difference
list_probs <- c(0.2, 0.2, 0.15, 0.25, 0.2, 0.18)
nobs <- c(20, 20, 20, 20, 20, 20)

# create a list to store the information
list_5_small_small <- list(list())
states <- matrix(ncol = 626, nrow = (nsim*length(list_probs)))

# run all nsim reps for all the information
start_time <- proc.time()
for (r in 1:nsim) {
  states[r, ] <- .Random.seed
  list_5_small_small[[r]] <- onerep(repetition = r, nobs = nobs, probs = list_probs)
}
end_time <- proc.time()
end_time - start_time

# save data from analysis
saveRDS(list_5_small_small, file = "list_5_small_small.rds")
saveRDS(states, file = "states_5_small_small.rds")



######### Problem with a repetition? Do:
#.Random.seed <- states[25, ]
#onerep(rep = 25, nobs = nobs, probs = c(0.01, 0.01, 0.01))

