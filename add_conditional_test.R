## Add the conditional exact test to the data and pseudo p_values for 
# asymptotic approach

# load packages
library(tidyverse)
library(rsimsum)
library(binMto)

# load the first list: equal proportions, big groups
list_3_no <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_no.rds")

# add the calculation of:
# and a dummy asymptotic p-value which is needed for further calculations and
# unconditional exact test
list_3_no <- lapply(list_3_no, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

# add the calculation of the conditional test from binMto 
# to a quarter of the data for faster computation
# split the information of the "observed data" into a smaller list
# there are 10.000 repetitions per type of setting
conditional_list_3_no <- lapply(list_3_no[c(1:2500,
                                            10001:12500,
                                            20001:22500,
                                            30001:32500,
                                            40001:42500,
                                            50001:52500)],
                                function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_no <- lapply(conditional_list_3_no, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_no, file = "./Simulation_results/list_3_no.rds")
saveRDS(conditional_list_3_no, file ="./Simulation_results/conditional_list_3_no.rds")

rm(list_3_no)
rm(conditional_list_3_no)

# load the second list: equal proportions, smaller sample size

list_3_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small.rds")

list_3_small <- lapply(list_3_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_small <- lapply(list_3_small[c(1:2500,
                                                  10001:12500,
                                                  20001:22500,
                                                  30001:32500,
                                                  40001:42500,
                                                  50001:52500)],
                                   function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small <- lapply(conditional_list_3_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small, file = "./Simulation_results/list_3_small.rds")
saveRDS(conditional_list_3_small, file ="./Simulation_results/conditional_list_3_small.rds")

rm(list_3_small)
rm(conditional_list_3_small)


# load the third list: equal proportions, much smaller sample size

list_3_small_equal <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small_equal.rds")

list_3_small_equal <- lapply(list_3_small_equal, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_small_equal <- lapply(list_3_small_equal[c(1:2500,
                                                              10001:12500,
                                                              20001:22500,
                                                              30001:32500,
                                                              40001:42500,
                                                              50001:52500)],
                                         function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small_equal <- lapply(conditional_list_3_small_equal, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small_equal, file = "./Simulation_results/list_3_small_equal.rds")
saveRDS(conditional_list_3_small_equal, file ="./Simulation_results/conditional_list_3_small_equal.rds")

rm(list_3_small_equal)
rm(conditional_list_3_small_equal)

# load the fourth list: small difference, big sample size

list_3_big_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_big_small.rds")

list_3_big_small <- lapply(list_3_big_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_big_small <- lapply(list_3_big_small[1:2500],
                                       function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_big_small <- lapply(conditional_list_3_big_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_big_small, file = "./Simulation_results/list_3_big_small.rds")
saveRDS(conditional_list_3_big_small, file ="./Simulation_results/conditional_list_3_big_small.rds")

rm(list_3_big_small)
rm(conditional_list_3_big_small)


# load the fifth list: small difference, medium sample size

list_3_medium_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_medium_small.rds")

list_3_medium_small <- lapply(list_3_medium_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_medium_small <- lapply(list_3_medium_small[1:5000],
                                          function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_medium_small <- lapply(conditional_list_3_medium_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_medium_small, file = "./Simulation_results/list_3_medium_small.rds")
saveRDS(conditional_list_3_medium_small, file ="./Simulation_results/conditional_list_3_medium_small.rds")

rm(list_3_medium_small)
rm(conditional_list_3_medium_small)


# load the sixth list: small difference, small sample size

list_3_small_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small_small.rds")

list_3_small_small <- lapply(list_3_small_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_small_small <- lapply(list_3_small_small,
                                         function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small_small <- lapply(conditional_list_3_small_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small_small, file = "./Simulation_results/list_3_small_small.rds")
saveRDS(conditional_list_3_small_small, file ="./Simulation_results/conditional_list_3_small_small.rds")

rm(list_3_small_small)
rm(conditional_list_3_small_small)


# load the seventh list: bigger difference, big sample size

list_3_big_big <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_big_big.rds")

list_3_big_big <- lapply(list_3_big_big, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_big_big <- lapply(list_3_big_big[1:2500],
                                     function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_big_big <- lapply(conditional_list_3_big_big, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_big_big, file = "./Simulation_results/list_3_big_big.rds")
saveRDS(conditional_list_3_big_big, file ="./Simulation_results/conditional_list_3_big_big.rds")

rm(list_3_big_big)
rm(conditional_list_3_big_big)


# load the eight list: bigger difference, medium sample size

list_3_medium_big <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_medium_big.rds")

list_3_medium_big <- lapply(list_3_medium_big, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_medium_big <- lapply(list_3_medium_big[1:5000],
                                        function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_medium_big <- lapply(conditional_list_3_medium_big, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_medium_big, file = "./Simulation_results/list_3_medium_big.rds")
saveRDS(conditional_list_3_medium_big, file ="./Simulation_results/conditional_list_3_medium_big.rds")

rm(list_3_medium_big)
rm(conditional_list_3_medium_big)


# load the ninth list: bigger difference, small sample size

list_3_small_big <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small_big.rds")

list_3_small_big <- lapply(list_3_small_big, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_small_big <- lapply(list_3_small_big,
                                       function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small_big <- lapply(conditional_list_3_small_big, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small_big, file = "./Simulation_results/list_3_small_big.rds")
saveRDS(conditional_list_3_small_big, file ="./Simulation_results/conditional_list_3_small_big.rds")

rm(list_3_small_big)
rm(conditional_list_3_small_big)


# load the tenth list: difference in different directions, big sample size

list_3_big_diff <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_big_diff.rds")

list_3_big_diff <- lapply(list_3_big_diff, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_big_diff <- lapply(list_3_big_diff[1:2500],
                                      function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_big_diff <- lapply(conditional_list_3_big_diff, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_big_diff, file = "./Simulation_results/list_3_big_diff.rds")
saveRDS(conditional_list_3_big_diff, file ="./Simulation_results/conditional_list_3_big_diff.rds")

rm(list_3_big_diff)
rm(conditional_list_3_big_diff)


# load the elevnth list: difference in different directions, medium sample size

list_3_medium_diff <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_medium_diff.rds")

list_3_medium_diff <- lapply(list_3_medium_diff, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_medium_diff <- lapply(list_3_medium_diff[1:5000],
                                         function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_medium_diff <- lapply(conditional_list_3_medium_diff, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_medium_diff, file = "./Simulation_results/list_3_medium_diff.rds")
saveRDS(conditional_list_3_medium_diff, file ="./Simulation_results/conditional_list_3_medium_diff.rds")

rm(list_3_medium_diff)
rm(conditional_list_3_medium_diff)


# load the twelvth list: difference in different directions, small sample size

list_3_small_diff <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small_diff.rds")

list_3_small_diff <- lapply(list_3_small_diff, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_small_diff <- lapply(list_3_small_diff,
                                        function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small_diff <- lapply(conditional_list_3_small_diff, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small_diff, file = "./Simulation_results/list_3_small_diff.rds")
saveRDS(conditional_list_3_small_diff, file ="./Simulation_results/conditional_list_3_small_diff.rds")

rm(list_3_small_diff)
rm(conditional_list_3_small_diff)


# load the thirteenth list: different group sizes, big sample size

list_3_big_size <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_big_size.rds")

list_3_big_size <- lapply(list_3_big_size, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_3_big_size <- lapply(list_3_big_size[1:5000],
                                      function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_big_size <- lapply(conditional_list_3_big_size, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_big_size, file = "./Simulation_results/list_3_big_size.rds")
saveRDS(conditional_list_3_big_size, file ="./Simulation_results/conditional_list_3_big_size.rds")

rm(list_3_big_size)
rm(conditional_list_3_big_size)


# load the fourteenth list: different group sizes, small sample size
list_3_small_size <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_3_small_size.rds")

list_3_small_size <- lapply(list_3_small_size, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
})

conditional_list_3_small_size <- lapply(list_3_small_size,
                                        function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_3_small_size <- lapply(conditional_list_3_small_size, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_3_small_size, file = "./Simulation_results/list_3_small_size.rds")
saveRDS(conditional_list_3_small_size, file ="./Simulation_results/conditional_list_3_small_size.rds")

rm(list_3_small_size)
rm(conditional_list_3_small_size)


# load the sixteenth list: 2 groups, big difference, equal sample size
list_2_small_big <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_2_small_big.rds")

list_2_small_big <- lapply(list_2_small_big, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_2_small_big <- lapply(list_2_small_big,
                                       function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_2_small_big <- lapply(conditional_list_2_small_big, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_2_small_big, file = "./Simulation_results/list_2_small_big.rds")
saveRDS(conditional_list_2_small_big, file ="./Simulation_results/conditional_list_2_small_big.rds")

rm(list_2_small_big)
rm(conditional_list_2_small_big)


# load the seventeenth list: 2 groups, big difference, equal sample size
list_2_small_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_2_small_small.rds")

list_2_small_small <- lapply(list_2_small_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_2_small_small <- lapply(list_2_small_small,
                                         function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_2_small_small <- lapply(conditional_list_2_small_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_2_small_small, file = "./Simulation_results/list_2_small_small.rds")
saveRDS(conditional_list_2_small_small, file ="./Simulation_results/conditional_list_2_small_small.rds")

rm(list_2_small_small)
rm(conditional_list_2_small_small)


# load the seventeenth list: 5 groups, small difference, equal small sample size
list_5_small_small <- readRDS("/nfsmb/koll/marie.kraft/Masters_thesis/maximally_selected_statistics/Simulation_results/list_5_small_small.rds")

list_5_small_small <- lapply(list_5_small_small, function(entry) {
  entry$asymptotic_p_value <- ifelse(is.null(entry$asymptotic_p_value), NA, entry$asymptotic_p_value)
  return(entry)
})

conditional_list_5_small_small <- lapply(list_5_small_small[1:100],
                                         function(x) x[c("repetition", "probabilities", "nobs", "obs_data")])
# run the conditional function on all those datasets
conditional_list_5_small_small <- lapply(conditional_list_5_small_small, function(entry) {
  entry$conditional_p_value <- ec.mto(n = entry$nobs, x = entry$obs_data[1, ], alternative = "two.sided")
  entry$conditional_decision <- ifelse(entry$conditional_p_value <= 0.05, "reject", "cant reject")
  
  return(entry)
})

saveRDS(list_5_small_small, file = "./Simulation_results/list_5_small_small.rds")
saveRDS(conditional_list_5_small_small, file ="./Simulation_results/conditional_list_5_small_small.rds")

rm(list_5_small_small)
rm(conditional_list_5_small_small)
