#### Script for the real world example

# load packages
library(dplyr)
library(tidyr)
library(binMto)

# source the implemented function
source("./full function.R")
source("./unconditional_exact_test.R")

# insert the information of Charytan et al. (2019)
nobs_itt <- c(51, 27, 26, 25)
nobs_ast <- c(52, 31, 21, 25)

hyperkalemia_itt <- c(9, 4, 4, 8)
hyperkalemia_ast <- c(10, 5, 2, 8)

hypotension_itt <- c(0, 2, 0, 3)
hypotension_ast <- c(0, 2, 0, 3)

# get the four contingency tables
data_hyper_itt <- matrix(c(hyperkalemia_itt, nobs_itt - hyperkalemia_itt),
                         nrow = 2, byrow = TRUE)
data_hyper_ast <- matrix(c(hyperkalemia_ast, nobs_ast - hyperkalemia_ast),
                         nrow = 2, byrow = TRUE)
data_hypo_itt <- matrix(c(hypotension_itt, nobs_itt - hypotension_itt),
                        nrow = 2, byrow = TRUE)
data_hypo_ast <- matrix(c(hypotension_ast, nobs_ast - hypotension_ast),
                        nrow = 2, byrow = TRUE)

# initialize a result dataframe

results <- data.frame(method = character(6),
                      hyperkalemia_itt = numeric(6),
                      hyperkalemia_ast = numeric(6),
                      hypotension_itt = numeric(6),
                      hypotension_ast = numeric(6))

results$method <- list("naive",
                    "bonferroni",
                    "k-group",
                    "flexible",
                    "conditional",
                    "unconditional")

# do the analysis for all dataframes and tests

temp_result <- max_chisq_test(data_hyper_itt, alternative = "two.sided")
results$hyperkalemia_itt <- list(temp_result[2],
                                 temp_result[2] * 3,
                                 chisq.test(data_hyper_itt, correct = FALSE)$p.value,
                                 temp_result[1],
                                 ec.mto(n = nobs_itt, x = hyperkalemia_itt, alternative = "two.sided"),
                                 calculate_exact_uncond_p_value(data_hyper_itt))
# asymptotic gives CI's and no p-value
asymptotic_hyper_itt <- binMto(n = nobs_itt, x = hyperkalemia_itt)$conf.int


temp_result <- max_chisq_test(data_hyper_ast, alternative = "two.sided")
results$hyperkalemia_ast <- list(temp_result[2],
                                 temp_result[2] * 3,
                                 chisq.test(data_hyper_ast, correct = FALSE)$p.value,
                                 temp_result[1],
                                 ec.mto(n = nobs_ast, x = hyperkalemia_ast, alternative = "two.sided"),
                                 calculate_exact_uncond_p_value(data_hyper_ast))
asymptotic_hyper_ast <- binMto(n = nobs_ast, x = hyperkalemia_ast)$conf.int


temp_result <- max_chisq_test(data_hypo_itt, alternative = "two.sided")
results$hypotension_itt <- list(temp_result[2],
                                 temp_result[2] * 3,
                                 chisq.test(data_hypo_itt, correct = FALSE)$p.value,
                                 temp_result[1],
                                 ec.mto(n = nobs_itt, x = hypotension_itt, alternative = "two.sided"),
                                 calculate_exact_uncond_p_value(data_hypo_itt))
asymptotic_hypo_itt <- binMto(n = nobs_itt, x = hypotension_itt)$conf.int

temp_result <- max_chisq_test(data_hypo_ast, alternative = "two.sided")
results$hypotension_ast <- list(temp_result[2],
                                temp_result[2] * 3,
                                chisq.test(data_hypo_ast, correct = FALSE)$p.value,
                                temp_result[1],
                                ec.mto(n = nobs_ast, x = hypotension_ast, alternative = "two.sided"),
                                calculate_exact_uncond_p_value(data_hypo_ast))
asymptotic_hypo_ast <- binMto(n = nobs_ast, x = hypotension_ast)$conf.int

# Show results

results
asymptotic_hyper_itt
asymptotic_hyper_ast
asymptotic_hypo_itt
asymptotic_hypo_ast

## For both Hyperkalemia settings, all methods cant reject the null
## For Hypotension ITT analysis just all methods reject the null,
# except for unconditional exact test
## For Hypotension AST analysis the naive, bonferroni-adjusted and
# k-group approaches reject the null


### Results for new flexible design

design <- matrix(c(1, 2, 0, 0,
                   1, 0, 2, 0,
                   1, 0, 0 ,2,
                   0, 1, 2, 0,
                   0, 1, 0, 2), nrow = 5, byrow = TRUE)

max_chisq_test(data_hyper_itt, designmatrix = design)
max_chisq_test(data_hyper_ast, designmatrix = design)
max_chisq_test(data_hypo_itt, designmatrix = design)
max_chisq_test(data_hypo_ast, designmatrix = design)



