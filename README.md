# Maximally selected chi-square statistics for K-group comparison: a flexible exact approach
Repository for my Thesis about maximally selected chi-square statistics.

## Abstract
This thesis introduces a flexible exact approach for modeling the distribution of
the maximally selected chi-square statistic (χ2_max) in K-group comparisons, effectively
addressing challenges related to small sample sizes and multiple dependent
tests. The method’s adaptability allows for dynamic customization of the alternative
hypothesis by incorporating or omitting tests within the standard control-versus-
treatment framework.
An extensive simulation study compares this approach with several established
statistical methods to evaluate their ability to control the type-I error rate while
maximizing statistical power. The methods examined include the naive χ2-test,
the Bonferroni-adjusted χ2-test, the K-group χ2-test, the newly introduced flexible
exact χ2
max-test, the conditional exact test, the unconditional exact test, and the
asymptotic simultaneous confidence interval method.
To demonstrate its practical relevance, the proposed method is applied to a real-
world multi-dose trial. Simulation results show that conditional methods based on
the flexible exact framework outperform other approaches in accurately analyzing
studies with multiple treatment arms and binary outcomes, making them particu-
larly suitable for multi-dose and dose-finding trials with small sample sizes.


## How to work with these files:
full function.R provides the function max.chisq.test to get the p_value for the flexible
conditional χ2_max
uncontiional_exact_test is the implemented Version of the unconditional exact test for
a test statistic with unpooled variance estimation by Koch and Hothorn (1999)

To recreate the simulation, run the runsim.R document, but leave the conditional exact test out
Afterwards run add_conditional_test.R
Alternatively the simulation results are safed in Simulation_results/
To get the results and plots run results.R and get_results.R

To recreate the results of the real-world example, run real-world.R


Koch, H.-F. and Hothorn, L. A. (1999). Exact unconditional distributions for dichoto-
mous data in many-to-one comparisons, Journal of Statistical Planning and Inference
82(1): 83–99.
URL: https://www.sciencedirect.com/science/article/pii/S0378375899000336
