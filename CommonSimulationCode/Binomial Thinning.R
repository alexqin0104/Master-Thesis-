rm(list = ls())
setwd("C:/Users/eloua/Desktop/master aangewezen/code uploaden/CommonSimulationCode")
#setwd("C:/Users/eloua/Desktop/master aangewezen/extension/extension 1 treatment effect/")
source("Binomial Thinning Simulation Functions New.R")
set.seed(1)
plan(sequential)


iterations <- 20000

#n is the maximum length of the baseline
n <- 101
treatment_length = 15
mu <- c(5, 15, 25)
phi <- c(0, 0.2, 0.4)
effect_size <- 0.5 #expressed in standard deviations


params <- expand.grid(iterations = iterations, n = n,treatment_length = treatment_length,
                      mu = mu, phi = phi, effect_size = effect_size)

binomial_thinning_results <- evaluate_by_row(params = params, sim_function = sim_driver)

#splitting the fixed results from the adaptive results
#Commenting fixed49 out, since it's only used for the Power analysis.
fixed_results <- binomial_thinning_results %>%
  filter(Algorithm == "fixed5")

adaptive_results <- binomial_thinning_results %>%
  filter(Algorithm != "fixed5")


save(fixed_results, file = "binomthinFixed.Rdata")
save(adaptive_results, file = "binomthinAdaptive.Rdata")