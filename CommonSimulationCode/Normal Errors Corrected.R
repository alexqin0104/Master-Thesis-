rm(list = ls())
setwd("C:/Users/eloua/Desktop/master aangewezen/code uploaden/CommonSimulationCode")
source("Normal Errors Simulation Functions Corrected New.R")

library(furrr)
plan(multisession)

iterations <- 20000
n <- 101
treatment_length = 15
phi <- c(0, 0.20, 0.40)
sigma_squared <- c(0.25, 1, 2.25)
effect_size <- 0.5

params <- expand.grid(
  iterations = iterations,
  n = n,
  treatment_length = treatment_length,
  phi = phi,
  sigma_squared = sigma_squared,
  effect_size = effect_size
)

set.seed(1)

normal_errors_results <-
  evaluate_by_row(params = params, sim_function = sim_driver)



#two separate datasets are now created. Fixed contains the reference (the one with the baseline length fixed),
#adaptive the data with the algorithms applied 
fixed_results <- normal_errors_results %>%
  filter(Algorithm == "fixed5")

adaptive_results <- normal_errors_results %>%
  filter(Algorithm != "fixed5")


save(adaptive_results, file = "normCorrectedAdaptive.Rdata")
save(fixed_results, file = "normCorrectedFixedTEST.Rdata")
