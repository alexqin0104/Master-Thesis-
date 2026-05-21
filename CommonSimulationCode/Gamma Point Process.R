rm(list = ls())
setwd("C:/Users/eloua/Desktop/master aangewezen/code uploaden/CommonSimulationCode")
source("Gamma Point Process Simulation Functions New.R")

library(furrr)
set.seed(1)

plan(multisession)

plan()

iterations <- 20000
n <- 101
treatment_length = 15
mu_a <- c(5, 15, 25)
shape <- c(2/5, 2/3, 1, 3/2, 5/2) #just added "1" to the list.
effect_size <- 0.5 #expressed in standard deviations

params <- expand.grid(iterations = iterations, n = n, treatment_length = treatment_length,
                      mu_a = mu_a, shape = shape, effect_size = effect_size)


gamma_point_process_results <- evaluate_by_row(params = params, sim_function = sim_driver)

#two separate datasets are now created. Fixed contains the reference (the one with the baseline length fixed),
#adaptive the data with the algorithms applied 
fixed_results <- gamma_point_process_results %>%
  filter(alg == "fixed5")

adaptive_results <- gamma_point_process_results %>%
  filter(alg != "fixed5")


save(adaptive_results, file = "gppAdaptive.Rdata")
save(fixed_results, file = "gppFixed49.Rdata")

plan(sequential)