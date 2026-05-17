rm(list = ls())
source("Binomial Thinning Simulation.R")
library(future)
library(furrr)
library(tidyverse)
set.seed(1)
plan(multisession)

# ============================================================
# Parameter settings
# ============================================================
iterations       <- 2000
n                <- 40    # max_bl
treatment_length <- 15
mu               <- c(5, 15, 25)
phi              <- c(0, 0.2, 0.4)
effect_size      <- c(0, 1, 2)

params <- expand.grid(
  iterations       = iterations,
  n                = n,
  treatment_length = treatment_length,
  mu               = mu,
  phi              = phi,
  effect_size      = effect_size
)

cat("==================================================\n")
cat("Running multilevel-based MBD simulation\n")
cat(nrow(params), "parameter combinations,",
    iterations, "iterations each\n")
cat("==================================================\n")

# ============================================================
binomial_mbd_results <- evaluate_by_row(
  params       = params,
  sim_function = sim_driver
)

save(binomial_mbd_results, file = "binom_MBD_lme_full2000.Rdata")
write.csv(binomial_mbd_results,
          file = "binom_MBD_lme_full2000.csv",
          row.names = FALSE)

plan(sequential)
cat("Simulation complete. Results saved.\n\n")