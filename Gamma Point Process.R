rm(list = ls())
source("Gamma Point Process Simulation.R")
library(future)
library(furrr)
library(tidyverse)
set.seed(1)
plan(multisession)

# ============================================================
# Gamma Point Process MBD simulation
# ============================================================
iterations       <- 2000
n                <- 40    # max_bl
treatment_length <- 15
mu               <- c(5, 15, 25)
shape            <- c(2/5, 2/3, 1, 3/2, 5/2)
effect_size      <- c(0, 1, 2)

params <- expand.grid(
  iterations       = iterations,
  n                = n,
  treatment_length = treatment_length,
  mu               = mu,
  shape            = shape,
  effect_size      = effect_size
)

cat("==================================================\n")
cat("Running multilevel-based MBD simulation (gamma)\n")
cat(nrow(params), "parameter combinations,",
    iterations, "iterations each\n")
cat("==================================================\n")

# ============================================================
gamma_mbd_results <- evaluate_by_row(
  params       = params,
  sim_function = sim_driver
)

save(gamma_mbd_results, file = "gamma_MBD_lme_full2000.Rdata")
write.csv(gamma_mbd_results,
          file = "gamma_MBD_lme_full2000.csv",
          row.names = FALSE)

plan(sequential)
cat("Simulation complete. Results saved.\n\n")