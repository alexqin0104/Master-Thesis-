rm(list = ls())
source("Normal Errors Simulation.R")
library(future)
library(furrr)
library(tidyverse)
set.seed(1)
plan(multisession)

# ============================================================
# Normal Errors MBD simulation
# ============================================================
iterations       <- 2000
n                <- 40    # max_bl
treatment_length <- 15
marginal_var     <- c(0.25, 1, 2.25)
phi              <- c(0, 0.20, 0.40)
effect_size      <- c(0, 1, 2)

params <- expand.grid(
  iterations       = iterations,
  n                = n,
  treatment_length = treatment_length,
  marginal_var     = marginal_var,
  phi              = phi,
  effect_size      = effect_size
)

cat("==================================================\n")
cat("Running multilevel-based MBD simulation (normal)\n")
cat(nrow(params), "parameter combinations,",
    iterations, "iterations each\n")
cat("==================================================\n")

# ============================================================
normal_mbd_results <- evaluate_by_row(
  params       = params,
  sim_function = sim_driver
)

save(normal_mbd_results, file = "normal_MBD_lme_full2000.Rdata")
write.csv(normal_mbd_results,
          file = "normal_MBD_lme_full2000.csv",
          row.names = FALSE)

plan(sequential)
cat("Simulation complete. Results saved.\n\n")
