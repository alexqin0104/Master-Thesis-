setwd("C:/Users/eloua/Desktop/master aangewezen/extension/typeII")

library(tidyverse)
library(ARPobservation)
library(furrr)

source("Analyst Algorithms.R")

plan(multisession)

# -----------------------------
# GPP generator
# -----------------------------
generate_gpp_series <- function(n, mu_a, shape) {
  event_counting(
    r_behavior_stream(
      n = n,
      mu = 0,
      lambda = 1 / mu_a,
      F_event = F_exp(),
      F_interim = F_gam(shape),
      stream_length = 1
    )
  )
}

generate_treatment_gpp <- function(length, mu_a, shape, delta) {
  generate_gpp_series(
    n = length,
    mu_a = mu_a * exp(delta),
    shape = shape
  )
}

# -----------------------------
# SETTINGS
# -----------------------------
set.seed(1)

iterations <- 20000
n <- 101
treatment_length <- 15

mu_a_vals <- c(5, 15, 25)
shape_vals <- c(2/5, 2/3, 1, 3/2, 5/2)

delta <- log(1.5)

algorithms <- list(
  fixed5 = fixed5#,
  #kaz10 = kaz10,
  #kaz15 = kaz15,
  #vva = vva,
  #gl_full = gl_full,
  #gl_final = gl_final,
  #gl_abs = gl_abs,
  #gl_rel = gl_rel
)

param_grid <- expand.grid(
  mu_a = mu_a_vals,
  shape = shape_vals
)

# -----------------------------
# SINGLE PARAM FUNCTION
# -----------------------------
run_one_param <- function(p, mu_a, shape,
                          iterations, n, treatment_length,
                          delta, algorithms) {
  source("Analyst Algorithms.R")
  
  results_local <- vector("list", iterations * length(algorithms))
  k <- 1
  
  for (iter in 1:iterations) {
    
    y_full <- generate_gpp_series(n + treatment_length, mu_a, shape)
    
    for (alg_name in names(algorithms)) {
      
      alg <- algorithms[[alg_name]]
      
      len <- NA
      for (i in 5:n) {
        if (alg(y_full[1:i])) {
          len <- i
          break
        }
      }
      if (is.na(len)) len <- n
      
      baseline <- y_full[1:len]
      
      treatment <- generate_treatment_gpp(
        length = treatment_length,
        mu_a = mu_a,
        shape = shape,
        delta = delta
      )
      
      y <- c(baseline, treatment)
      phase <- c(rep(0, len), rep(1, treatment_length))
      
      model <- tryCatch(
        glm(y ~ phase, family = quasipoisson(link = "log")),
        error = function(e) NULL
      )
      
      pval <- if (is.null(model)) NA else
        summary(model)$coefficients["phase", "Pr(>|t|)"]
      
      results_local[[k]] <- data.frame(
        iter = iter,
        mu_a = mu_a,
        shape = shape,
        Algorithm = alg_name,
        p_value = pval
      )
      
      k <- k + 1
    }
  }
  
  bind_rows(results_local)
}

# -----------------------------
# PARALLEL EXECUTION
# -----------------------------
results <- future_map(
  1:nrow(param_grid),
  .options = furrr_options(seed = TRUE),
  ~ run_one_param(
    p = .x,
    mu_a = param_grid$mu_a[.x],
    shape = param_grid$shape[.x],
    iterations = iterations,
    n = n,
    treatment_length = treatment_length,
    delta = delta,
    algorithms = algorithms
  )
) %>% bind_rows()

# -----------------------------
# FINAL SUMMARY
# -----------------------------
power_results <- results %>%
  group_by(mu_a, shape, Algorithm) %>%
  summarise(
    R_eff = sum(!is.na(p_value)),
    
    power = mean(p_value < 0.05, na.rm = TRUE),
    beta  = mean(p_value >= 0.05, na.rm = TRUE),
    
    MCSE_power = sqrt(power * (1 - power) / R_eff),
    MCSE_beta  = sqrt(beta * (1 - beta) / R_eff),
    
    power_lower = pmax(0, power - 1.96 * MCSE_power),
    power_upper = pmin(1, power + 1.96 * MCSE_power),
    
    beta_lower = pmax(0, beta - 1.96 * MCSE_beta),
    beta_upper = pmin(1, beta + 1.96 * MCSE_beta),
    
    .groups = "drop"
  )

print(power_results)

save(power_results, file = "TypeIIGamma49Observations.Rdata")