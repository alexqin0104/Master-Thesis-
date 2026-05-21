setwd("C:/Users/eloua/Desktop/master aangewezen/extension/typeII")

library(tidyverse)
library(furrr)

source("Analyst Algorithms.R")

plan(multisession)

# -----------------------------
# Normal data generator (AR(1))
# -----------------------------
generate_normal_series <- function(n, phi, sigma_squared) {
  sd <- sqrt(sigma_squared)
  
  if (phi == 0) {
    y <- 5 + rnorm(n, sd = sd)
  } else {
    y <- 5 + as.numeric(
      arima.sim(
        n = n,
        model = list(ar = phi),
        sd = sd
      )
    )
  }
  
  y
}

# -----------------------------
# SETTINGS
# -----------------------------
set.seed(1)

iterations <- 20000
n <- 101
treatment_length <- 15

phi_vals <- c(0, 0.2, 0.4)
sigma_squared_vals <- c(0.25, 1, 2.25)
effect_size <- 0.5

algorithms <- list(
  fixed5 = fixed5#,
 # kaz10 = kaz10,
 # kaz15 = kaz15,
  #vva = vva,
  #gl_full = gl_full,
  #gl_final = gl_final,
  #gl_abs = gl_abs,
  #gl_rel = gl_rel
)

# -----------------------------
# PARAM GRID
# -----------------------------
param_grid <- expand.grid(
  phi = phi_vals,
  sigma_squared = sigma_squared_vals
)

# -----------------------------
# SINGLE PARAM FUNCTION
# -----------------------------
run_one_param <- function(p,
                          phi,
                          sigma_squared,
                          iterations,
                          n,
                          treatment_length,
                          effect_size,
                          algorithms) {
  
  source("Analyst Algorithms.R")
  
  absolute_effect <- effect_size * sqrt(sigma_squared)
  
  out <- vector("list", iterations * length(algorithms))
  k <- 1
  
  for (iter in 1:iterations) {
    
    y_full <- generate_normal_series(n + treatment_length, phi, sigma_squared)
    
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
      treatment_raw <- y_full[(len + 1):(len + treatment_length)]
      treatment <- treatment_raw + absolute_effect
      
      y <- c(baseline, treatment)
      phase <- c(rep(0, len), rep(1, treatment_length))
      
      model <- tryCatch(lm(y ~ phase), error = function(e) NULL)
      
      pval <- if (is.null(model)) NA else
        summary(model)$coefficients["phase", "Pr(>|t|)"]
      
      out[[k]] <- data.frame(
        iter = iter,
        phi = phi,
        sigma_squared = sigma_squared,
        Algorithm = alg_name,
        p_value = pval
      )
      k <- k + 1
    }
  }
  
  bind_rows(out)
}

# -----------------------------
# PARALLEL EXECUTION
# -----------------------------
results <- future_map(
  1:nrow(param_grid),
  .options = furrr_options(seed = TRUE),
  ~ run_one_param(
    p = .x,
    phi = param_grid$phi[.x],
    sigma_squared = param_grid$sigma_squared[.x],
    iterations = iterations,
    n = n,
    treatment_length = treatment_length,
    effect_size = effect_size,
    algorithms = algorithms
  )
) %>% bind_rows()

# -----------------------------
# POWER + TYPE II
# -----------------------------
power_results <- results %>%
  group_by(phi, sigma_squared, Algorithm) %>%
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

save(power_results, file = "TypeIINormal49Observations.Rdata")