library(tidyverse)
library(furrr)

source("Analyst Algorithms.R")

plan(multisession)

# -----------------------------
# AR(1) generator
# -----------------------------
r_pois_AR1 <- function(n, mu, phi) {
  y <- numeric(n)
  y[1] <- rpois(1, mu)
  
  for (t in 2:n) {
    y[t] <- rbinom(1, size = y[t-1], prob = phi) +
      rpois(1, lambda = (1 - phi) * mu)
  }
  
  y
}

# -----------------------------
# Treatment generator
# -----------------------------
generate_treatment_AR1 <- function(last_y, mu, phi, delta, length) {
  mu_treat <- mu * exp(delta)
  
  y <- numeric(length)
  y[1] <- rbinom(1, last_y, phi) +
    rpois(1, (1 - phi) * mu_treat)
  
  for (t in 2:length) {
    y[t] <- rbinom(1, y[t-1], phi) +
      rpois(1, (1 - phi) * mu_treat)
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

mu_vals  <- c(5, 15, 25)
phi_vals <- c(0, 0.2, 0.4)

effect_size <- 0.5
delta <- log(1.5)

algorithms <- list(
  fixed5 = fixed5#,
 # kaz10 = kaz10,
  #kaz15 = kaz15,
  #vva = vva,
  #gl_full = gl_full,
  #gl_final = gl_final,
  #gl_abs = gl_abs,
  #gl_rel = gl_rel
)

param_grid <- expand.grid(
  mu = mu_vals,
  phi = phi_vals
)

# -----------------------------
# SINGLE PARAM RUN
# -----------------------------
run_one_param <- function(p, mu, phi, iterations, n, treatment_length,
                          effect_size, delta, algorithms) {
  #source("Analyst Algorithms.R")
  
  absolute_effect <- effect_size * sqrt(mu)
  
  out <- vector("list", iterations * length(algorithms))
  k <- 1
  
  for (iter in 1:iterations) {
    
    y_full <- r_pois_AR1(n + treatment_length, mu, phi)
    
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
        mu = mu,
        phi = phi,
        Algorithm = alg_name,
        len = len,
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
    mu = param_grid$mu[.x],
    phi = param_grid$phi[.x],
    iterations = iterations,
    n = n,
    treatment_length = treatment_length,
    effect_size = effect_size,
    delta = delta,
    algorithms = algorithms
  )
) %>% bind_rows()

# -----------------------------
# POWER + TYPE II
# -----------------------------
power_results <- results %>%
  group_by(mu, phi, Algorithm) %>%
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

save(power_results, file = "TypeIIBinomial49Observations.Rdata")