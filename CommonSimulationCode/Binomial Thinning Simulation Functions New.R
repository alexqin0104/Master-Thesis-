source("Analyst Algorithms.R")
source("common sim functions.R")

r_pois_AR1 <- function(replicate_n, mu, phi) {
  
  len <- length(mu)
  ar_vec <- ifelse(mu[-1] - phi * mu[-len] >= 0, phi, mu[-1] / mu[-len])
  lambda <- pmax(0, mu[-1] - ar_vec * mu[-len])
  
  y <- vector(mode = "numeric", length = len)
  y[1] <- rpois(1, lambda = mu[1])
  for (i in 1:(len-1)) y[i+1] <- rbinom(1, size = y[i], prob = ar_vec[i]) + rpois(1, lambda = lambda[i])
  data.frame(replicate = replicate_n, outcome = y)
}



generate_treatment_AR1 <- function(last_y, mu_base, phi, delta, length) {
  
  mu_treat <- mu_base * exp(delta)
  
  y <- numeric(length)
  y[1] <- rbinom(1, size = last_y, prob = phi) + 
    rpois(1, lambda = (1 - phi) * mu_treat)
  
  for (i in 2:length) {
    y[i] <- rbinom(1, size = y[i-1], prob = phi) + 
      rpois(1, lambda = (1 - phi) * mu_treat)
  }
  
  y
}




sim_driver <- function(iterations, n, treatment_length, mu, phi, effect_size, seed = NULL){
  
  set.seed(1)
  
  alg_vec <- c(
    "fixed5" = fixed5, #this is the "algorithm" used to generate the reference data
#    "fixed49" = fixed49, #this one is commented out here, since it's only used for the Power analysis.
    "kaz10" = kaz10,
    "kaz15" = kaz15,
    "vva" = vva,
    "gl_full" = gl_full,
    "gl_final" = gl_final,
    "gl_abs" = gl_abs,
    "gl_rel" = gl_rel
  )
  
  total_length = n + treatment_length
  # changed the total number of data points generated
  
  dat <- map_dfr(1:iterations, r_pois_AR1,
                 mu = rep(mu, total_length), phi = phi) %>%
    mutate(mu = mu, phi = phi)  #TYPE II
  
  lengths <- dat %>%
    group_by(replicate) %>%
    group_modify(~ map_dfr(.x = alg_vec,
                           .f = min_length,
                           outcomes = .x$outcome,
                           phase_length = n,
                           begin_at = 3)) %>%
    arrange(replicate) %>%
    gather(key = "Algorithm", value = "length", fixed5:gl_rel)
  
  
  # Adding the effect to the treatment data set
  # We do this by simply adding the absolute_effect to the treatment_length last observations
  absolute_effect = effect_size * sqrt(mu)

  
  
  
  # Adding the effect to the treatment data.

  # outcome = original generated data
  # outcome_effect = treatment data with added effect

  dat_effect <- lengths %>%
    full_join(dat) %>%
    group_by(replicate, Algorithm) %>%
    mutate(
      len = first(length),
      outcome_effect = ifelse(
        row_number() >= (len + 1) & row_number() <= (len + treatment_length),
        outcome + absolute_effect,
        outcome
      )
    ) %>%
    ungroup()
  
  
  # Data with phase indicator for type I error calculation
  dat_reg <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, Algorithm) %>%
    mutate(
      len = first(length),
      time = row_number(),
      phase = ifelse(time <= len, 0, 1)
    ) %>%
    ungroup()
  
  
  
  #Need to be very careful to consider which data is being used for the calculations.
  
  
  # Baseline stats
  baseline_stats <- lengths %>%
    full_join(dat) %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length),
      base_mean = ifelse(len < 101, mean(outcome[1:len]), NA),
      base_var  = ifelse(len < 101, var(outcome[1:len]), NA)
    )
  
  # Treatment stats, without effect added to the treatment data
  treatment_stats <- lengths %>%
    full_join(dat) %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length), #this is the baseline length, still needed for the selection of the length
      treat_mean =  mean(outcome[(len+1):(len + treatment_length)]), #we do not need to check any condition anymore for the treatment, since we are not using an algo anyways
      treat_var  = var(outcome[(len+1):(len + treatment_length)])
    )
  
  
  # Treatment stats WITH effect added to the treatment data
  treatment_stats_effect <- dat_effect %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length), #check if this is correct, if it does not need to be first(len)
      treat_mean = mean(outcome_effect[(len+1):(len + treatment_length)]),
      treat_var  = var(outcome_effect[(len+1):(len + treatment_length)])
    )
  
  
  phase_stats <- baseline_stats %>%
    full_join(treatment_stats, by = c("replicate", "Algorithm", "len"))
  
  
  # phase_stats with the treatment effect data included
  phase_stats_effect <- baseline_stats %>%
    full_join(treatment_stats_effect, by = c("replicate", "Algorithm", "len")) 
  
  
  
  # ---------------------------------------------------------
  # Vectorized effect size calculation WITH effect added to the treatment
  # ---------------------------------------------------------
  phase_stats_effect <- phase_stats_effect %>%
    mutate(
      pooled_sd = sqrt(((len - 1) * base_var + (treatment_length - 1) * treat_var) / (len + treatment_length - 2)),
      SMD     = ifelse(pooled_sd > 0, (treat_mean - base_mean) / pooled_sd, NA),
      nu      = len + treatment_length - 2,
      J_nu    = 1 - 3 / (4 * nu - 1),
      g       = J_nu * SMD,
      SE = sqrt((len + treatment_length) / (len * treatment_length) + g^2 / (2 * (len + treatment_length))),
      p_value = 2 * pnorm(-abs(g / SE))
    )
  
  # ---------------------------------------------------------
  # Type I error calculation
  # ---------------------------------------------------------
  reg_results <- dat_reg %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      p_value = {
        if (length(unique(phase)) < 2) {
          NA_real_
        } else {
          model <- glm(outcome ~ phase, family = quasipoisson(link = "log"))
          coefs <- summary(model)$coefficients
          coefs["phase", "Pr(>|t|)"]
        }
      },
      .groups = "drop"
    )
  
  type1_error_stats <- reg_results %>%
    group_by(Algorithm) %>%
    summarize(
      R_eff = sum(!is.na(p_value)),  # effective number of replications
      
      # Proportion of simulations where a valid phase split occurred
      prop_valid = R_eff / n(),
      
      type1_error = mean(p_value < 0.05, na.rm = TRUE),
      
      # MCSE
      MCSE = sqrt(type1_error * (1 - type1_error) / R_eff),
      
      # 95% Monte Carlo CI
      MC_lower = pmax(0, type1_error - 1.96 * MCSE), #making sure that it actually falls in a [0,1] interval
      MC_upper = pmin(1, type1_error + 1.96 * MCSE),
      
      .groups = "drop"
    )
  
  
  
  
  

  
  
  # ---------------------------------------------------------
  # Summary stats including new metrics (carefully consider if you are using the one with or without effect size)
  # ---------------------------------------------------------
  #These are the calculations done without the effect added to the treatment
  summary_stats <- phase_stats %>%
    group_by(Algorithm) %>%
    summarize(
      # Baseline
      base_mean_mean     = mean(base_mean, na.rm = TRUE),
      base_mean_var      = var(base_mean, na.rm = TRUE),
      base_var_mean      = mean(base_var, na.rm = TRUE),
      base_var_var       = var(base_var, na.rm = TRUE),
      
      # Treatment
      treat_mean_mean    = mean(treat_mean, na.rm = TRUE),
      treat_mean_var     = var(treat_mean, na.rm = TRUE),
      treat_var_mean     = mean(treat_var, na.rm = TRUE),
      treat_var_var      = var(treat_var, na.rm = TRUE),
      
      # Relative bias
      base_mean_bias     = mean(base_mean - mu, na.rm = TRUE) / mu,
      base_var_mean_bias = mean(base_var - mu, na.rm = TRUE) / mu,
      
      treat_mean_bias    = mean(treat_mean - mu, na.rm = TRUE) / mu,
      treat_var_mean_bias = mean(treat_var - mu, na.rm = TRUE) / mu

    )
  
  # These are the calculations done using the treatment data WITH an effect added to them
  summary_stats_effect <- phase_stats_effect %>%
    group_by(Algorithm) %>%
    summarize(
      mean_g   = mean(g, na.rm = TRUE),
      mean_SE  = mean(SE, na.rm = TRUE),
      true_SE  = sd(g, na.rm = TRUE),
      SE_bias  = mean(SE, na.rm = TRUE) - sd(g, na.rm = TRUE),
      ci_coverage = mean(abs(g - effect_size) / SE <= 1.96, na.rm = TRUE)
    )
  
  phase_stats %>%
    group_by(Algorithm) %>%
    nest() %>%
    full_join(summary_stats, by = "Algorithm") %>%
    full_join(summary_stats_effect, by = "Algorithm") %>%
    full_join(type1_error_stats, by = "Algorithm")
  

}