setwd("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/")
source("Analyst Algorithms.R")
source("common sim functions.R")

# Generate normal data
gen_norm_dat <- function(replicate_n, n, treatment_length, phi, sigma_squared){
  #for consistency, this is written out again
  total_length = n + treatment_length
  
  sd <- sqrt(sigma_squared)
  
  if (phi == 0) {
    dat <- 5 + rnorm(total_length, sd = sd)
  } else {
    dat <- 5 + as.numeric(
      arima.sim(
        n = total_length,
        model = list(ar = phi),
        sd = sd
      )
    )
  }
  
  data.frame(replicate = replicate_n, outcome = dat)
}

sim_driver <- function(iterations, n, treatment_length, phi, sigma_squared, effect_size, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)

  #the absolute effect which we will add to the treatment group
  absolute_effect = effect_size * sqrt(sigma_squared)
  
  alg_vec <- c(
    "fixed5" = fixed5, #this is the "algorithm" used to generate the reference data
#    "fixed49" = fixed49#,
    kaz10   = kaz10,
    kaz15   = kaz15,
    vva     = vva,
    gl_full = gl_full,
    gl_final= gl_final,
    gl_abs  = gl_abs,
    gl_rel  = gl_rel
  )
  
  total_length = n + treatment_length
  
 
  dat <- map_dfr(
    1:iterations,
    gen_norm_dat,
    n = n,
    treatment_length = treatment_length,
    phi = phi,
    sigma_squared = sigma_squared
  )
  
  
  lengths <- dat %>%
    group_by(replicate) %>%
    group_modify(
      ~ map_dfr(
        .x = alg_vec,
        .f = min_length,
        outcomes = .x$outcome,
        phase_length = n,
        begin_at = 3
      )
    ) %>%
    arrange(replicate) %>%
    gather(key = "Algorithm", value = "length", fixed5:gl_rel)
  
  
  
  #Creating data used for type I error calculation. Adding an indicator variable,
  #which is required to calculate the type I error
  dat_reg <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, Algorithm) %>%
    mutate(
      len = first(length),
      time = row_number(),
      phase = ifelse(time <= len, 0, 1)
    ) %>%
    ungroup()
  
  # Data with the effect added
  dat_effect <- lengths %>%
    full_join(dat, by = "replicate") %>%
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
  
  
  
  
  # Baseline stats
  baseline_stats <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length),
      base_mean = ifelse(len < 101, mean(outcome[1:len]), NA),
      base_var  = ifelse(len < 101, var(outcome[1:len]), NA),
      .groups = "drop"
    )
  
  # Treatment stats
  treatment_stats <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length),
      treat_mean = mean(outcome[(len+1):(len + treatment_length)]),
      treat_var  = var(outcome[(len+1):(len + treatment_length)]),
      .groups = "drop"
    )
  
  
  # Treatment stats WITH effect added
  treatment_stats_effect <- dat_effect %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      len = first(length),
      treat_mean = mean(outcome_effect[(len + 1):(len + treatment_length)]),
      treat_var  = var(outcome_effect[(len + 1):(len + treatment_length)]),
      .groups = "drop"
    )
  
  
  
  phase_stats <- baseline_stats %>%
    full_join(treatment_stats, by = c("replicate", "Algorithm", "len"))
  
  phase_stats_effect <- baseline_stats %>%
    full_join(treatment_stats_effect, by = c("replicate", "Algorithm", "len"))
  
  # ---------------------------------------------------------
  # Vectorized effect size calculation WITH effect added to it
  # ---------------------------------------------------------
  phase_stats_effect <- phase_stats_effect %>%
    mutate(
      pooled_sd = sqrt(((len - 1) * base_var + (treatment_length - 1) * treat_var) / (len + treatment_length - 2)),
      SMD       = ifelse(pooled_sd > 0, (treat_mean - base_mean) / pooled_sd, NA),
      nu        = len + treatment_length - 2,
      J_nu      = 1 - 3 / (4 * nu - 1),
      g         = J_nu * SMD,
      SE = sqrt((len + treatment_length) / (len * treatment_length) + g^2 / (2 * (len + treatment_length))),
      p_value   = 2 * pnorm(-abs(g / SE))
    )
 
  # ---------------------------------------------------------
  #Type I error calculations
  # ---------------------------------------------------------
  
  #Calculate the regression coefficient
  reg_results <- dat_reg %>%
    group_by(replicate, Algorithm) %>%
    summarize(
      p_value = summary(lm(outcome ~ phase))$coefficients["phase", "Pr(>|t|)"], #t-test is already included in lm() summary, so that makes it easy for us
      .groups = "drop"
    )


  
  
    
  #Calculate the proportion of rejections, i.e. the type I error (the number of times the 
  #regression coefficient in our model is considered significant)
  type1_error_stats <- reg_results %>%
    group_by(Algorithm) %>%
    summarize(
      R_eff = sum(!is.na(p_value)),  # effective number of replications
      
      # Proportion of simulations which is valid
      prop_valid = R_eff / n(),
      
      type1_error = mean(p_value < 0.05, na.rm = TRUE),
      
      # MCSE
      MCSE = sqrt(type1_error * (1 - type1_error) / R_eff),
      
      # 95% Monte Carlo CI
      MC_lower = pmax(0, type1_error - 1.96 * MCSE), #making sure it falls in a [0,1] interval
      MC_upper = pmin(1, type1_error + 1.96 * MCSE),
      
      .groups = "drop"
    )
  
  
  
  
  
  
   
  # ---------------------------------------------------------
  # Summary stats without effect
  # ---------------------------------------------------------
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
      base_mean_bias     = mean(base_mean - 5, na.rm = TRUE) / 5,
      base_var_mean_bias = mean(base_var - sigma_squared, na.rm = TRUE),
      
      treat_mean_bias    = mean(treat_mean - 5, na.rm = TRUE) / 5,
      treat_var_mean_bias = mean(treat_var - sigma_squared, na.rm = TRUE),

      .groups = "drop"
    )
  
  
  # Effect size metrics
  summary_stats_effect <- phase_stats_effect %>%
    group_by(Algorithm) %>%
    summarize(
      # Effect size metrics
      mean_g             = mean(g, na.rm = TRUE),
      mean_SE            = mean(SE, na.rm = TRUE),
      
      # SE bias and CI coverage
      true_SE     = sd(g, na.rm = TRUE),
      SE_bias     = mean(SE, na.rm = TRUE) - sd(g, na.rm = TRUE),
      ci_coverage = mean(abs(g - effect_size) / SE <= 1.96, na.rm = TRUE)
    )
  
  
  
  phase_stats %>%
    group_by(Algorithm) %>%
    nest() %>%
    full_join(summary_stats, by = "Algorithm") %>%
    full_join(summary_stats_effect, by = "Algorithm") %>%
    full_join(type1_error_stats, by = "Algorithm")
}