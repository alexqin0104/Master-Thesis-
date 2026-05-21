source("Analyst Algorithms.R")
source("common sim functions.R")

library(ARPobservation)


generate_gpp <- function(replicate_n, mu_a, shape, phase_length){
  data.frame(
    replicate = replicate_n,
    outcome = c(event_counting(
      r_behavior_stream(
        n = phase_length,
        mu = 0,
        lambda = 1/mu_a,
        F_event = F_exp(),
        F_interim = F_gam(shape),
        stream_length = 1
      )
    ))
  )
}

sim_driver <- function(iterations, n, treatment_length, mu_a, shape, effect_size){
  
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
  dat <- map_dfr(
    1:iterations,
    generate_gpp,
    mu_a = mu_a,
    shape = shape,
    phase_length = total_length
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
    gather(key = "alg", value = "length", fixed5:gl_rel)
  
  

  #This is the variance as calculated in the original code. Hence, it was not changed.
  var_i = (mu_a/shape) + (shape^2 - 1)/(6 * shape^2)
  absolute_effect = effect_size * sqrt(var_i)
  
  #data with the effect added
  dat_effect <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, alg) %>%
    mutate(
      len = first(length),
      outcome_effect = ifelse(
        row_number() >= (len + 1) & row_number() <= (len + treatment_length),
        outcome + absolute_effect,
        outcome
      )
    ) %>%
    ungroup()
  
  
  #data with indicator, for type I calculation
  dat_reg <- lengths %>%
    full_join(dat, by = "replicate") %>%
    group_by(replicate, alg) %>%
    mutate(
      len = first(length),
      time = row_number(),
      phase = ifelse(time <= len, 0, 1)
    ) %>%
    ungroup()
  
  
  # baseline stats
  baseline_stats <- lengths %>%
    full_join(dat) %>%
    group_by(replicate, alg) %>%
    summarize(
      len = first(length),
      base_mean = ifelse(len < 101, mean(outcome[1:len]), NA),
      base_var  = ifelse(len < 101, var(outcome[1:len]), NA)
    )
  
  # treatment stats
  treatment_stats <- lengths %>%
    full_join(dat) %>%
    group_by(replicate, alg) %>%
    summarize(
      len = first(length),
      treat_mean = mean(outcome[(len+1):(len + treatment_length)]),
      treat_var  = var(outcome[(len+1):(len + treatment_length)])
    )
  
  # Treatment stats WITH effect added
  treatment_stats_effect <- dat_effect %>%
    group_by(replicate, alg) %>%
    summarize(
      len = first(length), #check if this is correct, if it does not need to be first(len)
      treat_mean = mean(outcome_effect[(len + 1):(len + treatment_length)]),
      treat_var  = var(outcome_effect[(len + 1):(len + treatment_length)])
    )
  
  
  
  phase_stats <- baseline_stats %>%
    full_join(treatment_stats, by = c("replicate", "alg", "len"))
  
  #Phase stats with the effect added
  phase_stats_effect <- baseline_stats %>%
    full_join(treatment_stats_effect, by = c("replicate", "alg", "len"))
  
  # ---------------------------------------------------------
  # Vectorized effect size calculation, WITH treatment effect added
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
  
  
  #---------------------------------------------------------
  #Type I error calculation
  #---------------------------------------------------------
  reg_results <- dat_reg %>%
    group_by(replicate, alg) %>%
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
    group_by(alg) %>%
    summarize(
      R_eff = sum(!is.na(p_value)),  # effective number of replications
      
      #From R_eff we can easily derive the proportion of baselines where the algorithm actually considers it stable.
      #This is because we created a phase indicator, and based on this phase indicator we conducted the testing. Now,
      #we simply need check the proportion of valid tests, this is the same as the proportion of baselines considered
      #stable by the algorithm.
      #We ended up not discussing the proportion of valid outcomes, since either all or the overwhelming number of outcomes
      #turned out stable, so there is no real difference between the different algorithm/data generating setups
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
    group_by(alg) %>%
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
      
      # relative bias
      base_mean_bias     = mean(base_mean - mu_a, na.rm = TRUE) / mu_a,
      base_var_mean_bias = mean(base_var - var_i, na.rm = TRUE) / var_i,
      
      treat_mean_bias    = mean(treat_mean - mu_a, na.rm = TRUE) / mu_a,
      treat_var_mean_bias = mean(treat_var - var_i, na.rm = TRUE) / var_i
      
    )
  
  # Effect size metrics
  
  summary_stats_effect <- phase_stats_effect %>%
    group_by(alg) %>%
    summarize(
      # Effect size metrics
      mean_g             = mean(g, na.rm = TRUE),
      # mean_SMD           = mean(g, na.rm = TRUE),
      mean_SE            = mean(SE, na.rm = TRUE),
      
      # SE bias and CI coverage
      true_SE     = sd(g, na.rm = TRUE),
      SE_bias     = mean(SE, na.rm = TRUE) - sd(g, na.rm = TRUE),
      ci_coverage = mean(abs(g - effect_size) / SE <= 1.96, na.rm = TRUE)
    )
  
  
  
  phase_stats %>%
    group_by(alg) %>%
    nest() %>%
    full_join(summary_stats, by = "alg") %>%
    full_join(summary_stats_effect, by = "alg") %>%
    full_join(type1_error_stats, by = "alg")
  
}

