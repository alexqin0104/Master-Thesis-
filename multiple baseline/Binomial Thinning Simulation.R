library(tidyverse)
library(nlme)
library(scdhlm)
source("Analyst Algorithms.R")
source("common sim functions.R")

# ============================================================
# 1. Generate one Poisson binomial thinning observation (point-by-point)
# ============================================================
gen_one_pois <- function(prev_value, mu, phi) {
  if (is.na(prev_value) || prev_value < 0) prev_value <- 0
  prev_value <- round(prev_value)
  
  if (prev_value == 0) return(rpois(1, lambda = mu))
  
  ar_vec <- ifelse(mu - phi * prev_value >= 0, phi, mu / prev_value)
  ar_vec <- pmin(pmax(ar_vec, 0), 1)
  lambda <- pmax(0, mu - ar_vec * prev_value)
  
  rbinom(1, size = prev_value, prob = ar_vec) + rpois(1, lambda = lambda)
}

# ============================================================
# 2. Check whether treatment response has occurred (Joo et al., 2018)
# ============================================================
has_treatment_response <- function(baseline_data, treatment_data) {
  if (length(treatment_data) < 3) return(FALSE)
  
  n_b <- length(baseline_data)
  n_t <- length(treatment_data)
  
  v_b <- var(baseline_data)
  v_t <- var(treatment_data)
  if (is.na(v_t)) v_t <- 0
  
  pooled_sd <- sqrt(((n_b - 1) * v_b + (n_t - 1) * v_t) / (n_b + n_t - 2))
  
  if (is.na(pooled_sd) || pooled_sd == 0) return(FALSE)
  
  cond1 <- (mean(treatment_data) - mean(baseline_data)) > 2 * pooled_sd
  cond2 <- (mean(tail(treatment_data, 3)) - mean(baseline_data)) > 2 * pooled_sd
  
  return(cond1 & cond2)
}

# ============================================================
# 3. Generate data for a single case (treatment data always generated
#    even when baseline ceiling is hit)
# ============================================================
generate_case_data <- function(alg, mu_i, phi, absolute_effect,
                               min_start, max_bl = 40, treatment_length = 15) {
  
  # Generate initial baseline
  data <- numeric(min_start)
  data[1] <- rpois(1, lambda = mu_i)
  if (min_start > 1) {
    for (t in 2:min_start) {
      data[t] <- gen_one_pois(data[t - 1], mu_i, phi)
    }
  }
  
  # Baseline stability check (response-guided only)
  bl_len <- min_start
  hit_ceiling <- FALSE
  reached_stability <- FALSE
  
  if (!is.character(alg) || alg != "fixed") {
    repeat {
      is_stable <- do.call(alg, list(outcomes = data[1:bl_len]))
      if (is_stable) {
        reached_stability <- TRUE
        break
      }
      if (bl_len >= max_bl) {
        hit_ceiling <- TRUE
        break
      }
      bl_len <- bl_len + 1
      data[bl_len] <- gen_one_pois(data[bl_len - 1], mu_i, phi)
    }
  } else {
    reached_stability <- TRUE  # fixed design does not require stability check
  }
  
  # === Treatment data always generated regardless of ceiling hit ===
  tx_data <- numeric(0)
  responded <- FALSE
  effect_achieved_day <- NA
  
  for (t in 1:treatment_length) {
    new_pt <- gen_one_pois(data[length(data)], mu_i, phi) + absolute_effect
    tx_data <- c(tx_data, new_pt)
    data <- c(data, new_pt)
    
    # Response check: skipped for fixed design, applied for response-guided
    if (!is.character(alg) || alg != "fixed") {
      if (t >= 3 && has_treatment_response(data[1:bl_len], tx_data)) {
        responded <- TRUE
        effect_achieved_day <- bl_len + t
        break
      }
    }
  }
  
  # Finalize effect_achieved_day
  if (is.character(alg) && alg == "fixed") {
    effect_achieved_day <- bl_len + treatment_length
    responded <- TRUE
  } else if (!responded) {
    effect_achieved_day <- bl_len + length(tx_data)
  }
  
  list(
    baseline = data[1:bl_len],
    tx_data = tx_data,
    bl_len = bl_len,
    tx_len = length(tx_data),
    responded = responded,
    hit_ceiling = hit_ceiling,
    reached_stability = reached_stability,
    effect_achieved_day = effect_achieved_day
  )
}

# ============================================================
# 4. Fixed condition: baseline lengths staggered by design (5, 8, 11, 14)
# ============================================================
generate_fixed_case_data <- function(case_id, mu_i, phi, absolute_effect,
                                     min_bl = 5, min_stagger = 3,
                                     treatment_length = 15) {
  
  bl_len <- min_bl + (case_id - 1) * min_stagger
  
  data <- numeric(bl_len)
  data[1] <- rpois(1, lambda = mu_i)
  if (bl_len > 1) {
    for (t in 2:bl_len) {
      data[t] <- gen_one_pois(data[t - 1], mu_i, phi)
    }
  }
  
  tx_data <- numeric(treatment_length)
  prev_v <- data[bl_len]
  for (t in 1:treatment_length) {
    new_v <- gen_one_pois(prev_v, mu_i, phi) + absolute_effect
    tx_data[t] <- new_v
    prev_v <- new_v
  }
  
  list(
    baseline = data,
    tx_data = tx_data,
    bl_len = bl_len,
    tx_len = treatment_length,
    responded = NA,
    hit_ceiling = FALSE,
    reached_stability = TRUE,
    effect_achieved_day = bl_len + treatment_length
  )
}

# ============================================================
# 5. Simulate one replicate with 4 cases and random intercepts
# ============================================================
sim_one_replicate_mbd <- function(alg, mu, phi, absolute_effect,
                                  n_cases = 4, min_bl = 5, min_stagger = 3,
                                  max_bl = 40, treatment_length = 15) {
  
  # Case-level random intercepts (sigma_r0 = sqrt(mu))
  sigma_r0 <- sqrt(mu)
  case_intercepts <- rnorm(n_cases, 0, sigma_r0)
  mu_i_vec <- pmax(0.5, mu + case_intercepts)
  
  # Fixed condition
  if (is.character(alg) && alg == "fixed") {
    case_results <- lapply(1:n_cases, function(i) {
      generate_fixed_case_data(i, mu_i_vec[i], phi, absolute_effect,
                               min_bl, min_stagger, treatment_length)
    })
    for (i in seq_along(case_results)) case_results[[i]]$mu_i <- mu_i_vec[i]
    return(case_results)
  }
  
  # Response-guided condition
  case_results <- vector("list", n_cases)
  prev_effect_day <- 0
  
  for (case_id in 1:n_cases) {
    current_min_bl <- max(min_bl, prev_effect_day + min_stagger)
    # Cap current_min_bl at max_bl to avoid exceeding ceiling
    current_min_bl <- min(current_min_bl, max_bl)
    
    cr <- generate_case_data(alg, mu_i_vec[case_id], phi, absolute_effect,
                             current_min_bl, max_bl, treatment_length)
    cr$mu_i <- mu_i_vec[case_id]
    
    case_results[[case_id]] <- cr
    prev_effect_day <- cr$effect_achieved_day
  }
  
  case_results
}

# ============================================================
# 6. Convert case results to long format
# ============================================================
prepare_lme_data <- function(case_results) {
  long_data <- bind_rows(lapply(seq_along(case_results), function(i) {
    cr <- case_results[[i]]
    data.frame(
      case = factor(i),
      time = 1:(cr$bl_len + cr$tx_len),
      phase = c(rep("baseline", cr$bl_len), rep("treatment", cr$tx_len)),
      outcome = c(cr$baseline, cr$tx_data)
    )
  }))
  
  long_data$phase <- factor(long_data$phase, levels = c("baseline", "treatment"))
  long_data
}

# ============================================================
# 7. LMM analysis (REML + AR(1) + fallback chain) with g_REML
# ============================================================
analyze_with_lme <- function(case_results, true_effect, baseline_sd) {
  
  long_data <- prepare_lme_data(case_results)
  
  # Data sufficiency check (only triggered in extreme cases)
  if (nrow(long_data) < 10 || 
      sum(long_data$phase == "baseline") < 5 ||
      sum(long_data$phase == "treatment") < 5) {
    return(list(
      estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, fit_method = "insufficient_data",
      g_est = NA, g_se = NA
    ))
  }
  
  # === Try 1: REML + AR(1) ===
  model <- tryCatch({
    suppressWarnings(
      lme(outcome ~ phase,
          random = ~ 1 | case,
          correlation = corAR1(form = ~ time | case),
          data = long_data,
          method = "REML",
          control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100,
                               returnObject = TRUE))
    )
  }, error = function(e) NULL)
  fit_method <- "REML_AR1"
  has_ar1 <- TRUE
  
  # === Try 2: REML without AR(1) ===
  if (is.null(model)) {
    model <- tryCatch({
      suppressWarnings(
        lme(outcome ~ phase,
            random = ~ 1 | case,
            data = long_data,
            method = "REML",
            control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100,
                                 returnObject = TRUE))
      )
    }, error = function(e) NULL)
    fit_method <- "REML_noAR"
    has_ar1 <- FALSE
  }
  
  # === Try 3: ML + AR(1) ===
  if (is.null(model)) {
    model <- tryCatch({
      suppressWarnings(
        lme(outcome ~ phase,
            random = ~ 1 | case,
            correlation = corAR1(form = ~ time | case),
            data = long_data,
            method = "ML",
            control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100,
                                 returnObject = TRUE))
      )
    }, error = function(e) NULL)
    fit_method <- "ML_AR1"
    has_ar1 <- TRUE
  }
  
  # === Try 4: ML without AR(1) ===
  if (is.null(model)) {
    model <- tryCatch({
      suppressWarnings(
        lme(outcome ~ phase,
            random = ~ 1 | case,
            data = long_data,
            method = "ML",
            control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100,
                                 returnObject = TRUE))
      )
    }, error = function(e) NULL)
    fit_method <- "ML_noAR"
    has_ar1 <- FALSE
  }
  
  if (is.null(model)) {
    return(list(
      estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, fit_method = "all_failed",
      g_est = NA, g_se = NA
    ))
  }
  
  # Extract fixed effect estimate (gamma_11)
  ft <- summary(model)$tTable
  if (!"phasetreatment" %in% rownames(ft)) {
    return(list(
      estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, fit_method = fit_method,
      g_est = NA, g_se = NA
    ))
  }
  
  estimate <- ft["phasetreatment", "Value"]
  se <- ft["phasetreatment", "Std.Error"]
  df <- ft["phasetreatment", "DF"]
  
  t_crit <- qt(0.975, df)
  ci_lower <- estimate - t_crit * se
  ci_upper <- estimate + t_crit * se
  
  # === g_REML: design-comparable SMD ===
  # Requires REML estimation with AR(1) correlation structure
  # r_const = c(1, 0, 1): random intercept var + AR phi (=0) + level-1 sigma^2
  # Consistent with the official test (Laski example)
  g_est <- NA
  g_se  <- NA
  
  if (fit_method == "REML_AR1") {
    g_result <- tryCatch({
      suppressWarnings(
        scdhlm::g_REML(
          m_fit   = model,
          p_const = c(0, 1),
          r_const = c(1, 0, 1),
          returnModel = FALSE
        )
      )
    }, error = function(e) NULL)
    
    if (!is.null(g_result) && !is.null(g_result$g_AB)) {
      g_est <- g_result$g_AB
      g_se  <- sqrt(g_result$V_g_AB)
    }
  }
  
  list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    converged = TRUE,
    fit_method = fit_method,
    g_est = g_est,
    g_se = g_se
  )
}

# ============================================================
# 8. Simulation driver: 4 cases with ceiling/stability/g_REML tracking
# ============================================================
sim_driver <- function(iterations, n, treatment_length, mu, phi, effect_size) {
  
  alg_vec <- list(
    "fixed"    = "fixed",
    "kaz10"    = kaz10,
    "kaz15"    = kaz15,
    "vva"      = vva,
    "gl_full"  = gl_full,
    "gl_final" = gl_final,
    "gl_abs"   = gl_abs,
    "gl_rel"   = gl_rel
  )
  
  abs_effect <- effect_size * sqrt(mu)
  true_effect <- abs_effect
  baseline_sd <- sqrt(mu)  # reference SD for g standardization
  
  res_list <- map_dfr(names(alg_vec), function(alg_name) {
    
    alg_func <- alg_vec[[alg_name]]
    
    reps <- map_dfr(1:iterations, function(rep_i) {
      
      case_results <- sim_one_replicate_mbd(
        alg = alg_func, mu = mu, phi = phi,
        absolute_effect = abs_effect,
        n_cases = 4,
        max_bl = n, treatment_length = treatment_length
      )
      
      result <- analyze_with_lme(case_results, true_effect, baseline_sd)
      
      bl_lens <- sapply(case_results, function(x) x$bl_len)
      responded <- sapply(case_results, function(x) x$responded)
      hit_ceiling_per_case <- sapply(case_results, function(x) x$hit_ceiling)
      reached_stability_per_case <- sapply(case_results, function(x) x$reached_stability)
      
      data.frame(
        estimate = result$estimate,
        se = result$se,
        ci_lower = result$ci_lower,
        ci_upper = result$ci_upper,
        converged = result$converged,
        fit_method = result$fit_method,
        g_est = result$g_est,
        g_se = result$g_se,
        bl1 = bl_lens[1], bl2 = bl_lens[2], bl3 = bl_lens[3], bl4 = bl_lens[4],
        responded_rate = mean(responded, na.rm = TRUE),
        any_hit_ceiling = any(hit_ceiling_per_case),
        all_hit_ceiling = all(hit_ceiling_per_case),
        prop_hit_ceiling = mean(hit_ceiling_per_case),
        prop_reached_stability = mean(reached_stability_per_case),
        # Stagger preservation metrics:
        # Ideal stagger: case 4 baseline is 9 longer than case 1 (= 3 * min_stagger)
        # stagger_ratio = 1.0: perfect stagger (Fixed design)
        # stagger_ratio = 0.0: stagger fully collapsed (all cases same baseline length)
        stagger_diff  = bl_lens[4] - bl_lens[1],
        stagger_ratio = (bl_lens[4] - bl_lens[1]) / 9,
        stagger_collapsed = (bl_lens[4] - bl_lens[1] < 1)
      )
    })
    
    reps$Algorithm <- alg_name
    reps
  })
  
  # Summarize across replicates
  res_list %>%
    group_by(Algorithm) %>%
    summarize(
      convergence_rate = mean(converged, na.rm = TRUE),
      
      reml_ar1_rate    = mean(fit_method == "REML_AR1", na.rm = TRUE),
      reml_noar_rate   = mean(fit_method == "REML_noAR", na.rm = TRUE),
      ml_fallback_rate = mean(fit_method %in% c("ML_AR1", "ML_noAR"), na.rm = TRUE),
      
      # Algorithm viability indicators
      any_hit_ceiling_rate    = mean(any_hit_ceiling, na.rm = TRUE),
      all_hit_ceiling_rate    = mean(all_hit_ceiling, na.rm = TRUE),
      mean_prop_hit_ceiling   = mean(prop_hit_ceiling, na.rm = TRUE),
      mean_prop_stability     = mean(prop_reached_stability, na.rm = TRUE),
      
      # Stagger preservation (multiple-baseline design integrity)
      mean_stagger_diff       = mean(stagger_diff, na.rm = TRUE),
      mean_stagger_ratio      = mean(stagger_ratio, na.rm = TRUE),
      stagger_collapsed_rate  = mean(stagger_collapsed, na.rm = TRUE),
      
      # Treatment effect estimates (gamma_11)
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - true_effect, na.rm = TRUE),
      relative_bias = if (true_effect == 0) NA_real_ else mean((estimate - true_effect) / true_effect, na.rm = TRUE),
      rmse = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
      mean_se = mean(se, na.rm = TRUE),
      empirical_se = sd(estimate, na.rm = TRUE),
      se_bias = mean(se, na.rm = TRUE) - sd(estimate, na.rm = TRUE),
      ci_coverage = mean(
        ci_lower <= true_effect & true_effect <= ci_upper,
        na.rm = TRUE
      ),
      reject_null = mean(
        ci_lower > 0 | ci_upper < 0,
        na.rm = TRUE
      ),
      
      # Design-comparable SMD (g_REML)
      g_rate_computed = mean(!is.na(g_est)),
      mean_g          = mean(g_est, na.rm = TRUE),
      bias_g          = mean(g_est - effect_size, na.rm = TRUE),  # effect_size is the true SMD
      rmse_g          = sqrt(mean((g_est - effect_size)^2, na.rm = TRUE)),
      mean_g_se       = mean(g_se, na.rm = TRUE),
      empirical_g_se  = sd(g_est, na.rm = TRUE),
      
      # Baseline lengths
      mean_bl1 = mean(bl1, na.rm = TRUE),
      mean_bl2 = mean(bl2, na.rm = TRUE),
      mean_bl3 = mean(bl3, na.rm = TRUE),
      mean_bl4 = mean(bl4, na.rm = TRUE),
      
      response_rate = mean(responded_rate, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    select(Algorithm, convergence_rate,
           reml_ar1_rate, reml_noar_rate, ml_fallback_rate,
           any_hit_ceiling_rate, all_hit_ceiling_rate,
           mean_prop_hit_ceiling, mean_prop_stability,
           mean_stagger_diff, mean_stagger_ratio, stagger_collapsed_rate,
           response_rate,
           mean_bl1, mean_bl2, mean_bl3, mean_bl4,
           mean_estimate, bias, relative_bias, rmse,
           mean_se, empirical_se, se_bias, ci_coverage, reject_null,
           g_rate_computed, mean_g, bias_g, rmse_g, mean_g_se, empirical_g_se)
}
