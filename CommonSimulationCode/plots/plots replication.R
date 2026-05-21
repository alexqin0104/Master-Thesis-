rm(list = ls())
load("C:/Users/eloua/Desktop/master aangewezen/replication/n5/binomthin.Rdata")
load("C:/Users/eloua/Desktop/master aangewezen/replication/n5/gpp.Rdata")
load("C:/Users/eloua/Desktop/master aangewezen/replication/n5/normCorrected.Rdata")

library(tidyverse)


alg_order <- c(
  "kaz10",
  "kaz15",
  "gl_full",
  "gl_final",
  "gl_abs",
  "gl_rel",
  "vva"
)

binomial_thinning_results <- binomial_thinning_results %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_order))

gamma_point_process_results <- gamma_point_process_results %>%
  mutate(alg = factor(alg, levels = alg_order))

normal_errors_results <- normal_errors_results %>%
  mutate(Algorithm = factor(Algorithm, levels = alg_order))


#everything is already there, nothing needs to be calculated, only the plots need to be created

head(binomial_thinning_results)




#forgot to add the label used in the original paper to the plots, so had to add it afterwards.
#therefore the folllowing is used:
alg_title_labeller <- labeller(
  Algorithm = c(
    "gl_abs"   = "GSAbs",
    "gl_final" = "GSFinal",
    "gl_full"  = "GSFull",
    "gl_rel"   = "GSRel",
    "kaz10"    = "Kaz10",
    "kaz15"    = "Kaz15",
    "vva"      = "VVA"
  ),
  alg = c(
    "gl_abs"   = "GSAbs",
    "gl_final" = "GSFinal",
    "gl_full"  = "GSFull",
    "gl_rel"   = "GSRel",
    "kaz10"    = "Kaz10",
    "kaz15"    = "Kaz15",
    "vva"      = "VVA"
  )
)




colnames(binomial_thinning_results)


#first graph (Poisson, phi = 0, varying mu)


#---------------------------------------------------------
# 1 — Unnest all simulation results
#---------------------------------------------------------
df <- binomial_thinning_results %>%
  select(mu, phi, Algorithm, data) %>%
  unnest(data)  # brings replicate, base_mean, base_var, len into the table



# 2. Keep only the independent Poisson case: phi == 0
df_indep <- df %>%
  filter(phi == 0)



#---------------------------------------------------------
# 2 — Compute cumulative percentages up to len = 20
#     BUT denominator is total replicates for the condition
#---------------------------------------------------------
df_cum <- df_indep %>%
  group_by(Algorithm, mu) %>%
  mutate(total_reps = n()) %>%        # denominator = ALL replicates
  # Count how many reached stability at or before each len value
  filter(len <= 20) %>%               # only for plotting
  count(len, total_reps, .groups = "drop_last") %>%   # n = count at each len
  
  # Make sure ALL len = 1:20 appear, even if no data
  # complete len sequence but KEEP total_reps
  tidyr::complete(len = 1:20,
                  fill = list(n = 0)) %>%
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  ungroup()




#some bars were missing due to no observations, so this needs to be added so that
#there is a bar for all values. This is the same for all other bar charts
#created here.


#---------------------------------------------------------
# 3 — Plot cumulative percentage bar charts
#---------------------------------------------------------
ggplot(df_cum,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(mu))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "mean"
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)



#####
#poisson, mu = 15, varying phi


#---------------------------------------------------------
# 1 — Unnest simulation data
#---------------------------------------------------------
df <- binomial_thinning_results %>%
  select(mu, phi, Algorithm, data) %>%
  unnest(data)


#---------------------------------------------------------
# 2 — Keep only mu = 15
#---------------------------------------------------------
df_mu15 <- df %>%
  filter(mu == 15)


#---------------------------------------------------------
# 3 — Compute cumulative percent for each Algorithm × phi
#---------------------------------------------------------
df_cum_phi <- df_mu15 %>%
  group_by(Algorithm, phi) %>%
  
  # denominator = ALL replicates for this condition
  mutate(total_reps = n()) %>%
  
  # only plot up to len=20 but still count full denominator
  filter(len <= 20) %>%
  
  count(len, total_reps, name = "n") %>%
  
  # ensure all len values appear
  complete(len = 1:20, fill = list(n = 0)) %>%
  
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  ungroup()

ggplot(df_cum_phi,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(phi))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "phi"
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)



#####
#gamma point, fixing dispersion, varying mu (this was originally broken, since dispersion = 1 didn't exist berfore)

#---------------------------------------------------------
# 1 — Unnest gamma point-process simulation results
#---------------------------------------------------------
df_gamma <- gamma_point_process_results %>%
  select(mu_a, shape, alg, data) %>%
  unnest(data)    # brings replicate, base_mean, base_var, len


#---------------------------------------------------------
# 2 — Fix dispersion = 1  (i.e. shape = 1)
#---------------------------------------------------------
df_gamma_disp1 <- df_gamma %>%
  filter(shape == 1, mu_a %in% c(5, 15, 25))


#---------------------------------------------------------
# 3 — Compute cumulative percentage (len ≤ 20, full denominator)
#---------------------------------------------------------
df_gamma_cum <- df_gamma_disp1 %>%
  group_by(alg, mu_a) %>%
  mutate(total_reps = n()) %>%
  filter(len <= 20) %>%
  count(len, total_reps, name = "n") %>%
  complete(len = 1:20, fill = list(n = 0)) %>%
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  ungroup()


#---------------------------------------------------------
# 4 — Plot (grouped by mean level)
#---------------------------------------------------------
ggplot(df_gamma_cum,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(mu_a))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "Mean level (mu_a)"
  ) +
  
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ alg, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)


#####
#gamma point, varying dispersion
#---------------------------------------------------------
# 1 — Unnest gamma point-process simulation results
#---------------------------------------------------------
df_gamma <- gamma_point_process_results %>%
  select(mu_a, shape, alg, data) %>%
  unnest(data)    # brings replicate, base_mean, base_var, len


#---------------------------------------------------------
# 2 — Fix mean = 15
#---------------------------------------------------------
df_gamma_mu15 <- df_gamma %>%
  filter(mu_a == 15) %>%
  mutate(dispersion = 1 / shape)


#---------------------------------------------------------
# 3 — Compute cumulative percentage (len ≤ 20, full denominator)
#---------------------------------------------------------
df_gamma_cum <- df_gamma_mu15 %>%
  group_by(alg, dispersion) %>%
  mutate(total_reps = n()) %>%
  filter(len <= 20) %>%
  count(len, total_reps, name = "n") %>%
  complete(len = 1:20, fill = list(n = 0)) %>%
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  ungroup()


ggplot(
  df_gamma_cum,
  aes(
    x = factor(len),
    y = cum_pct,
    fill = factor(dispersion)
  )
) +
  geom_col(position = "dodge") +
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "Dispersion"
  ) +
  scale_x_discrete(
    breaks = as.character(seq(1, 20, by = 2))
  ) +
  scale_fill_discrete(
    labels = function(x) sprintf("%.2f", as.numeric(x))
  ) +
  facet_wrap(
    ~ alg,
    ncol = 3,
    labeller = alg_title_labeller
  ) +
  theme_minimal(base_size = 14)




#normal, sigma squared = 1, varying autocorrelation
#####


#---------------------------------------------------------
# 1 — Unnest normal-errors simulation results
#---------------------------------------------------------
df_norm <- normal_errors_results %>%
  select(phi, sigma_squared, Algorithm, data) %>%
  unnest(data)   # adds replicate, base_mean, base_var, len


#---------------------------------------------------------
# 2 — Fix variance = 1
#---------------------------------------------------------
df_norm_var1 <- df_norm %>%
  filter(sigma_squared == 1)


#---------------------------------------------------------
# 3 — Compute cumulative percentage (len ≤ 20, full denominator)
#---------------------------------------------------------
df_norm_cum <- df_norm_var1 %>%
  group_by(Algorithm, phi) %>%
  
  # total replicates for the condition
  mutate(total_reps = n()) %>%
  
  # keep len ≤ 20 for plotting but denominator stays full
  filter(len <= 20) %>%
  
  count(len, total_reps, name = "n") %>%
  
  # ensure lengths 1 through 20 appear
  complete(len = 1:20, fill = list(n = 0)) %>%
  
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  
  ungroup()


ggplot(df_norm_cum,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(phi))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "phi"
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)


#normal, varying sigma squared, fixing phi
######


#---------------------------------------------------------
# 1 — Unnest normal-errors simulation results
#---------------------------------------------------------
df_norm <- normal_errors_results %>%
  select(phi, sigma_squared, Algorithm, data) %>%
  unnest(data)   # brings replicate, base_mean, base_var, len


#---------------------------------------------------------
# 2 — Fix autocorrelation = 0
#---------------------------------------------------------
df_norm_phi0 <- df_norm %>%
  filter(phi == 0)


#---------------------------------------------------------
# 3 — Compute cumulative percentage (len ≤ 20, full denominator)
#---------------------------------------------------------
df_norm_cum <- df_norm_phi0 %>%
  group_by(Algorithm, sigma_squared) %>%
  
  mutate(total_reps = n()) %>%            # denominator = ALL replicates
  
  filter(len <= 20) %>%                   # only using these in plot
  
  count(len, total_reps, name = "n") %>%  # number with len == this value
  
  complete(len = 1:20, fill = list(n = 0)) %>%  # ensure bars exist for 1–20
  
  arrange(len) %>%
  fill(total_reps, .direction = "downup") %>%
  
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  
  ungroup()

ggplot(df_norm_cum,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(sigma_squared))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "Variance"
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)



#relative bias, poisson distribution
#####

#---------------------------------------------------------
# 1 — Prepare data
#    Unnest + keep Poisson results only
#---------------------------------------------------------
df_pois <- binomial_thinning_results %>%
  select(mu, phi, Algorithm, base_mean_bias) %>% 
  mutate(
    relative_bias = base_mean_bias   # it is actually already calculated
  )


#---------------------------------------------------------
# 2 — Summarize by condition (Algorithm × mu × phi)
#    We need the mean relative bias over replicates
#---------------------------------------------------------
df_bias <- df_pois %>%
  group_by(Algorithm, mu, phi) %>%
  summarise(
    rel_bias_mean = mean(relative_bias, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_bias,
       aes(x = mu,
           y = rel_bias_mean,
           color = factor(phi),
           group = phi)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  # dashed ±5% reference lines
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "autocorrelation",
    values = c(
      "0"   = "#2E1A47",
      "0.2" = "#1BA098",
      "0.4" = "#F2D849"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_bias$mu))) +
  
  labs(
    x = "Mean Level",
    y = "Relative Mean Bias"
  ) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)




#gamma point, with the missing dispersion = 1 added
#####


#---------------------------------------------------------
# 1 — Keep relevant gamma results
#---------------------------------------------------------
df_gamma <- gamma_point_process_results %>%
  select(mu_a, shape, alg, base_mean_bias) %>%
  mutate(
    relative_bias = base_mean_bias,
    shape_chr = case_when(
      shape == 2/3 ~ "2/3",
      TRUE ~ as.character(shape)
    )
  )

#---------------------------------------------------------
# 2 — Summarize bias by alg × mean × dispersion
#---------------------------------------------------------
df_gamma_bias <- df_gamma %>%
  group_by(alg, mu_a, shape_chr) %>%
  summarise(
    rel_bias_mean = mean(relative_bias, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_gamma_bias,
       aes(x = mu_a,
           y = rel_bias_mean,
           color = shape_chr,
           group = shape_chr)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "dispersion (shape)",
    values = c(
      "0.4" = "#2E1A47",
      "2/3" = "#1BA098",
      "1"   = "red",
      "1.5" = "#F2D849",
      "2.5" = "#E07A5F"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_gamma_bias$mu_a))) +
  
  labs(
    x = "Mean Level",
    y = "Relative Mean Bias"
  ) +
  
  facet_wrap(~ alg, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)





#normal distribution
#####
library(dplyr)
library(ggplot2)

#---------------------------------------------------------
# 1 — Keep relevant normal error results
#---------------------------------------------------------
df_normal <- normal_errors_results %>%
  select(phi, sigma_squared, Algorithm, base_mean_bias) %>%
  mutate(
    relative_bias = base_mean_bias,
    sigma_chr = as.character(sigma_squared)  # Convert variance to character for color mapping
  )

#---------------------------------------------------------
# 2 — Summarize bias by Algorithm × autocorrelation × variance
#---------------------------------------------------------
df_normal_bias <- df_normal %>%
  group_by(Algorithm, phi, sigma_chr) %>%
  summarise(
    rel_bias_mean = mean(relative_bias, na.rm = TRUE),
    .groups = "drop"
  )

#---------------------------------------------------------
# 3 — Plot
#---------------------------------------------------------
ggplot(df_normal_bias,
       aes(x = phi,
           y = rel_bias_mean,
           color = sigma_chr,
           group = sigma_chr)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "Variance (sigma^2)",
    values = c(
      "0.25" = "#2E1A47",
      "1"   = "#1BA098",
      "2.25"   = "#F2D849"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_normal_bias$phi))) +
  
  labs(
    x = "Autocorrelation (phi)",
    y = "Relative Mean Bias"
  ) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)











#bias in variance
#
#relative bias, poisson distribution
#####

#---------------------------------------------------------
# 1 — Prepare data
#    Unnest + keep Poisson results only
#---------------------------------------------------------
df_pois <- binomial_thinning_results %>%
  select(mu, phi, Algorithm, base_var_mean_bias) %>% 
  mutate(
    relative_bias = base_var_mean_bias   # it is actually already calculated
  )


#---------------------------------------------------------
# 2 — Summarize by condition (Algorithm × mu × phi)
#    We need the mean relative bias over replicates
#---------------------------------------------------------
df_bias <- df_pois %>%
  group_by(Algorithm, mu, phi) %>%
  summarise(
    rel_bias_mean = mean(relative_bias, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_bias,
       aes(x = mu,
           y = rel_bias_mean,
           color = factor(phi),
           group = phi)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  # dashed ±5% reference lines
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "autocorrelation",
    values = c(
      "0"   = "#2E1A47",
      "0.2" = "#1BA098",
      "0.4" = "#F2D849"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_bias$mu))) +
  
  labs(
    x = "Mean Level",
    y = "Relative Variance Bias"
  ) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)




#gamma point, with the added missing dispersion = 1
#####


#---------------------------------------------------------
# 1 — Keep relevant gamma results
#---------------------------------------------------------
df_gamma <- gamma_point_process_results %>%
  select(mu_a, shape, alg, base_var_mean_bias) %>%
  mutate(
    relative_bias = base_var_mean_bias,
    shape_chr = case_when(
      shape == 2/3 ~ "2/3",
      TRUE ~ as.character(shape)
    )
  )

#---------------------------------------------------------
# 2 — Summarize bias by alg × mean × dispersion
#---------------------------------------------------------
df_gamma_bias <- df_gamma %>%
  group_by(alg, mu_a, shape_chr) %>%
  summarise(
    rel_bias_mean = mean(relative_bias, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_gamma_bias,
       aes(x = mu_a,
           y = rel_bias_mean,
           color = shape_chr,
           group = shape_chr)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "dispersion (shape)",
    values = c(
      "0.4" = "#2E1A47",
      "2/3" = "#1BA098",
      "1"   = "red",
      "1.5" = "#F2D849",
      "2.5" = "#E07A5F"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_gamma_bias$mu_a))) +
  
  labs(
    x = "Mean Level",
    y = "Relative Variance Bias"
  ) +
  
  facet_wrap(~ alg, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)





#normal distribution
#####


#---------------------------------------------------------
# 1 — Keep relevant normal error results
#---------------------------------------------------------
df_normal <- normal_errors_results %>%
  select(phi, sigma_squared, Algorithm, base_var_mean_bias) %>%
  mutate(
    relative_bias = base_var_mean_bias,
    sigma_chr = as.character(sigma_squared)  # Convert variance to character for color mapping
  )

#---------------------------------------------------------
# 2 — Summarize bias by Algorithm × autocorrelation × variance  DIT IS EIGENLIJK OVERBODIG?
#---------------------------------------------------------
df_normal_bias <- df_normal %>%
  group_by(Algorithm, phi, sigma_chr) %>%
  summarise(
    rel_bias_mean = relative_bias,#  mean(, na.rm = TRUE)
    .groups = "drop"
  )

#---------------------------------------------------------
# 3 — Plot
#---------------------------------------------------------
ggplot(df_normal_bias,
       aes(x = phi,
           y = rel_bias_mean,
           color = sigma_chr,
           group = sigma_chr)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  
  scale_color_manual(
    name = "Variance (sigma^2)",
    values = c(
      "0.25" = "#2E1A47",
      "1"   = "#1BA098",
      "2.25"   = "#F2D849"
    )
  ) +
  
  scale_x_continuous(breaks = sort(unique(df_normal_bias$phi))) +
  
  labs(
    x = "Autocorrelation (phi)",
    y = "Relative Variance Bias"
  ) +
  
  facet_wrap(~ Algorithm, ncol = 3, labeller = alg_title_labeller) +
  
  theme_minimal(base_size = 14)


#the missing plot (gamma point, varying mean level)
#####
#
df_gamma <- gamma_point_process_results %>%
  select(mu_a, shape, alg, data) %>%
  unnest(data)    # brings replicate, base_mean, base_var, len

df_gamma_disp15 <- df_gamma %>%
  filter(shape == 2/3, mu_a %in% c(5, 15, 25))

# Total reps per alg × mu_a
totals <- df_gamma_disp15 %>%
  group_by(alg, mu_a) %>%
  summarise(total_reps = n(), .groups = "drop")

# Counts for len ≤ 20
counts <- df_gamma_disp15 %>%
  filter(len <= 20) %>%
  count(alg, mu_a, len, name = "n")

# Full grid + cumulative percentage
df_gamma_cum <- expand_grid(
  alg = unique(totals$alg),
  mu_a = unique(totals$mu_a),
  len = 1:20
) %>%
  left_join(totals, by = c("alg", "mu_a")) %>%
  left_join(counts, by = c("alg", "mu_a", "len")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(alg, mu_a) %>%
  arrange(len) %>%
  mutate(cum_pct = cumsum(n) / total_reps * 100) %>%
  ungroup()

ggplot(df_gamma_cum,
       aes(x = factor(len),
           y = cum_pct,
           fill = factor(mu_a))) +
  
  geom_col(position = "dodge") +
  
  labs(
    x = "Baseline Length",
    y = "Cumulative Percentage",
    fill = "Mean level (mu_a)"
  ) +
  
  scale_x_discrete(breaks = as.character(seq(1, 20, by = 2))) +
  
  facet_wrap(~ alg, ncol = 3, labeller = alg_title_labeller) +
  theme_minimal(base_size = 14)