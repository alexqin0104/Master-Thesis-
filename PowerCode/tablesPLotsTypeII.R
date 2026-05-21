library(dplyr)
library(knitr)
library(kableExtra)



env1 <- new.env()
env2 <- new.env()
env3 <- new.env()



#data loading for adaptive data
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIBinomial.Rdata", envir = env1)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIGamma.Rdata", envir = env2)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIINormal.Rdata", envir = env3)
binomial_thinning_results <- env1[[ls(env1)]]
gamma_point_process_results <- env2[[ls(env2)]]
normal_errors_results <- env3[[ls(env3)]]



test = load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIBinomial49Observations.Rdata")

# fixed49 datasets
env4 <- new.env()
env5 <- new.env()
env6 <- new.env()

load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIBinomial49Observations.Rdata", envir = env4)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIGamma49Observations.Rdata", envir = env5)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIINormal49Observations.Rdata", envir = env6)

binomial_fixed49_results <- env4[[ls(env4)]]
gamma_fixed49_results    <- env5[[ls(env5)]]
normal_fixed49_results   <- env6[[ls(env6)]]



# Accidentally called "fixed5" as the name for the "fixed49"
#(this was because I simply changed the fixed5 algorithm instead of making a new algorithm, as a test at first, but
#forgot to then later make a completely separate fixed49 algorithm)
binomial_fixed49_results$Algorithm <- "fixed49"
gamma_fixed49_results$Algorithm    <- "fixed49"
normal_fixed49_results$Algorithm   <- "fixed49"





#appending the fixed49 results
binomial_thinning_results <- bind_rows(
  binomial_thinning_results,
  binomial_fixed49_results
)

gamma_point_process_results <- bind_rows(
  gamma_point_process_results,
  gamma_fixed49_results
)

normal_errors_results <- bind_rows(
  normal_errors_results,
  normal_fixed49_results
)



summarise_metrics <- function(df) {
  df %>%
    summarise(
      power = mean(power),
      beta  = mean(beta),
      
      MCSE_power = sqrt(mean(MCSE_power^2)),
      MCSE_beta  = sqrt(mean(MCSE_beta^2)),
      
      power_lower = mean(power_lower),
      power_upper = mean(power_upper),
      
      beta_lower = mean(beta_lower),
      beta_upper = mean(beta_upper),
      
      .groups = "drop"
    )
}




table_by_algorithm <- function(df) {
  df %>%
    group_by(Algorithm) %>%
    summarise(
      power = mean(power),
      beta  = mean(beta),
      
      MCSE_power = sqrt(mean(MCSE_power^2)),
      MCSE_beta  = sqrt(mean(MCSE_beta^2)),
      
      power_lower = mean(power_lower),
      power_upper = mean(power_upper),
      
      beta_lower = mean(beta_lower),
      beta_upper = mean(beta_upper),
      
      .groups = "drop"
    ) %>%
    arrange(desc(power))
}


table_by_parameters <- function(df, param_cols) {
  df %>%
    filter(!Algorithm %in% c("fixed5", "fixed49")) %>%
    group_by(across(all_of(param_cols))) %>%
    summarise(
      power = mean(power),
      beta  = mean(beta),
      
      MCSE_power = sqrt(mean(MCSE_power^2)),
      MCSE_beta  = sqrt(mean(MCSE_beta^2)),
      
      power_lower = mean(power_lower),
      power_upper = mean(power_upper),
      
      beta_lower = mean(beta_lower),
      beta_upper = mean(beta_upper),
      
      .groups = "drop"
    )
}

table_fixed_only <- function(df, param_cols) {
  df %>%
    filter(Algorithm == "fixed5") %>%
    select(
      all_of(param_cols),
      power, beta,
      MCSE_power, MCSE_beta,
      power_lower, power_upper,
      beta_lower, beta_upper
    )
}



to_latex_table <- function(df, caption, label) {
  df %>%
    select(-MCSE_power, -MCSE_beta) %>%   # remove MCSE columns
    format_table() %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      caption = caption,
      label = label,
      align = "lrrrrrr"
    ) %>%
    kable_styling(latex_options = c("hold_position"))
}



binomial_alg_table <- table_by_algorithm(binomial_thinning_results)

binomial_param_table <- table_by_parameters(
  binomial_thinning_results,
  param_cols = c("phi", "mu")
) %>%
  arrange(phi, mu)

binomial_fixed_table <- table_fixed_only(
  binomial_thinning_results,
  param_cols = c("mu", "phi")
)


gamma_alg_table <- table_by_algorithm(gamma_point_process_results)

gamma_param_table <- table_by_parameters(
  gamma_point_process_results,
  param_cols = c("mu_a", "shape")
)

gamma_fixed_table <- table_fixed_only(
  gamma_point_process_results,
  param_cols = c("mu_a", "shape")
)

normal_alg_table <- table_by_algorithm(normal_errors_results)

normal_param_table <- table_by_parameters(
  normal_errors_results,
  param_cols = c("phi", "sigma_squared")
) %>%
  arrange(phi, sigma_squared)

normal_fixed_table <- table_fixed_only(
  normal_errors_results,
  param_cols = c("phi", "sigma_squared")
)




format_table <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

format_table(binomial_alg_table)













latex_binomial_alg <- to_latex_table(
  binomial_alg_table,
  "Average power and Type II error by algorithm under binomial thinning.",
  "tab:BinomialPowerAlgorithms"
)

latex_binomial_param <- to_latex_table(
  binomial_param_table,
  "Average power and Type II error under binomial thinning by parameter setting.",
  "tab:BinomialPowerParameters"
)

latex_binomial_fixed <- to_latex_table(
  binomial_fixed_table,
  "Reference power and Type II error under binomial thinning (fixed baseline).",
  "tab:BinomialPowerReference"
)



latex_gamma_alg <- to_latex_table(
  gamma_alg_table,
  "Average power and Type II error by algorithm under the gamma point process.",
  "tab:GammaPowerAlgorithms"
)

latex_gamma_param <- to_latex_table(
  gamma_param_table,
  "Average power and Type II error under the gamma point process by parameter setting.",
  "tab:GammaPowerParameters"
)

latex_gamma_fixed <- to_latex_table(
  gamma_fixed_table,
  "Reference power and Type II error under the gamma point process.",
  "tab:GammaPowerReference"
)



latex_normal_alg <- to_latex_table(
  normal_alg_table,
  "Average power and Type II error by algorithm under normally distributed data.",
  "tab:NormalPowerAlgorithms"
)

latex_normal_param <- to_latex_table(
  normal_param_table,
  "Average power and Type II error under normal data by parameter setting.",
  "tab:NormalPowerParameters"
)

latex_normal_fixed <- to_latex_table(
  normal_fixed_table,
  "Reference power and Type II error under normal data.",
  "tab:NormalPowerReference"
)




cat(latex_binomial_alg)
cat(latex_binomial_param)
cat(latex_binomial_fixed)

cat(latex_gamma_alg)
cat(latex_gamma_param)
cat(latex_gamma_fixed)

cat(latex_normal_alg)
cat(latex_normal_param)
cat(latex_normal_fixed)





#plots (is just copy paste from the type I error rate plots)

library(ggplot2)


bin_ref_fixed5 <- binomial_thinning_results %>%
  filter(Algorithm == "fixed5") %>%
  group_by(phi, mu) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )

bin_ref_fixed49 <- binomial_thinning_results %>%
  filter(Algorithm == "fixed49") %>%
  group_by(phi, mu) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )




gamma_ref_fixed5 <- gamma_point_process_results %>%
  filter(Algorithm == "fixed5") %>%
  group_by(shape, mu_a) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = shape - 0.05,
    xmax = shape + 0.05
  )

gamma_ref_fixed49 <- gamma_point_process_results %>%
  filter(Algorithm == "fixed49") %>%
  group_by(shape, mu_a) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = shape - 0.05,
    xmax = shape + 0.05
  )





norm_ref_fixed5 <- normal_errors_results %>%
  filter(Algorithm == "fixed5") %>%
  group_by(phi, sigma_squared) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )

norm_ref_fixed49 <- normal_errors_results %>%
  filter(Algorithm == "fixed49") %>%
  group_by(phi, sigma_squared) %>%
  summarise(
    ref = mean(power),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )









norm_summary <- normal_errors_results %>%
  filter(!Algorithm %in% c("fixed5", "fixed49")) %>%
  group_by(phi, sigma_squared, Algorithm) %>%
  summarise(
    power = mean(power),
    power_lower = mean(power_lower),
    power_upper = mean(power_upper),
    .groups = "drop"
  ) %>%
  mutate(DGM = "Normal")

pd <- position_dodge(width = 0.05)

ggplot(norm_summary, aes(x = phi, y = power, color = Algorithm)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = norm_ref, aes(yintercept = ref), color = "red") +
  
  geom_errorbar(
    aes(ymin = power_lower, ymax = power_upper),
    position = pd,
    width = 0.08,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = Algorithm), position = pd, linewidth = 0.8) +
  facet_wrap(~sigma_squared) +
  labs(
    title = "Power (Normal errors)",
    x = expression(phi),
    y = "Power"
  ) +
  
  # fixed5 reference
  geom_segment(
    data = norm_ref_fixed5,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = norm_ref_fixed5,
    aes(x = phi, y = ref),
    color = "red",
    size = 3,
    shape = 18
  ) +
  
  # fixed49 reference
  geom_segment(
    data = norm_ref_fixed49,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "darkblue",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = norm_ref_fixed49,
    aes(x = phi, y = ref),
    color = "blue",
    size = 3,
    shape = 17
  ) +
  
  

  theme_minimal()








bin_summary <- binomial_thinning_results %>%
  filter(!Algorithm %in% c("fixed5", "fixed49")) %>%
  group_by(phi, mu, Algorithm) %>%
  summarise(
    power = mean(power),
    power_lower = mean(power_lower),
    power_upper = mean(power_upper),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.05)

ggplot(bin_summary, aes(x = phi, y = power, color = Algorithm)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = bin_ref, aes(yintercept = ref), color = "red") +
  geom_errorbar(
    aes(ymin = power_lower, ymax = power_upper),
    position = pd,
    width = 0.08,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = Algorithm), position = pd, linewidth = 0.8) +
  facet_wrap(~mu) +
  labs(
    title = "Power (Binomial thinning)",
    x = expression(phi),
    y = "Power"
  ) +

  
  
  
  
  # fixed5 reference
  geom_segment(
    data = bin_ref_fixed5,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = bin_ref_fixed5,
    aes(x = phi, y = ref),
    color = "red",
    size = 3,
    shape = 18
  ) +
  
  # fixed49 reference
  geom_segment(
    data = bin_ref_fixed49,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "darkblue",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = bin_ref_fixed49,
    aes(x = phi, y = ref),
    color = "blue",
    size = 3,
    shape = 17
  ) +
  
  theme_minimal()



gamma_summary <- gamma_point_process_results %>%
  filter(!Algorithm %in% c("fixed5", "fixed49")) %>%
  group_by(shape, mu_a, Algorithm) %>%
  summarise(
    power = mean(power),
    power_lower = mean(power_lower),
    power_upper = mean(power_upper),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.05)

ggplot(gamma_summary, aes(x = shape, y = power, color = Algorithm)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = gamma_ref, aes(yintercept = ref), color = "red") +
  
  geom_errorbar(
    aes(ymin = power_lower, ymax = power_upper),
    position = pd,
    width = 0.03,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = Algorithm), position = pd, linewidth = 0.8) +
  facet_wrap(~mu_a) +
  labs(
    title = "Power (Gamma point process)",
    x = expression(shape),
    y = "Power"
  ) +

  
  
  
  # fixed5 reference
  geom_segment(
    data = gamma_ref_fixed5,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = gamma_ref_fixed5,
    aes(x = shape, y = ref),
    color = "red",
    size = 3,
    shape = 18
  ) +
  
  # fixed49 reference
  geom_segment(
    data = gamma_ref_fixed49,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "darkblue",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = gamma_ref_fixed49,
    aes(x = shape, y = ref),
    color = "blue",
    size = 3,
    shape = 17
  ) +
  
  
  
  
  
  
  theme_minimal()





norm_summary <- normal_errors_results %>%
  filter(!Algorithm %in% c("fixed5", "fixed49")) %>%
  group_by(phi, sigma_squared, Algorithm) %>%
  summarise(
    power = mean(power),
    power_lower = mean(power_lower),
    power_upper = mean(power_upper),
    .groups = "drop"
  ) %>%
  mutate(DGM = "Normal")

pd <- position_dodge(width = 0.05)




#the plots of baseline length versus power

rm(list = ls())
load("C:/Users/eloua/Desktop/master aangewezen/extension/type I error results 20000/binomthin.Rdata")
load("C:/Users/eloua/Desktop/master aangewezen/extension/type I error results 20000/gpp.Rdata")
load("C:/Users/eloua/Desktop/master aangewezen/extension/type I error results 20000/normCorrected.Rdata")

env1 <- new.env()
env2 <- new.env()
env3 <- new.env()

#data loading for adaptive data
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIBinomial.Rdata", envir = env1)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIIGamma.Rdata", envir = env2)
load("C:/Users/eloua/Desktop/master aangewezen/extension/typeII/TypeIINormal.Rdata", envir = env3)
binomial_thinning_power <- env1[[ls(env1)]]
gamma_point_process_power <- env2[[ls(env2)]]
normal_errors_power <- env3[[ls(env3)]]




#Extracting the mean baseline lengths
library(dplyr)
library(tidyr)




bin_baseline_summary <- binomial_thinning_results %>%
  select(Algorithm, data) %>%
  unnest(data) %>%                 # expands the inner data frames
  group_by(Algorithm) %>%
  summarise(
    mean_baseline_length = mean(len, na.rm = TRUE),
    .groups = "drop"
  )



#fixed data baseline length is always 5, but we still need to add it to the tables.
bin_baseline_summary <- bin_baseline_summary %>%
  bind_rows(
    tibble(
      Algorithm = "fixed5",
      mean_baseline_length = 5
    )
  )


gpp_baseline_summary <- gamma_point_process_results %>%
  select(alg, data) %>%
  unnest(data) %>%                 # expands the inner data frames
  group_by(alg) %>%
  summarise(
    mean_baseline_length = mean(len, na.rm = TRUE),
    .groups = "drop"
  )

gpp_baseline_summary <- gpp_baseline_summary %>%
  bind_rows(
    tibble(
      alg = "fixed5",
      mean_baseline_length = 5
    )
  ) %>%
  rename(Algorithm = alg)


norm_baseline_summary <- normal_errors_results %>%
  select(Algorithm, data) %>%
  unnest(data) %>%                 # expands the inner data frames
  group_by(Algorithm) %>%
  summarise(
    mean_baseline_length = mean(len, na.rm = TRUE),
    .groups = "drop"
  )

norm_baseline_summary <- norm_baseline_summary %>%
  bind_rows(
    tibble(
      Algorithm = "fixed5",
      mean_baseline_length = 5
    )
  )





#scatterplots


#========================
# BINOMIAL THINNING PLOT
#========================

bin_baseline_summary$mean_baseline_length <- 
  round(bin_baseline_summary$mean_baseline_length, 1)

bin_power_summary <- binomial_thinning_power %>%
  group_by(Algorithm) %>%
  summarise(
    power = mean(power),
    .groups = "drop"
  )

bin_plot_data <- left_join(
  bin_power_summary,
  bin_baseline_summary,
  by = "Algorithm"
)


bin_plot_data <- bin_plot_data %>%
  filter(!Algorithm %in% c("fixed5", "fixed49"))



ggplot(
  bin_plot_data,
  aes(x = mean_baseline_length, y = power, color = Algorithm)
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = mean_baseline_length),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  
  # fixed49 reference
  geom_hline(
    yintercept = 0.409,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.409,
    label = "fixed49",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  # fixed5 reference
  geom_hline(
    yintercept = 0.208,
    linetype = "dashed",
    color = "grey40",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.208,
    label = "fixed5",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  labs(
    title = "Power vs Baseline Length (Binomial thinning)",
    x = "Mean baseline length",
    y = "Power"
  ) +
  theme_minimal()


#===========================
# GAMMA POINT PROCESS PLOT
#===========================

gpp_baseline_summary$mean_baseline_length <-
  round(gpp_baseline_summary$mean_baseline_length, 1)

gamma_power_summary <- gamma_point_process_power %>%
  group_by(Algorithm) %>%
  summarise(
    power = mean(power),
    .groups = "drop"
  )

gamma_plot_data <- left_join(
  gamma_power_summary,
  gpp_baseline_summary,
  by = "Algorithm"
)



gamma_plot_data <- gamma_plot_data %>%
  filter(!Algorithm %in% c("fixed5", "fixed49"))


ggplot(
  gamma_plot_data,
  aes(x = mean_baseline_length, y = power, color = Algorithm)
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = mean_baseline_length),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  
  # fixed49 reference
  geom_hline(
    yintercept = 0.952,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.952,
    label = "fixed49",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  # fixed5 reference
  geom_hline(
    yintercept = 0.725,
    linetype = "dashed",
    color = "grey40",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.725,
    label = "fixed5",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  labs(
    title = "Power vs Baseline Length (Gamma point process)",
    x = "Mean baseline length",
    y = "Power"
  ) +
  theme_minimal()



#=====================
# NORMAL ERRORS PLOT
#=====================

norm_baseline_summary$mean_baseline_length <-
  round(norm_baseline_summary$mean_baseline_length, 1)

norm_power_summary <- normal_errors_power %>%
  group_by(Algorithm) %>%
  summarise(
    power = mean(power),
    .groups = "drop"
  )

norm_plot_data <- left_join(
  norm_power_summary,
  norm_baseline_summary,
  by = "Algorithm"
)

norm_plot_data <- norm_plot_data %>%
  filter(!Algorithm %in% c("fixed5", "fixed49"))



ggplot(
  norm_plot_data,
  aes(x = mean_baseline_length, y = power, color = Algorithm)
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = mean_baseline_length),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  
  # fixed49 reference
  geom_hline(
    yintercept = 0.399,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.399,
    label = "fixed49",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  # fixed5 reference
  geom_hline(
    yintercept = 0.198,
    linetype = "dashed",
    color = "grey40",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = 0.198,
    label = "fixed5",
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  labs(
    title = "Power vs Baseline Length (Normal errors)",
    x = "Mean baseline length",
    y = "Power"
  ) +
  theme_minimal()



