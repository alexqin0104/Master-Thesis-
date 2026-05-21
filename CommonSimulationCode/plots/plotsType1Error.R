#####Reference data loading
#this is needed to make the plots for the fixed5 "algorithm"
#to change the names of each dataset
env1 <- new.env()
env2 <- new.env()
env3 <- new.env()

# load into those environments
#load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/binomthinFixed.Rdata", envir = env1)
#load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/gppFixed.Rdata", envir = env2)
#load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/normCorrectedFixed.Rdata", envir = env3)

# extract objects with new names
#binomial_thinning_results <- env1[[ls(env1)]]
#gamma_point_process_results <- env2[[ls(env2)]]
#normal_errors_results <- env3[[ls(env3)]]

#data loading for adaptive data
load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/binomthinAdaptive.Rdata", envir = env1)
load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/gppAdaptive.Rdata", envir = env2)
load("C:/Users/eloua/Desktop/master aangewezen/extension/fixed_baseline_21_04/normCorrectedAdaptive.Rdata", envir = env3)
binomial_thinning_results <- env1[[ls(env1)]]
gamma_point_process_results <- env2[[ls(env2)]]
normal_errors_results <- env3[[ls(env3)]]



typeIBin = c("mu", "phi", "Algorithm","type1_error", "MCSE", "MC_lower", "MC_upper")#binomial
typeIgam = c("mu_a", "shape", "alg", "type1_error", "MCSE", "MC_lower", "MC_upper")
typeInorm = c("phi", "sigma_squared", "Algorithm", "type1_error", "MCSE", "MC_lower", "MC_upper")

bintypeI = binomial_thinning_results[typeIBin]
gamtypeI = gamma_point_process_results[typeIgam]
normtypeI = normal_errors_results[typeInorm]



#the reference data
bin_ref <- data.frame(
  phi = c(0.0,0.0,0.0, 0.2,0.2,0.2, 0.4,0.4,0.4),
  mu  = c(5,15,25, 5,15,25, 5,15,25),
  ref = c(0.04760,0.04890,0.05145,
          0.08960,0.09395,0.09380,
          0.14985,0.15385,0.15275)
)
bin_ref <- bin_ref %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )



gamma_ref <- data.frame(
  shape = c(
    0.4,0.4,0.4,
    0.6667,0.6667,0.6667,
    1.0,1.0,1.0,
    1.5,1.5,1.5,
    2.5,2.5,2.5
  ),
  mu_a = rep(c(5,15,25), 5),
  ref = c(
    0.04550,0.04700,0.05035,
    0.04395,0.04760,0.04940,
    0.04780,0.05160,0.04845,
    0.05060,0.04945,0.04995,
    0.04945,0.05005,0.04835
  )
)
gamma_ref <- gamma_ref %>%
  mutate(
    xmin = shape - 0.05,
    xmax = shape + 0.05
  )


norm_ref <- data.frame(
  phi = c(0.0,0.0,0.0, 0.2,0.2,0.2, 0.4,0.4,0.4),
  sigma_squared = c(0.25,1.00,2.25, 0.25,1.00,2.25, 0.25,1.00,2.25),
  ref = c(0.04930,0.04765,0.05195,
          0.09790,0.09825,0.09520,
          0.15875,0.15245,0.15415)
)
norm_ref <- norm_ref %>%
  mutate(
    xmin = phi - 0.05,
    xmax = phi + 0.05
  )




#aggregate over the algorithms (too many rows, impossible to include in the report properly)
library(dplyr)
library(xtable)

norm_agg <- normtypeI %>%
  group_by(phi, sigma_squared) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  ) %>%
  arrange(phi, sigma_squared)

print(xtable(norm_agg, digits = 6),
      include.rownames = FALSE)


bin_agg <- bintypeI %>%
  group_by(phi, mu) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  ) %>%
  arrange(phi, mu)

print(xtable(bin_agg, digits = 6),
      include.rownames = FALSE)

gam_agg <- gamtypeI %>%
  group_by(shape, mu_a) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  ) %>%
  arrange(shape, mu_a)

print(xtable(gam_agg, digits = 6),
      include.rownames = FALSE)




#Now aggregate over the paramters, for completeness sake
norm_alg_agg <- normtypeI %>%
  group_by(Algorithm) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  )

bin_alg_agg <- bintypeI %>%
  group_by(Algorithm) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  )

gam_alg_agg <- gamtypeI %>%
  group_by(alg) %>%
  summarise(
    type1_error = mean(type1_error),
    MCSE = mean(MCSE),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  )

print(xtable(norm_alg_agg, digits = 6),
      include.rownames = FALSE)

print(xtable(bin_alg_agg, digits = 6),
      include.rownames = FALSE)

print(xtable(gam_alg_agg, digits = 6),
      include.rownames = FALSE)






#plots
library(ggplot2)

norm_summary <- normtypeI %>%
  group_by(phi, sigma_squared, Algorithm) %>%
  summarise(
    type1_error = mean(type1_error),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  ) %>%
  mutate(DGM = "Normal")

pd <- position_dodge(width = 0.05)

ggplot(norm_summary, aes(x = phi, y = type1_error, color = Algorithm)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = norm_ref, aes(yintercept = ref), color = "red") +
 
  geom_errorbar(
    aes(ymin = MC_lower, ymax = MC_upper),
    position = pd,
    width = 0.08,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = Algorithm), position = pd, linewidth = 0.8) +
  facet_wrap(~sigma_squared) +
  labs(
    title = "Type I Error (Normal errors)",
    x = expression(phi),
    y = "Type I Error"
  ) +
  geom_segment(
    data = norm_ref,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) + 
  geom_point(
    data = norm_ref,
    aes(x = phi, y = ref),
    color = "red",
    size = 3,
    shape = 18  # diamond (very visible)
  ) +
  theme_minimal()



bin_summary <- bintypeI %>%
  group_by(phi, mu, Algorithm) %>%
  summarise(
    type1_error = mean(type1_error),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.05)

ggplot(bin_summary, aes(x = phi, y = type1_error, color = Algorithm)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = bin_ref, aes(yintercept = ref), color = "red") +
  geom_errorbar(
    aes(ymin = MC_lower, ymax = MC_upper),
    position = pd,
    width = 0.08,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = Algorithm), position = pd, linewidth = 0.8) +
  facet_wrap(~mu) +
  labs(
    title = "Type I Error (Binomial thinning)",
    x = expression(phi),
    y = "Type I Error"
  ) +
  geom_segment(
    data = bin_ref,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) +
  geom_point(
    data = bin_ref,
    aes(x = phi, y = ref),
    color = "red",
    size = 3,
    shape = 18  # diamond (very visible)
  ) +
  theme_minimal()



gamma_summary <- gamtypeI %>%
  group_by(shape, mu_a, alg) %>%
  summarise(
    type1_error = mean(type1_error),
    MC_lower = mean(MC_lower),
    MC_upper = mean(MC_upper),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.05)

ggplot(gamma_summary, aes(x = shape, y = type1_error, color = alg)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  #geom_hline(data = gamma_ref, aes(yintercept = ref), color = "red") +
 
  geom_errorbar(
    aes(ymin = MC_lower, ymax = MC_upper),
    position = pd,
    width = 0.03,
    alpha = 0.9
  ) +
  geom_point(position = pd, size = 2.2) +
  geom_line(aes(group = alg), position = pd, linewidth = 0.8) +
  facet_wrap(~mu_a) +
  labs(
    title = "Type I Error (Gamma point process)",
    x = expression(shape),
    y = "Type I Error"
  ) +
  geom_segment(
    data = gamma_ref,
    aes(x = xmin, xend = xmax, y = ref, yend = ref),
    color = "black",
    linewidth = 1.5
  ) +
  geom_point(
    data = gamma_ref,
    aes(x = shape, y = ref),
    color = "red",
    size = 3,
    shape = 18  # diamond (very visible)
  ) +
  theme_minimal()







