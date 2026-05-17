library(tidyverse)

#kazdin 10 percent algo
kaz10 <- function(outcomes,increase = TRUE){
  if(sum(abs(outcomes)) == 0) return(TRUE)
  n <- length(outcomes)
  level <- mean(outcomes)
  
  last_3 <- outcomes[(n-2):n]/level
  
  all(last_3 <= 1.05 & last_3 >= .95)
}

# kazdin 15 percent algo
kaz15 <- function(outcomes,increase = TRUE){
  if(sum(abs(outcomes)) == 0) return(TRUE)
  n <- length(outcomes)
  level <- mean(outcomes)
  
  last_3 <- outcomes[(n-2):n]/level
  
  all(last_3 <= 1.075 & last_3 >= .925)
}

# VVA algo from Ferron et al.
vva<- function(outcomes,increase = TRUE){
  
  n <- length(outcomes)
  stdev <- sd(outcomes)
  
  observations <- 1:n
  
  full_lm <- lm(outcomes ~ observations)
  
  slope <- coef(full_lm)[[2]] 
  
  lm3 <- lm(outcomes[(n-2):n] ~ observations[(n-2):n])
  
  slope3 <- coef(lm3)[[2]]
  
  end_diff <- abs(mean(outcomes) - outcomes[n])
  
  halves_diff <- abs(mean(outcomes[1:floor(n/2)]) - mean(outcomes[ceiling(n/2 + 1):n]))
  
  
  if(abs(slope) <= 0.5 * stdev & abs(slope3) <= 0.5 * stdev & end_diff <=  2*stdev & halves_diff <= 1.5 * stdev){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#stability envelope calculation
stability_envelope <- function(outcomes){
  
  n <- length(outcomes)
  m_b<- median(outcomes)
  
  if(n < 5){
    n_b <- 5
  } else if(n > 20){
    n_b <- 20
  }else{
    n_b <- n
  }
  
  (.30 - n_b / 100) * m_b
}

#stability according to a split-middle envelope
split_middle_stability <- function(outcomes, w, increase, min_stable = 5){
  trend_stable <- FALSE
  
  n <- length(outcomes)
  
  if(n < min_stable) return(trend_stable)
  
  min_n <- ceiling(.8 * n)
  
  n_exclude <- n - min_n
  
  first_half <- floor(n/2)
  second_half <- ceiling(n/2 + 1)
  
  first_date <- median(1:first_half)
  second_date<- median(second_half:n)
  
  first_rate <- median(outcomes[1:first_half])
  second_rate <- median(outcomes[second_half:n])
  
  slope <- (second_rate - first_rate)/(second_date-first_date)
  
  detrend <- outcomes - (0:(n-1))*slope
  
  detrend <- sort(detrend)
  
  trend_dists <- vector(mode = "double")
  for(j in 1:(n_exclude+1)){
    trend_dists[j] <- detrend[(min_n + j-1)] - detrend[j]  
  }
  
  if(!increase) slope <- slope * -1
  
  if(min(trend_dists) <= w & slope <= 0){
    trend_stable <- TRUE
  }
  return(trend_stable)
}

#GL full stability
gl_full <- function(outcomes, increase = TRUE){
  level_stable <- FALSE
  trend_stable <- FALSE
  
  n <- length(outcomes)
  
  w <- stability_envelope(outcomes)
  
  min_n <- ceiling(.8 * n)
  
  sort_outcomes <- sort(outcomes)
  
  n_exclude <- n - min_n
  
  dists <- vector(mode = "double")
  for(i in 1:(n_exclude + 1)){
    dists[i] <- sort_outcomes[(min_n + i-1)] - sort_outcomes[i]  
  }
  
  if(min(dists) <= w){
    level_stable <- TRUE
  }
  
  trend_stable <- split_middle_stability(outcomes, w, increase)
  
  if(trend_stable | level_stable){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#stability according to the last 3 observations  
gl_final <- function(outcomes, increase = TRUE){
  
  level_stable <- FALSE
  trend_stable <- FALSE
  
  n <- length(outcomes)
  
  w <- stability_envelope(outcomes)
  
  if(max(outcomes[(n-2):n]) - min(outcomes[(n-2):n]) <= w) level_stable = TRUE
  
  trend_stable <- split_middle_stability(outcomes, w, increase)
  
  if(trend_stable | level_stable){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#stability according to ends
gl_abs <- function(outcomes, increase = TRUE){
  level_stable <- FALSE
  trend_stable <- FALSE
  
  n <- length(outcomes)
  
  w <- stability_envelope(outcomes)
  
  if(abs(outcomes[1] - outcomes[n]) <= 2 * w) level_stable <- TRUE
  
  trend_stable <- split_middle_stability(outcomes, w, increase)
  
  if(trend_stable | level_stable){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#stability according to median halves
gl_rel <- function(outcomes, increase = TRUE){
  level_stable <- FALSE
  trend_stable <- FALSE
  
  n <- length(outcomes)
  
  w <- stability_envelope(outcomes)
  
  if(abs(median(outcomes[1:floor(n/2)]) - median(outcomes[ceiling(n/2 + 1):n])) <= w) level_stable <- TRUE
  
  trend_stable <- split_middle_stability(outcomes, w, increase)
  
  if(trend_stable | level_stable){
    return(TRUE)
  }else{
    return(FALSE)
  }
}