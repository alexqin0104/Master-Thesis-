library(SingleCaseES)
library(future)
library(furrr)

calc_effects <- function(A_data, B_data, improvement){
  smd_ES <- tryCatch(calc_ES(A_data, B_data, ES = "SMD", format = "wide", confidence = NULL), error = function(e) NA)
  if(is.na(smd_ES[1])) smd_ES <- data.frame(SMD_Est = NA, SMD_SE = NA)
  other_es <- calc_ES(A_data, B_data, ES = c("LRRi", "NAP"), format = "wide", confidence = NULL, improvement = improvement)
  bind_cols(other_es, smd_ES)
}

#the begin_at parameter with default value 5; strange
min_length <- function(alg, outcomes, phase_length, begin_at = 5, increase = TRUE){
  stable <- FALSE
  for(i in begin_at:phase_length){
    stable <- do.call(alg, list(outcomes = outcomes[1:i], increase = increase[1]))
    if(stable) break
  }
  i
}

evaluate_by_row <- function(params, sim_function, ...,
                            .progress = FALSE, .options = furrr_options(),
                            system_time = TRUE) {
  
  sys_tm <- system.time(
    results_list <- furrr::future_pmap(
      params,
      .f = sim_function,
      ...,
      .progress = .progress,
      .options = .options
    )
  )
  
  if (system_time) print(sys_tm, "\n")
  
  params %>%
    dplyr::mutate(.results = results_list) %>%
    tidyr::unnest(.data$.results)
  
}