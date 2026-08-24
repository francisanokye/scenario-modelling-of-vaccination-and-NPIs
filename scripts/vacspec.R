
library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

calib <- rdsRead("calibrate")
vac_scenario <- pipeStar()

# ------------------------------------------------------------
# Booster vaccination trajectory
# ------------------------------------------------------------

new_changepoints <- as.integer(new_booster$day)
new_values <- as.numeric(new_booster$daily_rate)

if (!any(new_changepoints == 0)) {
  new_changepoints <- c(0L, new_changepoints)
  new_values <- c(0, new_values)
}

newspec <- mp_optimized_spec(calib, "modified")

newspec$integers$booster_vac_changepoints <- as.integer(new_changepoints)

newspec$default$booster_vac_values <- as.numeric(new_values)

# ------------------------------------------------------------
# Primary-series burn-in counterfactual
# ------------------------------------------------------------

if (exists("primary_series_target")) {
  
  target_day <- as.integer(primary_series_target_date - burn_in_start)
  
  # ----------------------------------------------------------
  # Initial V2 state
  # ----------------------------------------------------------
  
  default_V2_init <- as.numeric(newspec$default$V2_init)
  
  if (exists("use_new_V2_init") && isTRUE(use_new_V2_init)) {
    
    scenario_V2_init <- as.numeric(new_V2_init)
    
    newspec$default$V2_init <- scenario_V2_init
    
  } else {
    
    scenario_V2_init <- default_V2_init
  }
  
  # ----------------------------------------------------------
  # Original calibrated two-dose vaccination trajectory
  # ----------------------------------------------------------
  
  default_double_changepoints <- as.integer(newspec$integers$double_vac_changepoints)
  
  default_double_values <- as.numeric(newspec$default$double_vac_values)
  
  double_order <- order(default_double_changepoints)
  
  default_double_changepoints <- default_double_changepoints[double_order]
  
  default_double_values <- default_double_values[double_order]
  
  # ----------------------------------------------------------
  # Expand only the burn-in period to daily values
  # ----------------------------------------------------------
  
  burnin_days <- seq.int(0L, target_day)
  
  burnin_index <- findInterval(
    burnin_days
    , default_double_changepoints
  )
  
  burnin_daily_rate <- ifelse(
    burnin_index > 0
    , default_double_values[pmax(burnin_index, 1)]
    , 0
  )
  
  # ----------------------------------------------------------
  # Force primary-series coverage on 6 November
  # ----------------------------------------------------------
  
  required_burnin_doses <- primary_series_target - scenario_V2_init
  
  original_burnin_doses <- sum(burnin_daily_rate, na.rm = TRUE)
  
  burnin_scale_factor <- required_burnin_doses / original_burnin_doses
  
  # ----------------------------------------------------------
  # Correct for realised S1 -> V2 vaccination flow
  # ----------------------------------------------------------
  
  burnin_correction <- case_when(
    vac_scenario == "vax_national_rate" ~ 1.020840
    , vac_scenario == "nl_init_national_rate" ~ 1.025080
    , vac_scenario == "vax_actual_rate_can_init" ~ 1.004800
    , TRUE ~ 1
  )
  
  burnin_scale_factor <- burnin_scale_factor * burnin_correction
  
  burnin_daily_rate <- burnin_daily_rate * burnin_scale_factor
  
  # ----------------------------------------------------------
  # Scheduled cumulative primary-series total
  # ----------------------------------------------------------
  
  primary_series_target_check <- scenario_V2_init + sum(burnin_daily_rate, na.rm = TRUE)
  
  # ----------------------------------------------------------
  # Restore ORIGINAL double_vac_values after 6 November
  # ----------------------------------------------------------
  
  post_burnin_day <- target_day + 1L
  
  post_burnin_index <- findInterval(
    post_burnin_day
    , default_double_changepoints
  )
  
  post_burnin_value <- if (post_burnin_index > 0) {
    default_double_values[post_burnin_index]
  } else {
    0
  }
  
  post_burnin_default <- tibble(
    day = default_double_changepoints
    , daily_rate = default_double_values
  ) |>
    filter(day > post_burnin_day)
  
  scenario_double_vax <- bind_rows(
    tibble(day = burnin_days, daily_rate = burnin_daily_rate)
    , tibble(day = post_burnin_day, daily_rate = post_burnin_value)
    , post_burnin_default
  ) |>
    arrange(day) |>
    group_by(day) |>
    slice_tail(n = 1) |>
    ungroup()
  
  newspec$integers$double_vac_changepoints <- as.integer(
    scenario_double_vax$day
  )
  
  newspec$default$double_vac_values <- as.numeric(
    scenario_double_vax$daily_rate
  )
  
}

saveVars(newspec, vac_scenario)
