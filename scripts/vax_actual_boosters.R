
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Actual Boosters"
target_date   <- as.Date("2022-01-11")

scenario_daily_df_actual <- make_rollout_scenario(
  daily_base,
  speed = 1.00,
  shift_days = 0,
  label = scenario_name,
  scenario_start = booster_start
) |>
  filter(date >= booster_start, date <= analysis_end) |>
  arrange(date)

scenario_weekly_df_actual <- nl_actual |>
  transmute(
    scenario = scenario_name,
    week_end,
    weekly_rate = weekly_inc_star,
    cum_prop,
    cum_val = corrected_value,
    is_target = week_end == target_date,
    is_interp = week_end <= cutoff_date & week_end != target_date,
    is_phac = week_end > cutoff_date
  ) |>
  filter(week_end >= booster_start, week_end <= analysis_end) |>
  arrange(week_end)

scenario_booster_inputs_df_actual <- make_booster_input(
  scenario_daily_df_actual,
  scenario_name
)

print(
  validate_rollout(
    scenario_daily_df_actual,
    scenario_booster_inputs_df_actual,
    scenario_name
  )
)

actual_first_positive_date <- scenario_booster_inputs_df_actual |>
  filter(booster_vac_value > 0) |>
  summarise(first_date = min(date)) |>
  pull(first_date)

print(
  tibble(
    scenario = scenario_name,
    intended_start = booster_start,
    first_positive_input_date = actual_first_positive_date,
    lag_days = as.integer(actual_first_positive_date - booster_start)
  )
)

stopifnot(
  actual_first_positive_date >= booster_start,
  actual_first_positive_date <= booster_start + days(7)
)

saveVars(
  scenario_daily_df_actual,
  scenario_weekly_df_actual,
  scenario_booster_inputs_df_actual
)