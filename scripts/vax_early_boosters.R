suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

shift_days <- -35
scenario_name <- "Early Boosters"
scenario_start <- booster_start + days(shift_days)

target_final <- daily_base |>
  filter(date <= analysis_end) |>
  summarise(target = max(cum_value, na.rm = TRUE)) |>
  pull(target)

scenario_daily_df_early <- make_rollout_scenario(
  daily_base,
  speed = 1.00,
  shift_days = shift_days,
  label = scenario_name,
  scenario_start = scenario_start
) |>
  mutate(date = as.Date(date)) |>
  filter(date >= scenario_start, date <= analysis_end) |>
  arrange(date) |>
  mutate(
    cum_val = pmax(cum_val, 0),
    cum_val = cummax(cum_val)
  )

current_final <- max(scenario_daily_df_early$cum_val, na.rm = TRUE)

scale_factor <- if (is.finite(current_final) && current_final > 0) {
  target_final / current_final
} else {
  1
}

scenario_daily_df_early <- scenario_daily_df_early |>
  mutate(
    cum_val = cum_val * scale_factor,
    cum_val = pmin(cummax(cum_val), eligible_pop),
    daily_rate = cum_val - lag(cum_val, default = 0),
    daily_rate = pmax(daily_rate, 0),
    cum_val = cumsum(daily_rate),
    cum_prop = clip01(cum_val / NL_POP),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1),
    scenario = scenario_name
  )

scenario_weekly_df_early <- scenario_daily_df_early |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    cum_val = max(cum_val, na.rm = TRUE),
    cum_prop = max(cum_prop, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(is_star = FALSE) |>
  arrange(week_end)

scenario_booster_inputs_df_early <- make_booster_input(
  scenario_daily_df_early,
  scenario_name
)

print(
  validate_rollout(
    scenario_daily_df_early,
    scenario_booster_inputs_df_early,
    scenario_name
  )
)

first_positive_date <- scenario_booster_inputs_df_early |>
  filter(booster_vac_value > 0) |>
  summarise(first_date = min(date)) |>
  pull(first_date)

print(
  tibble(
    scenario = scenario_name,
    shift_days = shift_days,
    actual_booster_start = booster_start,
    scenario_start = scenario_start,
    first_positive_input_date = first_positive_date,
    start_day = as.integer(scenario_start - init_date),
    first_positive_day = as.integer(first_positive_date - init_date),
    analysis_start_day = as.integer(analysis_start - init_date)
  )
)

stopifnot(
  scenario_start == booster_start + days(shift_days),
  first_positive_date >= scenario_start,
  first_positive_date <= scenario_start + days(14)
)

saveVars(
  scenario_daily_df_early,
  scenario_weekly_df_early,
  scenario_booster_inputs_df_early
)