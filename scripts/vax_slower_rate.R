suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Slower Rate"

scenario_daily_df_slower <- make_rollout_scenario(
  daily_base,
  speed = 1 / 1.5,
  shift_days = 0,
  label = scenario_name
)

scenario_weekly_df_slower <- scenario_daily_df_slower |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = max(weekly_rate, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(is_star = FALSE)

scenario_booster_inputs_df_slower <- scenario_weekly_df_slower |>
  arrange(week_end) |>
  mutate(
    daydiff = as.numeric(week_end - lag(week_end, default = week_end[1] - 7)),
    interval_start = week_end - (daydiff - 1)
  ) |>
  tidyr::uncount(weights = daydiff, .remove = FALSE) |>
  group_by(week_end) |>
  mutate(
    date = seq(min(interval_start), max(week_end), by = "day"),
    booster_vac_value = weekly_rate / first(daydiff)
  ) |>
  ungroup() |>
  transmute(
    scenario = scenario_name,
    date,
    day = as.integer(date - init_date),
    booster_vac_value = as.numeric(booster_vac_value)
  ) |>
  arrange(day)

saveVars(scenario_daily_df_slower, scenario_weekly_df_slower, scenario_booster_inputs_df_slower)
