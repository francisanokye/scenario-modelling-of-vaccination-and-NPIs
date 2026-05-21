suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Slower Rate"

target_final <- daily_base |>
  filter(date <= analysis_end) |>
  summarise(target = max(cum_value, na.rm = TRUE)) |>
  pull(target)

scenario_daily_df_slower <- make_rollout_scenario(
  daily_base,
  speed = 1 / 1.5,
  shift_days = 0,
  label = scenario_name
) |>
  filter(date >= init_date, date <= analysis_end) |>
  arrange(date) |>
  mutate(
    cum_val = pmax(cum_val, 0),
    cum_val = cummax(cum_val)
  )

current_final <- max(scenario_daily_df_slower$cum_val, na.rm = TRUE)

scale_factor <- if (current_final > 0) {
  target_final / current_final
} else {
  1
}

scenario_daily_df_slower <- scenario_daily_df_slower |>
  mutate(
    cum_val = cum_val * scale_factor,
    cum_val = pmin(cummax(cum_val), eligible_pop),
    daily_rate = c(cum_val[1], diff(cum_val)),
    daily_rate = pmax(daily_rate, 0),
    cum_val = cumsum(daily_rate),
    cum_prop = clip01(cum_val / eligible_pop),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1),
    scenario = scenario_name
  )

scenario_weekly_df_slower <- scenario_daily_df_slower |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(is_star = FALSE) |>
  arrange(week_end)

scenario_booster_inputs_df_slower <- make_booster_input(
  scenario_daily_df_slower,
  scenario_name
)

print(
  validate_rollout(
    scenario_daily_df_slower,
    scenario_booster_inputs_df_slower,
    scenario_name
  )
)

saveVars(
  scenario_daily_df_slower,
  scenario_weekly_df_slower,
  scenario_booster_inputs_df_slower
)