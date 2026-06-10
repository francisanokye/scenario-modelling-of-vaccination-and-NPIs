suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Faster Boosters"
speed_factor  <- 1.5

# build the faster booster rollout
daily_scenarios_df <- make_rollout_scenario(
  daily_base,
  speed          = speed_factor,
  shift_days     = 0,
  label          = scenario_name,
  scenario_start = actual_booster_start
) |>
  mutate(
    days_since_init = as.integer(date - init_date)
  ) |>
  filter(days_since_init >= 0)

# weekly curve for plotting
scenario_weekly_df_faster <- daily_scenarios_df |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(is_star = FALSE) |>
  arrange(week_end) |>
  filter(week_end >= actual_booster_start + 1)

# daily booster input for the model
scenario_booster_inputs_df_faster <- daily_scenarios_df |>
  transmute(
    scenario,
    date,
    day               = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  arrange(scenario, day) |>
  filter(date >= actual_booster_start + 1)

cat(sprintf(
  "\n== Faster Boosters: final coverage = %.1f%% ==\n",
  max(cumsum(scenario_booster_inputs_df_faster$booster_vac_value)) / NL_POP * 100
))

saveVars(
  scenario_weekly_df_faster,
  scenario_booster_inputs_df_faster
)