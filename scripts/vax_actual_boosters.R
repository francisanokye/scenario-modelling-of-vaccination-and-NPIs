# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(lubridate)
#   library(shellpipes)
# })
# 
# loadEnvironments()
# 
# scenario_name <- "Actual Boosters"
# speed_factor  <- 1.0
# 
# # build the actual booster rollout
# daily_scenarios_df <- bind_rows(
#   make_rollout_scenario(
#     daily_base,
#     speed      = speed_factor,
#     shift_days = 0,
#     label      = scenario_name
#   )
# ) |>
#   mutate(
#     days_since_init = as.integer(date - init_date)
#   ) |>
#   filter(days_since_init >= 0)
# 
# # weekly curve for plotting
# scenario_weekly_df_actual <- daily_scenarios_df |>
#   group_by(scenario, week_end) |>
#   summarise(
#     weekly_rate = sum(daily_rate, na.rm = TRUE),
#     cum_val     = max(cum_val, na.rm = TRUE),
#     cum_prop    = max(cum_prop, na.rm = TRUE),
#     .groups     = "drop"
#   ) |>
#   mutate(is_star = FALSE) |>
#   arrange(week_end) |>
#   filter(week_end >= actual_booster_start + 1)
# 
# # daily booster input for the model
# scenario_booster_inputs_df_actual <- daily_scenarios_df |>
#   transmute(
#     scenario,
#     date,
#     day               = as.integer(date - init_date),
#     booster_vac_value = as.numeric(daily_rate)
#   ) |>
#   arrange(scenario, day) |>
#   filter(date >= actual_booster_start + 1)
# 
# saveVars(
#   scenario_weekly_df_actual,
#   scenario_booster_inputs_df_actual
# )






library(tidyverse)
library(lubridate)
library(shellpipes)

loadEnvironments()

analysis_start <- as.Date("2021-11-08")
analysis_end   <- as.Date("2022-05-22")

init_date        <- analysis_start - days(offset)


actual_boosters <- boosterframe |>
  arrange(date) |>
  filter(date >= analysis_start, date <= analysis_end) |>
  mutate(
    day = as.numeric(date - init_date),
    scenario = "Actual Boosters",
    daily_rate = daily_rate,
    cumulative_total = cumsum(daily_rate)
  ) |>
  select(date, day, scenario, daily_rate, cumulative_total)


rdsSave(actual_boosters)

saveEnvironment()