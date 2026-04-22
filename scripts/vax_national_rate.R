suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(tidyr)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "National Rate"

scenario_weekly_df_national <- can_weekly_nl_inputs |>
  arrange(week_end) |>
  filter(week_end >= booster_start, week_end <= study_end) |>
  mutate(
    weekly_rate   = replace_na(weekly_rate, 0),
    weekly_rate   = pmax(weekly_rate, 0),
    cum_boosters  = cumsum(weekly_rate),
    cum_boosters  = cummax(cum_boosters),
    cum_boosters  = pmin(cum_boosters, eligible_pop),
    cum_prop      = clip01(cum_boosters / eligible_pop),
    cum_val       = cum_boosters,
    is_star       = FALSE,
    scenario      = scenario_name
  ) |>
  select(scenario, week_end, weekly_rate, cum_prop, cum_val, is_star)

scenario_daily_df_national <- scenario_weekly_df_national |>
  arrange(week_end) |>
  mutate(
    daydiff = as.numeric(week_end - lag(week_end, default = week_end[1] - 7)),
    interval_start = week_end - (daydiff - 1)
  ) |>
  tidyr::uncount(weights = daydiff, .remove = FALSE) |>
  group_by(week_end) |>
  mutate(
    date = seq(min(interval_start), max(week_end), by = "day"),
    daily_rate = weekly_rate / first(daydiff),
    days_since_init = as.integer(date - init_date)
  ) |>
  ungroup() |>
  arrange(date) |>
  mutate(
    cum_val = cumsum(daily_rate),
    cum_val = pmin(cummax(cum_val), eligible_pop),
    cum_prop = clip01(cum_val / eligible_pop)
  ) |>
  transmute(
    date,
    cum_val,
    daily_rate,
    cum_prop,
    week_end,
    scenario,
    weekly_rate,
    days_since_init
  )

scenario_booster_inputs_df_national <- scenario_weekly_df_national |>
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

saveVars(scenario_daily_df_national, scenario_weekly_df_national, scenario_booster_inputs_df_national)
