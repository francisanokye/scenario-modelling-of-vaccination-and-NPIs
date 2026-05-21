suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Early Boosters (-35 Days)"

# ------------------------------------------------------------
# 1. Create early rollout curve
# ------------------------------------------------------------

scenario_daily_df_early <- make_rollout_scenario(
  daily_base,
  speed = 1.00,
  shift_days = -35,
  label = scenario_name
) |>
  mutate(date = as.Date(date)) |>
  filter(date >= init_date, date <= analysis_end) |>
  arrange(date) |>
  mutate(
    cum_val  = pmax(cum_val, 0),
    cum_val  = cummax(cum_val),
    cum_prop = cum_val / eligible_pop
  )

# ------------------------------------------------------------
# 2. Convert shifted cumulative curve to weekly totals
# ------------------------------------------------------------

scenario_weekly_df_early <- scenario_daily_df_early |>
  mutate(
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(scenario, week_end) |>
  summarise(
    cum_val  = max(cum_val, na.rm = TRUE),
    cum_prop = max(cum_prop, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(week_end) |>
  mutate(
    prev_cum_val = lag(cum_val, default = 0),
    weekly_doses = pmax(cum_val - prev_cum_val, 0),
    weekly_rate = weekly_doses / 7,
    is_star = FALSE
  ) |>
  select(-prev_cum_val)

# ------------------------------------------------------------
# 3. Spread weekly totals evenly across days
# ------------------------------------------------------------

scenario_booster_inputs_df_early <- scenario_weekly_df_early |>
  mutate(
    interval_start = week_end - days(6)
  ) |>
  rowwise() |>
  mutate(
    date = list(seq(interval_start, week_end, by = "day"))
  ) |>
  unnest(date) |>
  ungroup() |>
  filter(date >= init_date, date <= analysis_end) |>
  group_by(scenario, week_end) |>
  mutate(
    n_days_in_interval = n(),
    booster_vac_value = weekly_doses / n_days_in_interval
  ) |>
  ungroup() |>
  mutate(
    day = as.integer(date - init_date)
  ) |>
  filter(day >= 0) |>
  transmute(
    scenario = scenario_name,
    date,
    day,
    booster_vac_value = as.numeric(booster_vac_value)
  ) |>
  arrange(day)

# ------------------------------------------------------------
# 4. Validation checks
# ------------------------------------------------------------

intended_total <- scenario_weekly_df_early |>
  summarise(total = max(cum_val, na.rm = TRUE)) |>
  pull(total)

model_input_total <- scenario_booster_inputs_df_early |>
  summarise(total = sum(booster_vac_value, na.rm = TRUE)) |>
  pull(total)

stopifnot(!any(is.na(scenario_booster_inputs_df_early$booster_vac_value)))
stopifnot(!any(scenario_booster_inputs_df_early$day < 0))
stopifnot(model_input_total > 0)
stopifnot(abs(model_input_total - intended_total) < 1e-6)

print(
  tibble(
    scenario = scenario_name,
    first_booster_date = min(
      scenario_booster_inputs_df_early$date[
        scenario_booster_inputs_df_early$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    last_booster_date = max(
      scenario_booster_inputs_df_early$date[
        scenario_booster_inputs_df_early$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    intended_total_doses = intended_total,
    model_input_total_doses = model_input_total,
    final_coverage = intended_total / eligible_pop,
    max_daily_doses = max(
      scenario_booster_inputs_df_early$booster_vac_value,
      na.rm = TRUE
    )
  )
)

saveVars(
  scenario_daily_df_early,
  scenario_weekly_df_early,
  scenario_booster_inputs_df_early
)