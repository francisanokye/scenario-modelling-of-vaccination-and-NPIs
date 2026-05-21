suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "National Rate"

# ------------------------------------------------------------
# 1. Build national-rate weekly cumulative curve
# ------------------------------------------------------------

scenario_weekly_df_national <- can_weekly_nl_inputs |>
  arrange(week_end) |>
  filter(week_end >= init_date, week_end <= analysis_end) |>
  mutate(
    weekly_rate  = replace_na(weekly_rate, 0),
    weekly_rate  = pmax(weekly_rate, 0),
    cum_val      = cumsum(weekly_rate),
    cum_val      = pmin(cummax(cum_val), eligible_pop),
    cum_prop     = clip01(cum_val / eligible_pop),
    is_star      = FALSE,
    scenario     = scenario_name
  ) |>
  select(
    scenario,
    week_end,
    weekly_rate,
    cum_prop,
    cum_val,
    is_star
  )

# ------------------------------------------------------------
# 2. Convert weekly national curve to one daily curve
# ------------------------------------------------------------

scenario_daily_df_national <- scenario_weekly_df_national |>
  arrange(week_end) |>
  mutate(
    daydiff = as.integer(
      week_end - lag(week_end, default = week_end[1] - days(7))
    ),
    interval_start = week_end - days(daydiff - 1)
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
    daily_rate = weekly_rate / n_days_in_interval
  ) |>
  ungroup() |>
  arrange(date) |>
  mutate(
    daily_rate = pmax(replace_na(daily_rate, 0), 0),
    cum_val = cumsum(daily_rate),
    cum_val = pmin(cummax(cum_val), eligible_pop),
    cum_prop = clip01(cum_val / eligible_pop),
    days_since_init = as.integer(date - init_date)
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

# ------------------------------------------------------------
# 3. Daily model input directly from daily_rate
# ------------------------------------------------------------

scenario_booster_inputs_df_national <- scenario_daily_df_national |>
  transmute(
    scenario = scenario_name,
    date,
    day = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  filter(day >= 0, date <= analysis_end) |>
  group_by(scenario, day) |>
  summarise(
    date = min(date),
    booster_vac_value = sum(booster_vac_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(day)

# ------------------------------------------------------------
# 4. Validation
# ------------------------------------------------------------

intended_total <- scenario_daily_df_national |>
  filter(date <= analysis_end) |>
  summarise(total = max(cum_val, na.rm = TRUE)) |>
  pull(total)

model_input_total <- scenario_booster_inputs_df_national |>
  summarise(total = sum(booster_vac_value, na.rm = TRUE)) |>
  pull(total)

stopifnot(!any(is.na(scenario_booster_inputs_df_national$booster_vac_value)))
stopifnot(!any(scenario_booster_inputs_df_national$day < 0))
stopifnot(abs(model_input_total - intended_total) < 1e-6)

print(
  tibble(
    scenario = scenario_name,
    first_booster_date = min(
      scenario_booster_inputs_df_national$date[
        scenario_booster_inputs_df_national$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    last_booster_date = max(
      scenario_booster_inputs_df_national$date[
        scenario_booster_inputs_df_national$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    intended_total_doses = intended_total,
    model_input_total_doses = model_input_total,
    final_coverage = intended_total / eligible_pop,
    max_daily_doses = max(
      scenario_booster_inputs_df_national$booster_vac_value,
      na.rm = TRUE
    )
  )
)

saveVars(
  scenario_daily_df_national,
  scenario_weekly_df_national,
  scenario_booster_inputs_df_national
)