suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Early Boosters"

# move the actual booster curve earlier
advance_days <- as.integer(actual_booster_start - early_start_date)

actual_cum_fun <- approxfun(
  x      = as.numeric(daily_base$date),
  y      = daily_base$cum_value,
  method = "linear",
  rule   = 2,
  ties   = "ordered"
)

early_cum_raw <- actual_cum_fun(as.numeric(daily_base$date) + advance_days)
early_cum_raw <- ifelse(daily_base$date < early_start_date, 0, early_cum_raw)
early_cum_raw <- pmax(early_cum_raw, 0)
early_cum_raw <- pmin(early_cum_raw, eligible_pop_NL)
early_cum_raw <- cummax(early_cum_raw)

# daily early-booster curve
daily_base_early <- daily_base |>
  mutate(
    cum_value = early_cum_raw,
    daily_inc = pmax(cum_value - lag(cum_value, default = 0), 0),
    week_end  = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(week_end) |>
  mutate(
    weekly_rate = sum(daily_inc, na.rm = TRUE)
  ) |>
  ungroup()

daily_scenarios_df <- daily_base_early |>
  mutate(
    cum_val         = cum_value,
    daily_rate      = daily_inc,
    cum_prop        = clip01(cum_value / NL_POP),
    scenario        = scenario_name,
    days_since_init = as.integer(date - init_date)
  ) |>
  filter(days_since_init >= 0)

# weekly curve for plotting
scenario_weekly_df_early <- daily_scenarios_df |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(is_star = FALSE) |>
  arrange(week_end) |>
  filter(week_end >= early_start_date)

# daily booster input for the model
scenario_booster_inputs_df_early <- daily_scenarios_df |>
  transmute(
    scenario,
    date,
    day               = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  arrange(scenario, day) |>
  filter(day >= 0, date >= early_start_date)

# check final coverage
final_cum <- max(cumsum(scenario_booster_inputs_df_early$booster_vac_value))

cat(sprintf("\n== Early Boosters validation ==\n"))
cat(sprintf("  advance_days:       %d\n", advance_days))
cat(sprintf(
  "  eligible_pop_NL:    %d  (%.1f%%)\n",
  eligible_pop_NL,
  eligible_pop_NL / NL_POP * 100
))
cat(sprintf(
  "  Final coverage:     %.0f  (%.1f%%)\n",
  final_cum,
  final_cum / NL_POP * 100
))

saveVars(
  scenario_weekly_df_early,
  scenario_booster_inputs_df_early
)