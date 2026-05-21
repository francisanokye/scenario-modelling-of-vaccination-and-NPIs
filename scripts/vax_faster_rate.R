suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Faster Rate"

# ------------------------------------------------------------
# 1. Create faster rollout curve
# ------------------------------------------------------------

scenario_daily_df_faster <- make_rollout_scenario(
  daily_base,
  speed = 1.50,
  shift_days = 0,
  label = scenario_name
) |>
  filter(date >= init_date, date <= analysis_end) |>
  arrange(date)

# ------------------------------------------------------------
# 2. Weekly summary for plotting only
# ------------------------------------------------------------

scenario_weekly_df_faster <- scenario_daily_df_faster |>
  mutate(
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(is_star = FALSE) |>
  arrange(week_end)

# ------------------------------------------------------------
# 3. Daily model input directly from daily_rate
# ------------------------------------------------------------

scenario_booster_inputs_df_faster <- make_booster_input(
  scenario_daily_df_faster,
  scenario_name
)

# ------------------------------------------------------------
# 4. Validation
# ------------------------------------------------------------

print(
  validate_rollout(
    scenario_daily_df_faster,
    scenario_booster_inputs_df_faster,
    scenario_name
  )
)

saveVars(
  scenario_daily_df_faster,
  scenario_weekly_df_faster,
  scenario_booster_inputs_df_faster
)