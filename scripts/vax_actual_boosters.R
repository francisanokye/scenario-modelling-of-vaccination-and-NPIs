suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
})

loadEnvironments()

scenario_name <- "Actual Boosters"
target_date   <- as.Date("2022-01-11")

# ------------------------------------------------------------
# 1. Create actual daily rollout curve
# ------------------------------------------------------------

scenario_daily_df_actual <- make_rollout_scenario(
  daily_base,
  speed = 1.00,
  shift_days = 0,
  label = scenario_name
) |>
  filter(date >= init_date, date <= analysis_end) |>
  arrange(date)

# ------------------------------------------------------------
# 2. Weekly observed/imputed curve for plotting only
# ------------------------------------------------------------

scenario_weekly_df_actual <- nl_actual |>
  transmute(
    scenario    = scenario_name,
    week_end,
    weekly_rate = weekly_inc_star,
    cum_prop,
    cum_val     = corrected_value,
    is_target   = week_end == target_date,
    is_interp   = week_end <= cutoff_date & week_end != target_date,
    is_phac     = week_end > cutoff_date
  ) |>
  filter(week_end >= init_date, week_end <= analysis_end) |>
  arrange(week_end)

# ------------------------------------------------------------
# 3. Daily model input directly from daily_rate
# ------------------------------------------------------------

scenario_booster_inputs_df_actual <- make_booster_input(
  scenario_daily_df_actual,
  scenario_name
)

# ------------------------------------------------------------
# 4. Validation
# ------------------------------------------------------------

print(
  validate_rollout(
    scenario_daily_df_actual,
    scenario_booster_inputs_df_actual,
    scenario_name
  )
)

saveVars(
  scenario_daily_df_actual,
  scenario_weekly_df_actual,
  scenario_booster_inputs_df_actual
)