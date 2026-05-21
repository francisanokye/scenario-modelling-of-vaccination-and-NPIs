suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(scales)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(grid)
  library(patchwork)
})

loadEnvironments()

# ============================================================
# Constants
# ============================================================

NL_POP       <- 510550
CA_POP       <- 36991981
eligible_pop <- 426198

init_date      <- as.Date("2021-07-24")
booster_start  <- as.Date("2021-11-08")
backlog_end    <- as.Date("2022-01-11")
cutoff_date    <- as.Date("2022-02-07")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
study_end      <- as.Date("2022-05-22")

earlier_vax_date <- as.Date("2021-10-04")
booster_offset   <- 150

target_date       <- as.Date("2022-01-11")
target_cumulative <- 106730
shape_power       <- 1.8

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

# ============================================================
# Load raw data
# ============================================================

dat <- csvRead()

# ============================================================
# Raw NL booster cumulative counts
# ============================================================

df_raw <- dat |>
  filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = if_else(
      is.na(numtotal_2nd_additional),
      0,
      numtotal_2nd_additional
    ),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  select(week_end, numtotal_additional) |>
  rename(value = numtotal_additional) |>
  filter(
    week_end >= earlier_vax_date,
    week_end <= as.Date("2022-06-19")
  ) |>
  arrange(week_end)

# ============================================================
# Correct cumulative booster series
# ============================================================

first_obs <- df_raw |>
  filter(!is.na(value)) |>
  arrange(week_end) |>
  slice(1)

first_obs_date  <- first_obs$week_end
first_obs_value <- first_obs$value

dates_a <- seq.Date(booster_start, target_date, by = "week")
t_idx   <- seq(0, 1, length.out = length(dates_a))

seg_a <- tibble(
  week_end = dates_a,
  corrected_value = target_cumulative * (t_idx ^ shape_power),
  point_type = "interpolated"
)

if (first_obs_date > target_date) {
  dates_b <- seq.Date(target_date, first_obs_date, by = "week")
  
  seg_b <- tibble(
    week_end = dates_b,
    corrected_value = seq(
      target_cumulative,
      first_obs_value,
      length.out = length(dates_b)
    ),
    point_type = "interpolated"
  ) |>
    filter(week_end > target_date)
  
} else {
  
  seg_b <- tibble(
    week_end = as.Date(character()),
    corrected_value = numeric(),
    point_type = character()
  )
}

df_imputed <- bind_rows(seg_a, seg_b) |>
  filter(week_end < first_obs_date) |>
  distinct(week_end, .keep_all = TRUE) |>
  arrange(week_end)

df_true <- df_raw |>
  filter(week_end >= first_obs_date) |>
  transmute(
    week_end,
    corrected_value = value,
    point_type = "observed"
  )

# ============================================================
# Redistribute artificial post-flat jumps
# ============================================================

redistribute_postflat_jumps <- function(df, backlog_end) {
  
  df <- df |>
    arrange(week_end)
  
  n <- nrow(df)
  i <- 2
  
  while (i <= n) {
    
    if (
      df$week_end[i] > backlog_end &&
      df$corrected_value[i] == df$corrected_value[i - 1]
    ) {
      
      plateau_start <- i - 1
      j <- i
      
      while (
        j <= n &&
        df$corrected_value[j] == df$corrected_value[plateau_start]
      ) {
        j <- j + 1
      }
      
      if (j <= n) {
        
        start_val <- df$corrected_value[plateau_start]
        end_val   <- df$corrected_value[j]
        
        if (end_val > start_val) {
          idx <- plateau_start:j
          df$corrected_value[idx] <- seq(
            start_val,
            end_val,
            length.out = length(idx)
          )
          df$point_type[idx[idx < j]] <- "interpolated"
        }
        
        i <- j
        
      } else {
        
        i <- n + 1
      }
      
    } else {
      
      i <- i + 1
    }
  }
  
  df
}

# ============================================================
# Process observed booster data
# ============================================================

nl_actual <- bind_rows(df_imputed, df_true) |>
  arrange(week_end) |>
  redistribute_postflat_jumps(backlog_end = backlog_end) |>
  mutate(
    corrected_value = pmax(corrected_value, 0),
    corrected_value = pmin(corrected_value, eligible_pop),
    corrected_value = cummax(corrected_value),
    cum_prop        = corrected_value / eligible_pop,
    weekly_inc_star = corrected_value - lag(corrected_value, default = 0),
    weekly_inc_star = pmax(weekly_inc_star, 0),
    daily_inc_star  = weekly_inc_star / 7
  )

# ============================================================
# Build baseline daily cumulative rollout
# ============================================================

daily_dates <- seq.Date(init_date, study_end, by = "1 day")

fit_func <- splinefun(
  x = nl_actual$week_end,
  y = nl_actual$corrected_value,
  method = "hyman"
)

last_obs_date <- max(nl_actual$week_end)
last_obs_cum <- nl_actual |>
  filter(week_end == last_obs_date) |>
  pull(corrected_value)

daily_base <- tibble(date = daily_dates) |>
  mutate(
    cum_value_raw = fit_func(date),
    cum_value_raw = pmax(cum_value_raw, 0),
    cum_value = case_when(
      date < booster_start ~ 0,
      date <= last_obs_date ~ cum_value_raw,
      TRUE ~ last_obs_cum
    ),
    cum_value = pmin(cum_value, eligible_pop),
    cum_value = cummax(cum_value),
    daily_inc = c(cum_value[1], diff(cum_value)),
    daily_inc = pmax(daily_inc, 0),
    cum_value = cumsum(daily_inc),
    cum_prop = clip01(cum_value / eligible_pop),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(week_end) |>
  mutate(
    weekly_rate = sum(daily_inc, na.rm = TRUE)
  ) |>
  ungroup()

# ============================================================
# Time-warp cumulative rollout
# ============================================================

make_rollout_scenario <- function(df, speed = 1, shift_days = 0, label) {
  
  df <- df |>
    arrange(date) |>
    mutate(
      date = as.Date(date),
      cum_value = pmax(cum_value, 0),
      cum_value = cummax(cum_value)
    )
  
  base_dates <- df$date
  base_num   <- as.numeric(base_dates)
  base_cum   <- df$cum_value
  
  cum_fun <- approxfun(
    x = base_num,
    y = base_cum,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  t0 <- min(base_num)
  
  scen_dates <- seq.Date(init_date, study_end, by = "1 day")
  scen_num   <- as.numeric(scen_dates)
  
  mapped_t <- t0 + speed * ((scen_num - shift_days) - t0)
  
  scen_cum <- cum_fun(mapped_t)
  scen_cum <- pmax(scen_cum, 0)
  scen_cum <- pmin(scen_cum, eligible_pop)
  scen_cum <- cummax(scen_cum)
  
  scenario_daily <- tibble(
    date = scen_dates,
    cum_val = scen_cum
  ) |>
    mutate(
      daily_rate = c(cum_val[1], diff(cum_val)),
      daily_rate = pmax(daily_rate, 0),
      cum_val = cumsum(daily_rate),
      cum_prop = clip01(cum_val / eligible_pop),
      week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1),
      scenario = label
    ) |>
    group_by(scenario, week_end) |>
    mutate(
      weekly_rate = sum(daily_rate, na.rm = TRUE)
    ) |>
    ungroup()
  
  scenario_daily
}

# ============================================================
# Helper: convert any scenario daily curve to model input
# ============================================================

make_booster_input <- function(scenario_daily_df, scenario_name) {
  
  scenario_daily_df |>
    arrange(date) |>
    transmute(
      scenario = scenario_name,
      date = as.Date(date),
      day = as.integer(date - init_date),
      booster_vac_value = as.numeric(daily_rate)
    ) |>
    filter(
      day >= 0,
      date <= analysis_end
    ) |>
    group_by(scenario, day) |>
    summarise(
      date = min(date),
      booster_vac_value = sum(booster_vac_value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(day)
}

# ============================================================
# Validation helper
# ============================================================

validate_rollout <- function(scenario_daily_df, booster_input_df, scenario_name) {
  
  intended_total <- scenario_daily_df |>
    filter(date <= analysis_end) |>
    summarise(total = max(cum_val, na.rm = TRUE)) |>
    pull(total)
  
  input_total <- booster_input_df |>
    summarise(total = sum(booster_vac_value, na.rm = TRUE)) |>
    pull(total)
  
  stopifnot(!any(is.na(booster_input_df$booster_vac_value)))
  stopifnot(!any(booster_input_df$day < 0))
  stopifnot(input_total >= 0)
  stopifnot(abs(input_total - intended_total) < 1e-5)
  
  tibble(
    scenario = scenario_name,
    first_booster_date = min(
      booster_input_df$date[booster_input_df$booster_vac_value > 0],
      na.rm = TRUE
    ),
    last_booster_date = max(
      booster_input_df$date[booster_input_df$booster_vac_value > 0],
      na.rm = TRUE
    ),
    intended_total_doses = intended_total,
    model_input_total_doses = input_total,
    final_coverage = intended_total / eligible_pop,
    max_daily_doses = max(booster_input_df$booster_vac_value, na.rm = TRUE)
  )
}

# ============================================================
# Save shared setup
# ============================================================

saveVars(
  NL_POP,
  CA_POP,
  eligible_pop,
  init_date,
  booster_start,
  backlog_end,
  cutoff_date,
  analysis_start,
  analysis_end,
  study_end,
  earlier_vax_date,
  booster_offset,
  target_date,
  target_cumulative,
  shape_power,
  clip01,
  nl_actual,
  daily_base,
  make_rollout_scenario,
  make_booster_input,
  validate_rollout
)