
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
  library(ggplot2)
})

loadEnvironments()

scenario_name <- "National Rate"

day_start <- as.Date("2021-11-08")

# ==============================================================
# Build can_weekly_nl_inputs: the Canada-wide weekly additional-dose
# series from the same raw vaccination-coverage-map.csv used for the NL
# (actual) scenario, rescaled from Canada's population to NL's population.
#
# Rescaling: at each week_end, Canada's cumulative additional-dose count
# is converted to a proportion of CA_POP, then that SAME proportion is
# applied to NL_POP to get an "NL-equivalent" cumulative count 
# ==============================================================
dat <- csvRead()

can_weekly_nl_inputs <- dat |>
  filter(prename == "Canada") |>
  select(week_end, contains("numtotal")) |>
  mutate(
    week_end = as.Date(week_end),
    numtotal_1additional   = ifelse(is.na(numtotal_1additional), 0, numtotal_1additional),
    numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional), 0, numtotal_2nd_additional),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional,
    # rescale Canada's cumulative proportion to NL's population
    cum_val_raw = numtotal_additional / CA_POP * NL_POP
  ) |>
  arrange(week_end) |>
  mutate(weekly_rate = pmax(cum_val_raw - lag(cum_val_raw, default = 0), 0)) |>
  select(week_end, weekly_rate)

# ==============================================================
# Replace one large weekly outlier and spread the excess forward
# (unchanged from the original national scenario script; addresses
# reporting-lag artifacts in the Canada-wide series, e.g. the
# 2022-01-22 spike caused by an 8-day reporting gap)
# ==============================================================
replace_weekly_outlier_preserve_total <- function(
    df,
    value_col            = "weekly_rate",
    date_col             = "week_end",
    outlier_date         = NULL,
    threshold_multiplier = 2.5,
    redistribute_weeks   = 4
) {
  df <- df |>
    arrange(.data[[date_col]]) |>
    mutate(
      row_id    = row_number(),
      value_raw = as.numeric(.data[[value_col]])
    )
  
  if (is.null(outlier_date)) {
    local_ref <- df |>
      mutate(
        prev_value          = lag(value_raw),
        next_value          = lead(value_raw),
        neighbour_mean      = rowMeans(cbind(prev_value, next_value), na.rm = TRUE),
        ratio_to_neighbours = value_raw / neighbour_mean
      )
    
    outlier_id <- local_ref |>
      filter(
        is.finite(ratio_to_neighbours),
        ratio_to_neighbours >= threshold_multiplier
      ) |>
      slice_max(ratio_to_neighbours, n = 1, with_ties = FALSE) |>
      pull(row_id)
    
  } else {
    outlier_id <- df |>
      filter(.data[[date_col]] == as.Date(outlier_date)) |>
      pull(row_id)
  }
  
  if (length(outlier_id) == 0 || is.na(outlier_id[1])) {
    return(
      df |>
        mutate(
          "{value_col}" := value_raw,
          outlier_corrected = FALSE
        ) |>
        select(-row_id, -value_raw)
    )
  }
  
  if (outlier_id <= 1 || outlier_id >= nrow(df)) {
    stop("Outlier cannot be first or last row because neighbours are required.")
  }
  
  replacement <- mean(
    c(df$value_raw[outlier_id - 1], df$value_raw[outlier_id + 1]),
    na.rm = TRUE
  )
  
  excess <- df$value_raw[outlier_id] - replacement
  
  df$value_raw[outlier_id] <- replacement
  df$outlier_corrected <- FALSE
  df$outlier_corrected[outlier_id] <- TRUE
  
  if (is.finite(excess) && excess > 0) {
    redistribute_ids <- seq(
      outlier_id + 1,
      min(outlier_id + redistribute_weeks, nrow(df))
    )
    
    redistribute_ids <- redistribute_ids[
      is.finite(df$value_raw[redistribute_ids])
    ]
    
    if (length(redistribute_ids) > 0) {
      df$value_raw[redistribute_ids] <-
        df$value_raw[redistribute_ids] + excess / length(redistribute_ids)
    }
  }
  
  df |>
    mutate(
      value_raw = pmax(value_raw, 0),
      "{value_col}" := value_raw
    ) |>
    select(-row_id, -value_raw)
}

# ==============================================================
# Weekly national-rate curve, rescaled to NL population.
# ==============================================================
scenario_weekly_df_national <- can_weekly_nl_inputs |>
  arrange(week_end) |>
  filter(
    week_end >= day_start,
    week_end <= analysis_end
  ) |>
  mutate(
    weekly_rate = replace_na(weekly_rate, 0),
    weekly_rate = pmax(weekly_rate, 0)
  ) |>
  replace_weekly_outlier_preserve_total(
    value_col            = "weekly_rate",
    date_col             = "week_end",
    outlier_date         = NULL,
    threshold_multiplier = 2.5,
    redistribute_weeks   = 4
  ) |>
  mutate(
    cum_val  = cummax(pmax(cumsum(weekly_rate), 0)),
    cum_prop = clip01(cum_val / NL_POP),
    scenario = scenario_name
  ) |>
  arrange(week_end) |>
  mutate(weekly_rate = pmax(cum_val - lag(cum_val, default = 0), 0)) |>
  mutate(
    is_target = FALSE,
    is_phac   = FALSE,
    is_interp = TRUE,
    source    = "national_rescaled"
  ) |>
  select(
    scenario, week_end, weekly_rate, cum_val, cum_prop,
    is_target, is_phac, is_interp, source, outlier_corrected
  )

# ==============================================================
# Spread each weekly total evenly across its reporting interval, to get
# a daily-resolution series
# ==============================================================
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
  filter(
    date >= day_start,
    date <= analysis_end
  ) |>
  group_by(scenario, week_end) |>
  mutate(
    n_days_in_interval = n(),
    daily_rate         = weekly_rate / n_days_in_interval
  ) |>
  ungroup() |>
  arrange(date) |>
  mutate(
    daily_rate = pmax(replace_na(daily_rate, 0), 0),
    cum_val    = cummax(pmax(cumsum(daily_rate), 0)),
    cum_prop   = clip01(cum_val / NL_POP),
    day        = as.integer(date - init_date)
  )

# ==============================================================
# Daily booster input for the model (combiner-compatible schema)
# ==============================================================
scenario_booster_inputs_df_national <- scenario_daily_df_national |>
  transmute(
    scenario          = scenario_name,
    date,
    day               = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate),
    cumulative        = cum_val,
    prop_total_pop    = cum_prop
  ) |>
  filter(day >= 0, date <= analysis_end) |>
  group_by(scenario, day) |>
  summarise(
    date              = min(date),
    booster_vac_value = sum(booster_vac_value, na.rm = TRUE),
    cumulative        = max(cumulative, na.rm = TRUE),
    prop_total_pop    = max(prop_total_pop, na.rm = TRUE),
    .groups           = "drop"
  ) |>
  arrange(day) |>
  select(scenario, date, day, booster_vac_value, cumulative, prop_total_pop)



# quick diagnostic plot
print(
  ggplot(scenario_daily_df_national, aes(week_end, daily_rate)) +
    geom_line() +
    geom_point(aes(shape = outlier_corrected), size = 2) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x       = "Week end",
      y       = "Daily booster-dose rate (Canada-scaled to NL)",
      shape   = "Corrected week",
      caption = sprintf(
        "Reference ceiling (not applied as a cap): %.1f%% of NL_POP (%d people)",
        eligible_pop_CA / NL_POP * 100, eligible_pop_CA
      )
    ) +
    theme_minimal()
)

# ==============================================================
saveVars(
  scenario_daily_df_national,
  scenario_weekly_df_national,
  scenario_booster_inputs_df_national
)
