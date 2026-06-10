suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
  library(ggplot2)
})

loadEnvironments()

scenario_name <- "National Rate"

# replace one large weekly outlier and spread the excess forward
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

# weekly canadian-rate curve capped at the canadian final uptake level
scenario_weekly_df_national <- can_weekly_nl_inputs |>
  arrange(week_end) |>
  filter(
    week_end >= init_date,
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
    cum_val  = cumsum(weekly_rate),
    cum_val  = pmin(cummax(cum_val), eligible_pop_CA),
    cum_prop = clip01(cum_val / NL_POP),
    is_star  = FALSE,
    scenario = scenario_name
  ) |>
  select(
    scenario,
    week_end,
    weekly_rate,
    cum_prop,
    cum_val,
    is_star,
    outlier_corrected
  ) |>
  mutate(
    weekly_rate = pmax(cum_val - lag(cum_val, default = 0), 0)
  )

print(
  scenario_weekly_df_national |>
    filter(outlier_corrected) |>
    select(week_end, weekly_rate, cum_val)
)

# spread each weekly total evenly across its reporting interval
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
    date >= init_date,
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
    daily_rate      = pmax(replace_na(daily_rate, 0), 0),
    cum_val         = cumsum(daily_rate),
    cum_prop        = clip01(cum_val / NL_POP),
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
    days_since_init,
    outlier_corrected
  )

# daily booster input for the model
scenario_booster_inputs_df_national <- scenario_daily_df_national |>
  transmute(
    scenario          = scenario_name,
    date,
    day               = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  filter(
    day >= 0,
    date <= analysis_end
  ) |>
  group_by(scenario, day) |>
  summarise(
    date              = min(date),
    booster_vac_value = sum(booster_vac_value, na.rm = TRUE),
    .groups           = "drop"
  ) |>
  arrange(day)

# checks
intended_total <- max(scenario_daily_df_national$cum_val, na.rm = TRUE)

model_input_total <- sum(
  scenario_booster_inputs_df_national$booster_vac_value,
  na.rm = TRUE
)

stopifnot(!any(is.na(scenario_booster_inputs_df_national$booster_vac_value)))
stopifnot(!any(scenario_booster_inputs_df_national$day < 0))
stopifnot(abs(model_input_total - intended_total) < 1e-6)

print(
  tibble(
    scenario                = scenario_name,
    first_booster_date      = min(
      scenario_booster_inputs_df_national$date[
        scenario_booster_inputs_df_national$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    last_booster_date       = max(
      scenario_booster_inputs_df_national$date[
        scenario_booster_inputs_df_national$booster_vac_value > 0
      ],
      na.rm = TRUE
    ),
    intended_total_doses    = round(intended_total),
    model_input_total_doses = round(model_input_total),
    final_coverage_pct      = round(intended_total / NL_POP * 100, 2),
    target_coverage_pct     = 46.70,
    max_daily_doses         = round(max(
      scenario_booster_inputs_df_national$booster_vac_value,
      na.rm = TRUE
    ))
  )
)

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
        "Target final coverage: 46.7%% of NL_POP (%d people)",
        eligible_pop_CA
      )
    ) +
    theme_minimal()
)

saveVars(
  scenario_daily_df_national,
  scenario_weekly_df_national,
  scenario_booster_inputs_df_national
)