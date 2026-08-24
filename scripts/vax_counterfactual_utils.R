
make_rollout_scenario <- function(
    boosterframe,
    burn_in_start,
    booster_max,
    rollout_speed = 1,
    shift_days = 0
) {
  if (
    length(rollout_speed) != 1
    || !is.finite(rollout_speed)
    || rollout_speed <= 0
  ) {
    stop("rollout_speed must be a positive number.")
  }
  
  observed <- (boosterframe
    |> transmute(
     date = as.Date(date),
     daily_rate = as.numeric(daily_rate)
    )
    |> filter(
     !is.na(date),
     is.finite(daily_rate),
     daily_rate >= 0
    )
    |> arrange(date)
    |> group_by(date)
    |> summarise(
     daily_rate = sum(daily_rate, na.rm = TRUE),
     .groups = "drop"
    )
    |> mutate(
     observed_cum = cumsum(daily_rate)
    )
  )
  
  observed_total <- sum(observed$daily_rate)
  
  target_total <- min(as.numeric(booster_max), observed_total)
  
  observed_time <- seq.int(0, nrow(observed) - 1)
  
  new_duration <- ceiling((nrow(observed) - 1) / rollout_speed)
  
  new_time <- seq.int(0, new_duration)
  
  mapped_time <- pmin(new_time * rollout_speed, max(observed_time))
  
  new_cum <- approx(
    x = observed_time,
    y = observed$observed_cum,
    xout = mapped_time,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
  
  new_cum <- pmin(cummax(new_cum), target_total)
  
  daily_rate <- c(new_cum[1], diff(new_cum))
  
  daily_rate <- pmax(daily_rate, 0)
  
  cumulative_before <- c(0, head(cumsum(daily_rate), -1))
  
  daily_rate <- pmin(daily_rate, pmax(target_total - cumulative_before, 0))
  
  residual <- target_total - sum(daily_rate)
  
  if (residual > 1e-8) {
    positive_index <- which(daily_rate > 0)
    
    final_index <- max(positive_index)
    
    daily_rate[final_index] <- (daily_rate[final_index] + residual)
  }
  
  positive_index <- which(daily_rate > 0)
  
  final_positive_index <- max(positive_index)
  
  daily_rate <- daily_rate[
    seq_len(final_positive_index)
  ]
  
  start_date <- (min(observed$date) - lubridate::days(shift_days))
  
  new_booster <- tibble(
    date = seq.Date(
      start_date,
      by = "day",
      length.out = length(daily_rate)
    ),
    daily_rate = daily_rate
  )
  
  terminal_zero_date <- (
    max(new_booster$date) +
      lubridate::days(1)
  )
  
  new_booster <- bind_rows(new_booster,tibble(date = terminal_zero_date,daily_rate = 0)) |>
    arrange(date) |>
    mutate(
      day = as.numeric(date - burn_in_start),
      new_cum = cumsum(daily_rate)
    )
  
  final_total <- sum(new_booster$daily_rate)

  
  new_booster
}