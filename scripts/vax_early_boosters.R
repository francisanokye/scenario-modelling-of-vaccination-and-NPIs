library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

offset <- 150

booster_start  <- as.Date("2021-11-08")
burn_in_start  <- omicron_start - days(offset)

booster_daily_vac <- (
  boosterframe 
  |> arrange(date) 
  |> mutate(day = as.numeric(date - booster_start))
)

start_date <- min(booster_daily_vac$date)
max_historical_day <- max(booster_daily_vac$day)

days_earlier     <- 35
slow_down_factor <- 0.5

plot_times_early <- seq(-days_earlier, max_historical_day - days_earlier, by = 1)

plot_times_main <- seq(0, ceiling(max_historical_day / slow_down_factor),by = 1)

make_rate_fun <- function(day, rate) {
  approxfun(
    x = day,
    y = rate,
    method = "linear",
    rule = 1,
    yleft = 0,
    yright = 0
  )
}

early_func <- make_rate_fun(booster_daily_vac$day - days_earlier, booster_daily_vac$daily_rate)

results_list <- list()

#####################################################################
## Early boosters: starts before actual rollout
#####################################################################

early_rates <- early_func(plot_times_early)
early_rates[is.na(early_rates)] <- 0
early_rates <- pmax(early_rates, 0)

results_list[["Early Boosters"]] <- (tibble(
  date = start_date + plot_times_early
  , scenario = "Early Boosters"
  , daily_rate = early_rates
  )
  )

vax_early_rate <- (
  bind_rows(results_list) 
  |> mutate(day = as.numeric(date - burn_in_start))
)

saveEnvironment()
