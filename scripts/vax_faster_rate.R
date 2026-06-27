library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

offset <- 150

booster_start  <- as.Date("2021-11-08")
burn_in_start  <- omicron_start - days(offset)

booster_daily_vac <- boosterframe |>
  arrange(date) |>
  mutate(day = as.numeric(date - booster_start))

start_date <- min(booster_daily_vac$date)
max_historical_day <- max(booster_daily_vac$day)

speed_up_factor  <- 1.5
slow_down_factor <- 0.5

plot_times_main <- seq(0, ceiling(max_historical_day / speed_up_factor),by = 1)

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

fast_func <- make_rate_fun(booster_daily_vac$day / speed_up_factor, booster_daily_vac$daily_rate * speed_up_factor)

results_list <- list()

#####################################################################
## Faster boosters: starts same date, finishes earlier
#####################################################################

faster_rates <- fast_func(plot_times_main)
faster_rates[is.na(faster_rates)] <- 0
faster_rates <- pmax(faster_rates, 0)

results_list[["Faster Boosters"]] <- (tibble(
  date = start_date + plot_times_main
  , scenario = "Faster Boosters"
  , daily_rate = faster_rates
  )
)

vax_faster_rate <- (bind_rows(results_list) 
  |> mutate(day = as.numeric(date - burn_in_start))
  )

# replace the first value 0.5 with 1
vax_faster_rate$daily_rate[1] <- 1

saveEnvironment()
