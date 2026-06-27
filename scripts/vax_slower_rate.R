
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

slow_down_factor <- 0.5

plot_times_main <- seq(0, ceiling(max_historical_day / slow_down_factor), by = 1)

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

slow_func <- make_rate_fun(booster_daily_vac$day / slow_down_factor, booster_daily_vac$daily_rate * slow_down_factor)

results_list <- list()

#####################################################################
## Slower boosters: starts same date, finishes later
#####################################################################

slower_rates <- slow_func(plot_times_main)
slower_rates[is.na(slower_rates)] <- 0
slower_rates <- pmax(slower_rates, 0)

results_list[["Slower Boosters"]] <- tibble(
  date = start_date + plot_times_main
  , scenario = "Slower Boosters"
  , daily_rate = slower_rates
  )

vax_slower_rate <- (bind_rows(results_list) 
                    |> mutate(day = as.numeric(date - burn_in_start))
                    )
# replace the first value 0.5 with 1
vax_slower_rate$daily_rate[1] <- 1


print(vax_slower_rate)

saveEnvironment()

