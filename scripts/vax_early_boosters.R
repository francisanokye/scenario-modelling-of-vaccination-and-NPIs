library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

days_earlier     <- 35
# slow_down_factor <- 0.5

booster_daily_vac <- (boosterframe 
	|> arrange(date) 
	|> mutate(day = as.numeric(date - booster_start))
)

max_historical_day <- max(booster_daily_vac$day)

plot_times_early <- seq(-days_earlier, max_historical_day - days_earlier, by = 1)

# plot_times_main <- seq(0, ceiling(max_historical_day / slow_down_factor),by = 1)

# print(plot_times_main)

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

print(early_func)

print(booster_daily_vac)

new_booster <- (booster_daily_vac
	|> mutate(NULL
  		, day2 = as.numeric(date - burn_in_start)
		, day = day - days_earlier
		, date = date + day
		, day = day2 - days_earlier
	)
	|> select(-day2) 
)

print(new_booster)

## Early boosters: starts before actual rollout
#####################################################################

# early_rates <- early_func(plot_times_early)
# early_rates[is.na(early_rates)] <- 0
# early_rates <- pmax(early_rates, 0)

#results_list[["Early Boosters"]] <- (tibble(
#  date = start_date + plot_times_early
#  , scenario = "Early Boosters"
#  , daily_rate = early_rates
#  )
#  )

#vax_early_rate <- (
#  bind_rows(results_list) 
#  |> mutate(day = as.numeric(date - burn_in_start))
#)

saveEnvironment()
