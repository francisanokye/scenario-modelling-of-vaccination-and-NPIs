library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

vscale_factor <- 1.5

new_booster <- (boosterframe 
	|> arrange(date) 
	|> mutate(NULL
		, day = as.numeric(date - burn_in_start)
		, daily_rate = daily_rate * vscale_factor 
		, new_cum = cumsum(daily_rate)
		, daily_rate = ifelse(new_cum > booster_max, 0, daily_rate)
	)
	|> filter(daily_rate > 0)
)

## To make sure we have the same boosterbox
lastdat <- data.frame(date = c(as.Date("2022-02-01"),as.Date("2022-02-02"))
	, daily_rate = c(1320,0)
	, day = c(198,199)
)

new_booster <- (new_booster
	|> bind_rows(lastdat)
)

print(new_booster)


#start_date <- min(booster_daily_vac$date)
#max_historical_day <- max(booster_daily_vac$day)

#speed_up_factor  <- 1.5
#slow_down_factor <- 0.5

#plot_times_main <- seq(0, ceiling(max_historical_day / speed_up_factor),by = 1)

#make_rate_fun <- function(day, rate) {
#  approxfun(
#    x = day,
#    y = rate,
#    method = "linear",
#    rule = 1,
#    yleft = 0,
#    yright = 0
#  )
#}

#fast_func <- make_rate_fun(booster_daily_vac$day / speed_up_factor, booster_daily_vac$daily_rate * speed_up_factor)
#
#results_list <- list()

#####################################################################
## Faster boosters: starts same date, finishes earlier
#####################################################################

#faster_rates <- fast_func(plot_times_main)
#faster_rates[is.na(faster_rates)] <- 0
#faster_rates <- pmax(faster_rates, 0)

#results_list[["Faster Boosters"]] <- (tibble(
#  date = start_date + plot_times_main
#  , scenario = "Faster Boosters"
#  , daily_rate = faster_rates
#  )
#)

#vax_faster_rate <- (bind_rows(results_list) 
#  |> mutate(day = as.numeric(date - burn_in_start))
#  )

## replace the first value 0.5 with 1
#vax_faster_rate$daily_rate[1] <- 1

saveEnvironment()
