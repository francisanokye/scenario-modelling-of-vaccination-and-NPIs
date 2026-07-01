library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

vscale_factor <- 0.5

new_booster <- (boosterframe 
	|> arrange(date) 
	|> mutate(NULL
		, day = as.numeric(date - burn_in_start)
		, daily_rate = daily_rate * vscale_factor 
		, new_cum = cumsum(daily_rate)
		, daily_rate = ifelse(new_cum > booster_max, 0, daily_rate)
	)
)

print(new_booster)
saveEnvironment()
