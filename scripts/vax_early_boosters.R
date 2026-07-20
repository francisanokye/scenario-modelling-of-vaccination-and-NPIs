library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

days_earlier     <- 35

vscale_factor <- 1

lastday <- data.frame(date = as.Date("2022-05-23"), daily_rate = 0)

actual_booster <- (boosterframe 
  |> arrange(date) 
  |> bind_rows(lastday)
  |> mutate(NULL
    , day = as.numeric(date - burn_in_start)
    , daily_rate = daily_rate * vscale_factor 
    , new_cum = cumsum(daily_rate)
    , daily_rate = ifelse(new_cum > booster_max, 0, daily_rate)
  )
)

new_booster <- (actual_booster
  |> mutate(NULL
    , date = date - days_earlier
    , day = as.numeric(date - burn_in_start)
    , new_cum = cumsum(daily_rate)
    # , daily_rate = ifelse(new_cum > booster_max, 0, daily_rate)
  )
)

saveEnvironment()
