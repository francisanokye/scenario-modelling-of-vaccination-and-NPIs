## What are we trying to process here?
## we need to fill the calibration offset period for boosters with the following features: 
## Booster doses started Nov 8
## we have a point estimate of ~100k (20.90% of NL pop) in Jan 11th 
## The PHAC booster vaccination data starts ~ Feb 13, 2022


library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

#####################################################################
## Load PHAC vaccination data for NL
#####################################################################

dat <- csvRead()

nl <- (dat 
	|> filter(grepl("ewfoundland", prename, ignore.case = TRUE)) 
	|> select(week_end, contains("numtotal")) 
	|> pivot_longer(-week_end, names_to = "type", values_to = "value")
)

# ------------------------------------------------------------
# Initial two-dose values at warm-up start
# ------------------------------------------------------------
# Pulled at the PHAC week_end nearest burn_in_start (2021-07-17, one
# day before the actual burn-in start, burn_in_start = 2021-07-18).
V2_init <- (nl 
	|> filter(type == "numtotal_fully", week_end == as.Date("2021-07-17")) 
	|> pull(value)
)


nl_vax <- (nl 
	|> filter(type %in% c("numtotal_additional", "numtotal_fully")) 
	|> filter(!is.na(value)) 
	|> arrange(type, week_end) 
	|> select(type, week_end, value) 
)

#####################################################################
## Process booster vaccinated doses
#####################################################################

filldat <- (data.frame(date = c(as.Date("2021-11-08"),as.Date("2022-01-11"))
	, value = c(1,round(0.209 * NL_POP))
	, type = "numtotal_additional"
	)
)

boosterdat <- (nl_vax
	|> filter(type == "numtotal_additional")
	|> transmute(NULL
	, type
	, date = week_end
	, value
	)
)

boosterframe <- (data.frame(date = as.Date(c(min(filldat$date):max(boosterdat$date))))
	|> mutate(type = "numtotal_additional")
	|> left_join(bind_rows(filldat,boosterdat))
)

boosterframe <- (boosterframe 
	|> mutate(cum_interp = round(na.approx(value))
		, daily_rate = pmax(cum_interp - lag(cum_interp, default = 0), 0)
	)
	|> filter(date >= booster_start & date <= analysis_end)
	|> select(date, cum_interp, daily_rate)
)

#####################################################################
## Process double vaccinated doses
#####################################################################

## We need to create a complete daily cumulative double-vaccination series
## from the weekly PHAC vaccination data.
## Missing days between weekly observations are filled using
## linear interpolation of the cumulative counts.
## Unlike boosters, real double-dose data exists back through the
## burn-in period, so this series is NOT truncated to booster_start --
## it extends back to whatever real PHAC double-dose data is available,
## down to burn_in_start.

# double-dose weekly cumulative data
doublevaxdat <- (nl_vax 
  |> filter(type == "numtotal_fully") 
  |> transmute(NULL
    , type
    , date = week_end
    , value
    ) 
  |> arrange(date) 
  |> mutate(value = cummax(value))
)

doublevaxframe <- (tibble(date = seq.Date(min(doublevaxdat$date), max(doublevaxdat$date), by = "day"), type = "numtotal_fully") 
  |> left_join(doublevaxdat, by = c("date", "type"))
  |> mutate(cum_interp = round(na.approx(value))
      , cum_interp = cummax(cum_interp)
      , daily_rate = pmax(cum_interp - lag(cum_interp, default = first(cum_interp)), 0)
      ) 
  |> filter(date >= burn_in_start, date <= analysis_end)
  |> select(date, cum_interp, daily_rate)
)

# print(doublevaxframe)

saveEnvironment()
