## What are we trying to process here?
## we need to fill the calibration offset period for boosters with the following features: 
## Booster doses started Nov 8
## we have a point estimate of ~100k in Jan 11th 
## The PHAC vaccination data starts ~ Feb


library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

#####################################################################
## Process booster vaccinated doses
#####################################################################

filldat <- (data.frame(date = c(as.Date("2021-11-08"),as.Date("2022-01-11"))
	, value = c(1,round(0.209 * 510550))
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

boosterframe$cum_interp <- round(na.approx(boosterframe$value))

boosterframe$daily_rate <- c(boosterframe$cum_interp[1], diff(boosterframe$cum_interp))

# boosterframe <- (boosterframe 
#                  |> filter(date >= as.Date("2021-11-08") & 
#                            date <= as.Date("2022-05-22")
#                            )
#                  )

#####################################################################
## Process double vaccinated doses
#####################################################################
## We need to create a complete daily cumulative double-vaccination series
## from the weekly PHAC vaccination data.
## Missing days between weekly observations are filled using
## linear interpolation of the cumulative counts.

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

mindate <- min(doublevaxdat$date)
maxdate <- max(doublevaxdat$date)

doublevaxframe <- (tibble(date = seq.Date(mindate, maxdate, by = "day")
                          , type = "numtotal_fully"
                          ) 
                   |> left_join(doublevaxdat, by = c("date", "type")) 
                   |> mutate(cum_interp = round(na.approx(value, na.rm = FALSE)))
)

doublevaxframe$daily_rate <- c(doublevaxframe$cum_interp[1], diff(doublevaxframe$cum_interp))

# doublevaxframe <- (doublevaxframe 
#                    |> filter(date >= as.Date("2021-07-24") & 
#                              date <= as.Date("2022-05-22")
#                              )
#                    )


saveEnvironment()