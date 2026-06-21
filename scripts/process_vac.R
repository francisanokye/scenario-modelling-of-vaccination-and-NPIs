## What are we trying to process here?
## we need to fill the calibration offset period for boosters with the following features: 
## Booster doses started Nov 8
## we have a point estimate of ~100k in Jan 11th 
## The PHAC vaccination data starts ~ Feb


library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

print(nl_vax,n=Inf)

filldat <- (data.frame(date = c(as.Date("2021-11-08"),as.Date("2022-01-11"))
	, value = c(1,106380)
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

print(boosterdat)

print(mindate <- min(filldat$date))
print(maxdate <- max(boosterdat$date))

boosterframe <- (data.frame(date = as.Date(c(min(filldat$date):max(boosterdat$date))))
	|> mutate(type = "numtotal_additional")
	|> left_join(bind_rows(filldat,boosterdat))
)


boosterframe$cum_interp <- round(na.approx(boosterframe$value))

print(boosterframe)


