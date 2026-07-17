
library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

# ------------------------------------------------------------
#  Load PHAC vaccination data
# ------------------------------------------------------------
dat <- csvRead()

Canada <- (dat 
	|> filter(grepl("Canada", prename, ignore.case = TRUE)) 
	|> select(week_end, contains("numtotal")) 
	|> pivot_longer(-week_end, names_to = "type", values_to = "value")
)

# ------------------------------------------------------------
# Initial two-dose values at warm-up start
# ------------------------------------------------------------
new_V2_init <- (Canada 
	|> filter(type == "numtotal_fully", week_end == as.Date("2021-07-17")) 
	|> pull(value)
)

new_V2_init <- round(new_V2_init * (NL_POP/CA_POP))


Canada <- (Canada 
	|> filter(type %in% c("numtotal_additional")) 
	|> arrange(type, week_end) 
	|> mutate(value = round(value * (NL_POP/CA_POP))) 
	|> select(type, week_end, value) 
	|> filter(week_end >= as.Date("2021-11-06"), week_end <= as.Date("2022-05-22") )
)


#####################################################################
## Process booster vaccinated doses
#####################################################################

filldat <- data.frame(date = c(as.Date("2021-11-06"),as.Date("2021-12-11"))
	, value = c(1, 27322)
	, type = "numtotal_additional"
)


can_boosterdat <- (Canada
	|> filter(type == "numtotal_additional")
	|> transmute(NULL
		, type
		, date = week_end
		, value
	)
)

can_boosterframe <- (data.frame(date = as.Date(c(min(filldat$date):max(can_boosterdat$date))))
	|> mutate(type = "numtotal_additional")
	|> left_join(bind_rows(filldat,can_boosterdat))
)


can_boosterframe <- (can_boosterframe 
	|> mutate(cum_interp = round(na.approx(value))
		, daily_rate = pmax(cum_interp - lag(cum_interp, default = 0), 0)
		)
	|> filter(date >= booster_start
		, date <= analysis_end
		, daily_rate > 0
		)
	|> select(date, cum_interp, daily_rate)
)


new_booster<- (can_boosterframe 
	|> mutate(day = as.numeric(date - burn_in_start)
	)
)

# print(new_booster |> ggplot(aes(date,daily_rate)) + geom_line() +geom_point())

saveEnvironment()



