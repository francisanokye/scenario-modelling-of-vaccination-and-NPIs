library(tidyverse)
library(shellpipes)

loadEnvironments()

# number of days before the first data point to start simulations
offset0 <- 150
booster_offset <- 107
NL_POP <- 510550

day0 <- as.Date("2021-12-15") - lubridate::days(offset0)

dat <- csvRead()

# ------------------------------------------------------------
# NL vaccination data
# ------------------------------------------------------------

nl <- (dat 
	|> filter(grepl("ewfoundland", prename, ignore.case = TRUE)) 
	|> select(week_end, contains("numtotal")) 
	|> mutate(NULL
		, numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional)
			, 0 , numtotal_2nd_additional
		)
		, numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  	) 
	|> pivot_longer(-week_end, names_to = "type", values_to = "value")
)

# ------------------------------------------------------------
# Initial two-dose values at warm-up start
# ------------------------------------------------------------

V2_init <- (nl 
	|> filter(type == "numtotal_fully", week_end == as.Date("2021-07-24")) 
	|> pull(value)
)

# ------------------------------------------------------------
# NL actual vaccination inputs
## MLi: Where are the documentations for this?
# ------------------------------------------------------------

nl_vax <- (nl 
	|> filter(type %in% c("numtotal_additional", "numtotal_fully")) 
	|> filter(!is.na(value)) 
	|> arrange(type, week_end) 
	|> group_by(type) 
	|> mutate(NULL
		, daydiff = as.numeric(diff(c(week_end, 0)))
		, boosterdiff = diff(c(value, 0))
		, daily_rate = boosterdiff / daydiff
		, days = as.numeric(week_end - day0)
  	) 
	|> ungroup() 
	|> select(type, days, daily_rate, week_end, value) 
#	|> filter(between(days, 0, 159 + offset0))
)

print(nl_vax,n=Inf)

saveEnvironment()
