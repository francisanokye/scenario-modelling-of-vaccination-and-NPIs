library(tidyverse)
library(shellpipes)

loadEnvironments()

# number of days before the first data point to start simulations
offset0 <- 150

dat <- csvRead()
day0 <- as.Date("2021-12-17") - lubridate::days(offset0)

nl <- (dat
	|> filter(grepl("ewfoundland",prename)) 
	|> select(week_end,contains("numtotal"))
	|> mutate( numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional),0,numtotal_2nd_additional)
		, additional_tot = numtotal_1additional + numtotal_2nd_additional
	)
	|> pivot_longer(-week_end,names_to="type",values_to="value")
)

gg <- (ggplot(nl,aes(week_end,value))
	+ geom_point()
	+ facet_wrap(~type,scale="free")
)

#print(gg)

V2_init <- (nl
	|> filter(type == "numtotal_fully")
	|> filter(week_end == as.Date("2021-07-24"))
	|> pull(value)
)

print(V2_init)

nl_additional <- (nl
	|> filter(type %in% c("numtotal_additional","numtotal_fully"))
	|> mutate(value = ifelse((week_end == as.Date("2021-12-18")) & (type == "numtotal_additional") ,10,value))
	|> filter(!is.na(value))
	|> arrange(type)
	|> group_by(type)
	|> mutate(daydiff = diff(c(week_end,0))
		, boosterdiff = diff(c(value,0))
		, daily_rate = boosterdiff/as.numeric(daydiff)
		, days = as.numeric(week_end - day0)
		, week_end = week_end
	)
	|> select(days,daily_rate,week_end)
	|> filter(between(days,0,159 + offset0))
#	|> filter(between(days,-20,159))

)

#print(nl_additional,n=Inf)

saveVars(nl_additional,V2_init, offset0)

