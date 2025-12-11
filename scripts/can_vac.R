library(tidyverse)
library(shellpipes)

loadEnvironments()

# number of days before the first data point to start simulations
offset0 <- 150

dat <- csvRead()
day0 <- as.Date("2021-12-17") - lubridate::days(offset0)

can <- (dat
        |> filter(grepl("Canada",prename)) 
        |> select(week_end,contains("numtotal"))
        |> mutate( numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional),0,numtotal_2nd_additional)
                   , additional_tot = numtotal_1additional + numtotal_2nd_additional
        )
        |> pivot_longer(-week_end,names_to="type",values_to="value")
        |> filter(type %in% c("numtotal_additional","numtotal_fully"))
        |> filter(!is.na(value))
        |> group_by(type)
        |> mutate(daydiff = diff(c(week_end,0))
                  , boosterdiff = diff(c(value,0))
                  , daily_rate = as.integer(boosterdiff/as.numeric(daydiff))
                  , days = as.numeric(week_end - day0)
                  , week_end = week_end
        )
        |> select(days,daily_rate,week_end)
        |> filter(between(days,0,159 + offset0))
        #	|> filter(between(days,-20,159))
)

saveVars(can)
