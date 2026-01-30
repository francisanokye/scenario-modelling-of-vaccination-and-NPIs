library(tidyverse)
library(shellpipes)

loadEnvironments()

# number of days before the first data point to start simulations
offset0 <- 150

dat <- csvRead()
day0 <- as.Date("2021-12-07") - lubridate::days(offset0)

# NL daily vaccination rate
nl <- (dat
        |> filter(grepl("ewfoundland",prename))
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
        ##|> mutate(daily_rate = as.integer(daily_rate * (510550 / 36991981)))
        |> select(days, week_end, daily_rate)
        |> filter(between(days,0,159 + offset0))
)

nl$daily_rate[nl$daily_rate < 0] = 0

# Canada national daily average vaccination rate
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
        |> mutate(daily_rate = as.integer(daily_rate * (510550 / 36991981)))
        |> select(days, week_end, daily_rate)
        |> filter(between(days,0,159 + offset0)) 
)

can$daily_rate[can$daily_rate < 0] = 0


# assuming early vaccination occurred at the national average daily
cutoff <- as.Date("2021-12-18")

vax_early_at_can_rate <- nl %>%
  left_join(
    can %>% select(type, week_end, can_daily_rate = daily_rate),
    by = c("type", "week_end")
  ) %>%
  mutate(
    daily_rate_new = if_else(
      type == "numtotal_fully" & week_end < cutoff,
      coalesce(can_daily_rate, daily_rate),  # use CAN if available; otherwise keep NL
      daily_rate
    )
  ) %>%
  select(-can_daily_rate)

# rename to daily_rate:
vax_early_at_can_rate <- vax_early_at_can_rate %>%
  mutate(daily_rate = daily_rate_new) %>%
  select(-daily_rate_new)

print(nl,n = Inf)
print(can,n = Inf)
print(vax_early_at_can_rate, n = Inf)

saveVars(vax_early_at_can_rate)
# saveVars(nl)

