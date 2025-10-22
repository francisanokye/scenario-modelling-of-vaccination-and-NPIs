# load libraries
library(shellpipes)
library(zoo)
library(tidyverse)

# load saved environments
loadEnvironments()

# define study period
start_date <- "2021-12-15"
last_date  <- "2022-07-15"  # "2022-05-26"

# read and format seroprevalence data
seroprevalence <- csvRead()
seroprevalence <- seroprevalence |>
  mutate(week_end = as.Date(week_end, format = "%Y-%m-%d")) |>
  # filter anti-nucleocapsid antibodies (infection-only)
  filter((geo == "NL") & (ab_estimate == "Anti-N estimate")) |>
  rename_at("week_end", ~"date") |>
  # keep only date and mean percentage
  select(date, pct_mean) |>
  # create index
  mutate(day = 1:n())

print(seroprevalence)

# convert to seroprevalence matrix format for modeling
serodat <- seroprevalence |>
  transmute(
    NULL,
    date,
    time   = as.numeric(date - as.Date(start_date)) + 1 + offset0,
    matrix = "serop",
    value  = pct_mean
  ) |>
  filter(between(date, as.Date(start_date), as.Date(last_date)))

# save formatted matrix
rdsSave(serodat)
