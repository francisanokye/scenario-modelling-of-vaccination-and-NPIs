library(tidyverse)
library(lubridate)
library(shellpipes)

loadEnvironments()

analysis_start <- as.Date("2021-11-08")
analysis_end   <- as.Date("2022-05-22")

init_date        <- analysis_start - days(offset)


actual_boosters <- boosterframe |>
  arrange(date) |>
  filter(date >= analysis_start, date <= analysis_end) |>
  mutate(
    day = as.numeric(date - init_date),
    scenario = "Actual Boosters",
    daily_rate = daily_rate,
    cumulative_total = cumsum(daily_rate)
  ) |>
  select(date, day, scenario, daily_rate, cumulative_total)


rdsSave(actual_boosters)

saveEnvironment()
