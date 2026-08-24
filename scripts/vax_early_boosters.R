
library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()
source("vax_counterfactual_utils.R")

days_earlier <- 35

new_booster <- make_rollout_scenario(
  boosterframe = boosterframe,
  burn_in_start = burn_in_start,
  booster_max = booster_max,
  rollout_speed = 1,
  shift_days = days_earlier
)

vax_early_boosters = new_booster

print(new_booster)
saveEnvironment()
