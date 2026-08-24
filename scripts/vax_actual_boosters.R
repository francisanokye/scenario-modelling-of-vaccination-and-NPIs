
library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()
source("vax_counterfactual_utils.R")

new_booster <- make_rollout_scenario(
  boosterframe = boosterframe,
  burn_in_start = burn_in_start,
  booster_max = booster_max,
  rollout_speed = 1,
  shift_days = 0
)

vax_actual_boosters = new_booster

print(new_booster)
saveEnvironment()
