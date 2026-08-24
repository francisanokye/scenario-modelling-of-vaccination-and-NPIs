
library(tidyverse);theme_set(theme_bw())
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

source("vax_counterfactual_utils.R")

primary_series_target_date <- as.Date("2021-11-06")
primary_series_target <- round(0.767 * NL_POP)

# ------------------------------------------------------------
# Canadian initial two-dose condition
# ------------------------------------------------------------

use_new_V2_init <- TRUE

new_V2_init <- 277707

if (
  length(new_V2_init) != 1
  || !is.finite(new_V2_init)
  || new_V2_init <= 0
  || new_V2_init >= NL_POP
) {
  stop(
    "Invalid Canadian new_V2_init: "
    , paste(
      new_V2_init
      , collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Actual NL booster rollout
# ------------------------------------------------------------

new_booster <- make_rollout_scenario(
  boosterframe = boosterframe,
  burn_in_start = burn_in_start,
  booster_max = booster_max,
  rollout_speed = 1,
  shift_days = 0
)

vax_actual_rate_can_init = new_booster

saveEnvironment()