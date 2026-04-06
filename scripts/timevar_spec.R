
library(macpan2)
library(dplyr)
library(zoo)
library(shellpipes)
library(tidyverse)
loadEnvironments()
booster_offset <- 150
day0 <- as.Date("2021-11-08") - lubridate::days(booster_offset)

double_daily_vac = nl_additional |> dplyr::filter(type == "numtotal_fully")
#booster_daily_vac = nl_additional |> dplyr::filter(type == "numtotal_additional")

booster_daily_vac <- weekly_scenarios_df |>
  dplyr::filter(scenario == "Actual Boosters") |>
  dplyr::select(week_end, scenario, weekly_rate) |>
  mutate(
    daydiff = as.numeric(week_end - lag(week_end, default = week_end[1] - 7)),
    daily_rate = weekly_rate / daydiff,
    days = as.numeric(week_end - day0)
  ) 

reporting_delay <- TRUE

initialize_state = list(
  A1 ~ 1, E1 ~ 1, I1 ~ 1, R1 ~ 1
  , A2 ~ 1, E2 ~ 1, I2 ~ 1, R2 ~ 1
  , A3 ~ 1, E3 ~ 1, I3 ~ 1, R3 ~ 1
  , V2 ~ 1 + V2_init
  , V3 ~ 1
  , S1 ~ N - A1 - E1 - I1 - R1 - V2 - A2 - E2 - I2 - R2 - V3 - A3 - E3 - I3 - R3
)

constant_computations = list(N ~ sum(S1, A1, E1, I1, R1, V2, A2, E2, I2, R2, V3, A3, E3, I3, R3))

spec <- mp_tmb_model_spec(
  before = c(initialize_state, constant_computations) 
  , during = flows
  , default = params
)

# we create a piece-wise time-varying vaccination rate based on true data
double_vac_changepoints = double_daily_vac$days - (offset0 %% 7) - 1
double_vac_values = double_daily_vac$daily_rate

booster_vac_changepoints = booster_daily_vac$days - (booster_offset %% 7) - 1
booster_vac_values = booster_daily_vac$daily_rate
if (!any(booster_vac_changepoints == 0)) {
  booster_vac_changepoints = c(0, booster_vac_changepoints)
  booster_vac_values = c(0, booster_vac_values)
}

# expressions
expr2 <- list(vac2 ~ time_var(double_vac_values, double_vac_changepoints))
expr3 <- list(vac3 ~ time_var(booster_vac_values, booster_vac_changepoints))

# update model specification to reflect the vaccination rates from data
timevar_spec <- spec |> mp_tmb_insert(
  phase = "during"
  , at = 1
  , expressions = c(expr2, expr3)
  , default = list(double_vac_values = double_vac_values, booster_vac_values = booster_vac_values)
  , integers = list(double_vac_changepoints = double_vac_changepoints, booster_vac_changepoints = booster_vac_changepoints)
)

timevar_spec = mp_tmb_insert(timevar_spec
                             , at = Inf
                             , expressions = list(newR ~ newR_total)
)

rdsSave(timevar_spec)