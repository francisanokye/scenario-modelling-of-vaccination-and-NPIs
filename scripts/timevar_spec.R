# library(macpan2)
# library(dplyr)
# library(zoo)
# library(shellpipes)
# library(tidyverse)
# 
# loadEnvironments()
# 
# # ============================================================
# init_date <- as.Date("2021-07-24")
# reporting_delay <- TRUE
# 
# # ------------------------------------------------------------
# # Vaccination inputs
# # ------------------------------------------------------------
# double_daily_vac <- nl_vax |>
#   dplyr::filter(type == "numtotal_fully") |>
#   dplyr::mutate(days = as.integer(days)) |>
#   dplyr::arrange(days)
# 
# booster_daily_vac <- booster_inputs_df |>
#   dplyr::filter(scenario == "Actual Boosters") |>
#   dplyr::transmute(days = as.integer(day), daily_rate = as.numeric(booster_vac_value)) |>
#   dplyr::arrange(days) |>
#   dplyr::group_by(days) |>
#   dplyr::summarise(
#     daily_rate = sum(daily_rate, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # ------------------------------------------------------------
# # Initial state
# # ------------------------------------------------------------
# initialize_state <- list(
#   A1 ~ 1, E1 ~ 1, I1 ~ 1, R1 ~ 1,
#   A2 ~ 1, E2 ~ 1, I2 ~ 1, R2 ~ 1,
#   A3 ~ 1, E3 ~ 1, I3 ~ 1, R3 ~ 1,
#   V2 ~ 1 + V2_init,
#   V3 ~ 1,
#   cum_double_vac ~ 0,
#   cum_booster ~ 0,
#   S1 ~ N - A1 - E1 - I1 - R1 - V2 - A2 - E2 - I2 - R2 - V3 - A3 - E3 - I3 - R3
# )
# 
# constant_computations <- list(
#   N ~ sum(S1, A1, E1, I1, R1, V2, A2, E2, I2, R2, V3, A3, E3, I3, R3)
# )
# 
# spec <- mp_tmb_model_spec(
#   before = c(initialize_state, constant_computations),
#   during = flows,
#   default = params
# )
# 
# # ------------------------------------------------------------
# # Time-varying vaccination rates
# # ------------------------------------------------------------
# double_vac_changepoints <- as.integer(double_daily_vac$days)
# double_vac_values <- as.numeric(double_daily_vac$daily_rate)
# 
# booster_vac_changepoints <- as.integer(booster_daily_vac$days)
# booster_vac_values <- as.numeric(booster_daily_vac$daily_rate)
# 
# if (!any(double_vac_changepoints == 0)) {
#   double_vac_changepoints <- c(0L, double_vac_changepoints)
#   double_vac_values <- c(0, double_vac_values)
# }
# 
# if (!any(booster_vac_changepoints == 0)) {
#   booster_vac_changepoints <- c(0L, booster_vac_changepoints)
#   booster_vac_values <- c(0, booster_vac_values)
# }
# 
# expr2 <- list(vac2 ~ time_var(double_vac_values, double_vac_changepoints))
# expr3 <- list(vac3 ~ time_var(booster_vac_values, booster_vac_changepoints))
# 
# timevar_spec <- spec |>
#   mp_tmb_insert(
#     phase = "during",
#     at = 1,
#     expressions = c(expr2, expr3),
#     default = list(
#       double_vac_values = double_vac_values,
#       booster_vac_values = booster_vac_values
#     ),
#     integers = list(
#       double_vac_changepoints = double_vac_changepoints,
#       booster_vac_changepoints = booster_vac_changepoints
#     )
#   )
# 
# timevar_spec <- mp_tmb_insert(
#   timevar_spec,
#   at = Inf,
#   expressions = list(
#     newR ~ newR_total
#   )
# )
# 
# rdsSave(timevar_spec)







library(macpan2)
library(dplyr)
library(zoo)
library(shellpipes)
library(tidyverse)

loadEnvironments()

# ============================================================
init_date        <- as.Date("2021-07-24")
reporting_delay  <- TRUE

# ------------------------------------------------------------
# Vaccination inputs
# ------------------------------------------------------------
double_daily_vac <- nl_vax |>
  dplyr::filter(type == "numtotal_fully") |>
  dplyr::mutate(days = as.integer(days)) |>
  dplyr::arrange(days)

booster_daily_vac <- booster_inputs_df |>
  dplyr::filter(scenario == "Actual Boosters") |>
  dplyr::transmute(days = as.integer(day), daily_rate = as.numeric(booster_vac_value)) |>
  dplyr::arrange(days) |>
  dplyr::group_by(days) |>
  dplyr::summarise(
    daily_rate = sum(daily_rate, na.rm = TRUE),
    .groups    = "drop"
  )

# ------------------------------------------------------------
# Initial state
# ------------------------------------------------------------
initialize_state <- list(
  A1 ~ 1, E1 ~ 1, I1 ~ 1, R1 ~ 1,
  A2 ~ 1, E2 ~ 1, I2 ~ 1, R2 ~ 1,
  A3 ~ 1, E3 ~ 1, I3 ~ 1, R3 ~ 1,
  V2 ~ 1 + V2_init,
  V3 ~ 1,
  cum_double_vac ~ 0,
  cum_booster    ~ 0,
  S1 ~ N - A1 - E1 - I1 - R1 - V2 - A2 - E2 - I2 - R2 - V3 - A3 - E3 - I3 - R3
)

constant_computations <- list(
  N ~ sum(S1, A1, E1, I1, R1, V2, A2, E2, I2, R2, V3, A3, E3, I3, R3)
)

spec <- mp_tmb_model_spec(
  before  = c(initialize_state, constant_computations),
  during  = flows,
  default = params
)

# ------------------------------------------------------------
# Time-varying vaccination rates
# ------------------------------------------------------------
double_vac_changepoints  <- as.integer(double_daily_vac$days)
double_vac_values        <- as.numeric(double_daily_vac$daily_rate)

booster_vac_changepoints <- as.integer(booster_daily_vac$days)
booster_vac_values       <- as.numeric(booster_daily_vac$daily_rate)

if (!any(double_vac_changepoints == 0)) {
  double_vac_changepoints <- c(0L, double_vac_changepoints)
  double_vac_values       <- c(0,  double_vac_values)
}

if (!any(booster_vac_changepoints == 0)) {
  booster_vac_changepoints <- c(0L, booster_vac_changepoints)
  booster_vac_values       <- c(0,  booster_vac_values)
}

expr2 <- list(vac2 ~ time_var(double_vac_values,  double_vac_changepoints))
expr3 <- list(vac3 ~ time_var(booster_vac_values, booster_vac_changepoints))

timevar_spec <- spec |>
  mp_tmb_insert(
    phase       = "during",
    at          = 1,
    expressions = c(expr2, expr3),
    default     = list(
      double_vac_values  = double_vac_values,
      booster_vac_values = booster_vac_values
    ),
    integers    = list(
      double_vac_changepoints  = double_vac_changepoints,
      booster_vac_changepoints = booster_vac_changepoints
    )
  )

timevar_spec <- mp_tmb_insert(
  timevar_spec,
  at          = Inf,
  expressions = list(newR ~ newR_total)
)

rdsSave(timevar_spec)
