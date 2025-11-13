library(macpan2)
library(dplyr)
library(zoo)
library(shellpipes)

loadEnvironments()


# ---------------------------------------
# Vaccination data
# ---------------------------------------
double_daily_vac  <- nl_additional |> dplyr::filter(type == "numtotal_fully")
booster_daily_vac <- nl_additional |> dplyr::filter(type == "numtotal_additional")

reporting_delay <- TRUE
# ---------------------------------------
# Initialize model
# ---------------------------------------
# initialize_state <- list(
#   A1 ~ 1, E1 ~ 1, I1 ~ 1, R1 ~ 1,
#   A2 ~ 1, E2 ~ 1, I2 ~ 1, R2 ~ 1,
#   A3 ~ 1, E3 ~ 1, I3 ~ 1, R3 ~ 1,
#   V2 ~ 1 + V2_init,
#   V3 ~ 1,
#   S1 ~ N - A1 - E1 - I1 - R1 - V2 - A2 - E2 - I2 - R2 - V3 - A3 - E3 - I3 - R3
# )
# 
# constant_computations <- list(
#   N ~ sum(S1, A1, E1, I1, R1, V2, A2, E2, I2, R2, V3, A3, E3, I3, R3)
# )

initialize_state <- list(
  # seeds for infection states (all cohorts)
  A1 ~ 1, E1 ~ 1, I1 ~ 1, R1 ~ 1,
  A2 ~ 1, E2 ~ 1, I2 ~ 1, R2 ~ 1,
  A3 ~ 1, E3 ~ 1, I3 ~ 1, R3 ~ 1,
  
  # vaccination cohort 2 (protected/waned split)
  V2_P ~ 1 + V2_init,
  V2_W ~ 0,
  
  # vaccination cohort 3 (protected/waned split)
  V3_P ~ 1,
  V3_W ~ 0,
  
  # susceptible baseline (fills the remainder)
  S1 ~ N - A1 - E1 - I1 - R1 - V2_P - V2_W - A2 - E2 - I2 - R2 - V3_P - V3_W - A3 - E3 - I3 - R3
)

constant_computations <- list(
  N ~ sum(S1, A1, E1, I1, R1, V2_P, V2_W, A2, E2, I2, R2, V3_P, V3_W, A3, E3, I3, R3)
)


spec <- mp_tmb_model_spec(
  before  = c(initialize_state, constant_computations),
  during  = flows,
  default = params
)

# ---------------------------------------
# Vaccination changepoints and scaling
# ---------------------------------------
double_vac_changepoints  <- as.integer(double_daily_vac$days  - (offset0 %% 7) - 1L)
double_vac_values        <- double_daily_vac$daily_rate

booster_vac_changepoints <- as.integer(booster_daily_vac$days - (offset0 %% 7) - 1L)
booster_vac_values       <- booster_daily_vac$daily_rate
if (!any(booster_vac_changepoints == 0)) {
  booster_vac_changepoints <- c(0, booster_vac_changepoints)
  booster_vac_values       <- c(0,  booster_vac_values)
}

expr2 <- list(vac2 ~ time_var(double_vac_values,  double_vac_changepoints))
expr3 <- list(vac3 ~ time_var(booster_vac_values, booster_vac_changepoints))

# ---------------------------------------
# insert time-varying vaccination and serop
# ---------------------------------------
timevar_spec <- spec |>
  mp_tmb_insert(
    phase = "during",
    at = 1,
    expressions = c(expr2, expr3),
    default  = list(
      double_vac_values  = double_vac_values,
      booster_vac_values = booster_vac_values
    ),
    integers = list(
      double_vac_changepoints  = double_vac_changepoints,
      booster_vac_changepoints = booster_vac_changepoints
    )
  ) |>
  mp_tmb_insert(
    at = Inf,
    expressions = list(serop ~ serop_total)
  )

rdsSave(timevar_spec)

