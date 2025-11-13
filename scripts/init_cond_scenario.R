library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

dat = rdsRead("seroprevdata")
cal = rdsRead("calibrate")


spec_fixed = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

print(spec_fixed$default$state)

library(outbreaks)
sir_spec = mp_tmb_library("starter_models"
                          , "sir"
                          , package = "macpan2"
)


sir_simulator = mp_simulator(sir_spec
                             , time_steps = 100
                             , outputs = c("S", "I", "R")
                             , default = list(N = 763, I = 10, beta = 1, gamma = 0.5)
)

traj_early <- sir_simulator |>
  mp_trajectory()

print(
        ggplot(traj_early, aes(x = time)) +
        geom_line(aes(y = value, color = matrix))
      + facet_wrap(~matrix)
      )

quit()


## --- PREP: make sure covariates exist as vectors of the right length ----
## time_steps must be your modeling horizon (days in st udy window)
## Example:
## time_steps <- as.integer(as.Date("2022-05-22") - as.Date("2021-12-15") + 1)

# 1) Ensure pristine series are present
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values

# 2) Coerce to numeric vectors of length = time_steps
as_vec_len <- function(x, L) {
  if (length(x) == 0L) return(rep(0, L))
  if (length(x) == 1L) return(rep(as.numeric(x), L))
  if (length(x) >= L)  return(as.numeric(x)[seq_len(L)])
  # if shorter than L, pad with zeros at the end
  c(as.numeric(x), rep(0, L - length(x)))
}

b_base <- as_vec_len(spec_fixed$default$booster_vac_values_base, time_steps)  # vac3
d_base <- as_vec_len(spec_fixed$default$double_vac_values_base,  time_steps)  # vac2

## --- CHOOSE A SAFE SHIFT -------------------------------------------------
k_req  <- 28L                           # desired early shift
k_days <- max(0L, min(k_req, length(b_base) - 1L))  # clamp to [0, n-1]

if (k_days == 0L) {
  warning("Early-booster shift reduced to 0 days because the series is too short.")
}

## --- BUILD THE EARLY SERIES (SHIFT LEFT BY k_days) -----------------------
## New day t gets what day t+k had; tail is zero-filled.
n <- length(b_base)
b_early <- c(b_base[(1 + k_days):n], rep(0, k_days))

## --- WARM-UP TO GET MODEL-IMPLIED PREWINDOW BOOSTS -----------------------
## If k_days==0, prewindow sum is zero and we skip warm-up automatically.
boost_prewindow <- 0
if (k_days > 0L) {
  spec_warm <- spec_fixed
  spec_warm$default$booster_vac_values <- b_early
  spec_warm$default$double_vac_values  <- d_base
  
  ## Request the booster flow as an output; name matches your flow: "booster_shot"
  traj_warm <- spec_warm |>
    mp_simulator(k_days, outputs = c("booster_shot")) |>
    mp_trajectory()

  ## Sum only the booster flow over the warm-up horizon
  boost_prewindow <- sum(traj_warm$value[traj_warm$matrix == "booster_shot"], na.rm = TRUE)
}

## --- PATCH INITIAL STATE (V2 -> V3 MOVE AT t0) ---------------------------
spec_early <- spec_fixed
spec_early$default$booster_vac_values <- b_early
spec_early$default$double_vac_values  <- d_base


quit()


pre_boost_name <- "V2"  # EDIT to your actual pre-booster compartment name(s)
boosted_name   <- "V3"  # EDIT to your actual boosted compartment name(s)

moveable <- min(boost_prewindow, state0[[pre_boost_name]])
state0[[pre_boost_name]] <- state0[[pre_boost_name]] - moveable
state0[[boosted_name]]   <- state0[[boosted_name]]   + moveable

spec_early$default$state0 <- state0

## --- RUN SCENARIO --------------------------------------------------------
traj_early <- spec_early |>
  mp_simulator(time_steps, "inc") |>
  mp_trajectory()
