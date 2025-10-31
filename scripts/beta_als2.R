library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

## We want to create a new specs without splines and simulate
## the trajectories by inputing our own time series of beta_t

stepkeep <- 149

dat = rdsRead("build_multipliers")
cal = rdsRead("calibrate")
spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

obs_serop = data.frame(
  time = cal$cal_spec$integers$obs_times_logit_serop
  , value = cal$cal_spec$default$obs_logit_serop
)
calibrated_beta = (spec
                   |> mp_simulator(time_steps, "beta")
                   |> mp_trajectory()
)

if (interactive()) {
  (calibrated_beta
   |> ggplot()
   + aes(time, value)
   + geom_line()
  )
}

new_beta = c(
  calibrated_beta$value[1:stepkeep]
  , dat$beta_cf_ALS2
  , calibrated_beta$value[(1+stepkeep+nrow(dat)):time_steps]
)
stopifnot(length(new_beta) == time_steps)

print(new_beta)

rdsSave(new_beta)
