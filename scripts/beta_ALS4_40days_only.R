library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

stepkeep <- 149

dat = rdsRead("build_multipliers")
cal = rdsRead("calibrate")
spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

obs_newR = data.frame(
  time = cal$cal_spec$integers$obs_times_log_newR
  , value = cal$cal_spec$default$obs_log_newR
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
  , dat$beta_ALS4_40days_only
  , calibrated_beta$value[(1+stepkeep+nrow(dat)):time_steps]
)
stopifnot(length(new_beta) == time_steps)
# 
print(new_beta)

rdsSave(new_beta)
