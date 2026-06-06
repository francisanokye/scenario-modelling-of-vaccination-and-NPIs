

library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

loadEnvironments()

stepkeep <- 150
als4_constant_value  <- 0.175

dat = df_est 
cal = rdsRead("calibrate")
spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

obs_newR = data.frame(time = cal$cal_spec$integers$obs_times_log_newR, 
                      value = cal$cal_spec$default$obs_log_newR
)

calibrated_beta = (spec
                   |> mp_simulator(time_steps, "beta")
                   |> mp_trajectory()
)

# uncomment to substitute the ALS-4 Only values
# calibrated_beta$value[stepkeep:length(calibrated_beta$value)] <- dat$beta_cf_ALS4

# uncomment to replace time 150 to 239 with a constant
calibrated_beta$value[stepkeep:239] <- als4_constant_value

new_beta <- calibrated_beta$value

stopifnot(length(new_beta) == time_steps)

print(new_beta)


rdsSave(new_beta)
