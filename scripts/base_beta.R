
library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

loadEnvironments()

new_constant_value <- 0.3675
cutoff_time <- 232

cal = rdsRead("calibrate")
spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

calibrated_beta = (spec
 |> mp_simulator(time_steps, "beta")
 |> mp_trajectory()
)

base_beta <- calibrated_beta$value

base_beta[cutoff_time:length(base_beta)] <- new_constant_value

stopifnot(length(base_beta) == time_steps)

rdsSave(base_beta)
