library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)

loadEnvironments()

cal = rdsRead("calibrate")
spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

calibrated_beta = (spec
	|> mp_simulator(time_steps, "beta")
	|> mp_trajectory()
)

base_beta <- calibrated_beta$value

stopifnot(length(base_beta) == time_steps)

print(base_beta)

rdsSave(base_beta)
