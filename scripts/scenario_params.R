library(shellpipes)

loadEnvironments()

## scenario params to constrain fitting

cutoff_time <- 232
intervention_removed <- as.Date("2022-03-14")

step_start <- as.numeric(omicron_start - burn_in_start) + 1
step_end   <- as.numeric(intervention_removed - burn_in_start) + 1

saveEnvironment()
