library(macpan2)
library(shellpipes)
rpcall("simulate.Rout simulate.R build_multipliers.rds calibrate.rds")

## We want to create a new specs without splines and simulate the trajectories by inputing our own time series of beta_t
