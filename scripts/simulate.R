library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)
rpcall("switch_als4_als2.simulate.Rout simulate.R switch_als4_als2.rds build_multipliers.rds calibrate.rds")
rpcall("policy_neutral.simulate.Rout simulate.R policy_neutral.rds build_multipliers.rds calibrate.rds")
rpcall("beta_als2.simulate.Rout simulate.R beta_als2.rds build_multipliers.rds calibrate.rds")
rpcall("beta_als4.simulate.Rout simulate.R beta_als4.rds build_multipliers.rds calibrate.rds")
rpcall("base_beta.simulate.Rout simulate.R base_beta.rds build_multipliers.rds calibrate.rds")

## We want to create a new specs without splines and simulate
## the trajectories by inputing our own time series of beta_t

dat = rdsRead("build_multipliers")
cal = rdsRead("calibrate")

## This is loading new_beta according to what explicit target we want
new_beta <- rdsRead(pipeStar())

spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

print(spec)

quit()



obs_newR = data.frame(
    time = cal$cal_spec$integers$obs_times_log_newR
  , value = cal$cal_spec$default$obs_log_newR
)

trajs = (spec
  |> mp_tmb_delete(phase = "before", at = 24L)
  |> mp_tmb_delete(phase = "before", at = 24L)
  |> mp_tmb_insert(default = list(timeseries_beta = new_beta))
  |> mp_simulator(time_steps, "log_newR")
  |> mp_trajectory()
)

  (trajs
    |> filter(matrix == "log_newR")
    |> ggplot()
    + aes(time, plogis(value))
    + geom_line()
    + geom_point(data = obs_newR)
)

newserop <- (trajs
	|> transmute(time
		, serop = plogis(value)
	)
	|> left_join(obs_newR)
	|> mutate(obs_newR = plogis(value))
	|> select(-value)
)

print(newserop)

