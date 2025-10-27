library(macpan2)
library(ggplot2)
library(dplyr)
library(shellpipes)
rpcall("policy_neutral.simulate.Rout simulate.R policy_neutral.rds build_multipliers.rds calibrate.rds")
rpcall("base_beta.simulate.Rout simulate.R base_beta.rds build_multipliers.rds calibrate.rds")

## We want to create a new specs without splines and simulate
## the trajectories by inputing our own time series of beta_t

dat = rdsRead("build_multipliers")
cal = rdsRead("calibrate")

## This is loading new_beta according to what explicit target we want
new_beta <- rdsRead(pipeStar())

spec = mp_optimized_spec(cal, spec_structure = "modified")
time_steps = cal$simulator$tmb_model$time_steps$time_steps

obs_serop = data.frame(
    time = cal$cal_spec$integers$obs_times_logit_serop
  , value = cal$cal_spec$default$obs_logit_serop
)

trajs = (spec
  |> mp_tmb_delete(phase = "before", at = 23L)
  |> mp_tmb_delete(phase = "before", at = 23L)
  |> mp_tmb_insert(
    default = list(timeseries_beta = new_beta)
  )
  |> mp_simulator(time_steps, "logit_serop")
  |> mp_trajectory()
)

  (trajs
    |> filter(matrix == "logit_serop")
    |> ggplot()
    + aes(time, plogis(value))
    + geom_line()
    + geom_point(data = obs_serop)
)

newserop <- (trajs
	|> transmute(time
		, serop = plogis(value)
	)
	|> left_join(obs_serop)
	|> mutate(obs_serop = plogis(value))
	|> select(-value)
)

print(newserop)

