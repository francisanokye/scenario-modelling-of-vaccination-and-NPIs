library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

calib <- rdsRead("calibrate")
newbeta <- rdsRead(pipeStar())

newspec <- (mp_optimized_spec(calib, "modified")
	|> mp_tmb_delete(phase = "before", at = 26L)
	|> mp_tmb_insert(default = list(timeseries_beta = newbeta))
)

newsim <- (newspec
    |> mp_simulator(time_steps = sim_horizon, outputs = model_outputs)
    |> mp_trajectory()
	 |> mutate(beta = pipeStar())
)

print(newsim)

rdsSave(newsim)
