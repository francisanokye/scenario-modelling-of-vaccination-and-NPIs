library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

# calib <- rdsRead("calibrate")
# newbeta <- rdsRead(pipeStar())


print(new_booster)

new_changepoints <- as.integer(new_booster$day)
new_values <- as.numeric(new_booster$daily_rate)


if (!any(new_changepoints == 0)) {
  new_changepoints <- c(0L, new_changepoints)
  new_values       <- c(0,  new_values)
}

# newspec <- mp_optimized_spec(calib, "modified")
newspec$integers$booster_vac_changepoints <- as.integer(new_changepoints)
newspec$default$booster_vac_values <- as.numeric(new_values)

newsim <- (newspec
    |> mp_simulator(time_steps = sim_horizon, outputs = model_outputs)
    |> mp_trajectory()
	 |> mutate(vac = pipeStar())
)

print(newsim)

rdsSave(newsim)
