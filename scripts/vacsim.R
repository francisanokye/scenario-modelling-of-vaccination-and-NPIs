library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

calib <- rdsRead("calibrate")
# newbeta <- rdsRead(pipeStar())


print(new_booster)

quit()
early_booster_changepoints <- as.integer(vax_early_rate$day)
early_booster_values <- as.numeric(vax_early_rate$daily_rate)

quit()

if (!any(early_booster_changepoints == 0)) {
  early_booster_changepoints <- c(0L, early_booster_changepoints)
  early_booster_values       <- c(0,  early_booster_values)
}
spec_early <- mp_optimized_spec(calibrator, "modified")
spec_early$integers$booster_vac_changepoints <- as.integer(early_booster_changepoints)
spec_early$default$booster_vac_values        <- as.numeric(early_booster_values)



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
