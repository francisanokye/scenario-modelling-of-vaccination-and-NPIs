library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

new_changepoints <- as.integer(new_booster$day)
new_values <- as.numeric(new_booster$daily_rate)

if (!any(new_changepoints == 0)) {
  new_changepoints <- c(0L, new_changepoints)
  new_values       <- c(0,  new_values)
}

newspec$integers$booster_vac_changepoints <- as.integer(new_changepoints)
newspec$default$booster_vac_values <- as.numeric(new_values)

# spec_natl <- spec_natl |> 
#   mp_tmb_insert(default = list(V2_init = new_V2_init))

newsim <- (newspec
  |> mp_simulator(time_steps = sim_horizon, outputs = model_outputs)
  |> mp_trajectory()
	|> mutate(vac = pipeStar()
	 	, beta = targetname()
		, beta = ifelse(grepl("base",beta),"base",beta)
		, beta = ifelse(grepl("strict_npis",beta),"strict_npis",beta)
		, beta = ifelse(grepl("no_npis",beta),"no_npis",beta)
	)
)

rdsSave(newsim)
