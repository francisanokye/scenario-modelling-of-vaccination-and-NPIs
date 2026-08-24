
library(macpan2)
library(tidyverse)
library(shellpipes)

loadEnvironments()

beta_scenario <- targetname()

newbeta <- if (grepl("strict_npis", beta_scenario)) {
  rdsRead("strict_npis")
} else if (grepl("no_npis", beta_scenario)) {
  rdsRead("no_npis")
} else {
  rdsRead("base")
}

newspec <- (newspec
  |> mp_tmb_delete(phase = "before", at = 26L)
  |> mp_tmb_insert(default = list(timeseries_beta = newbeta))
)

newsim <- (newspec
  |> mp_simulator(time_steps = sim_horizon , outputs = model_outputs)
  |> mp_trajectory()
  |> mutate(
   vac = vac_scenario
   , beta = beta_scenario
   , beta = ifelse(grepl("base", beta), "base", beta)
   , beta = ifelse(grepl("strict_npis", beta), "strict_npis", beta)
   , beta = ifelse(grepl("no_npis", beta), "no_npis", beta)
  )
)

rdsSave(newsim)







