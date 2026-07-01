suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(shellpipes)
  library(macpan2)
  library(lubridate)
})

options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)

loadEnvironments()

vaccination_snapshot_rows <- vector("list", length(vax_scenario_specs))

for (j in seq_along(vax_scenario_specs)) {
  
  vax_name <- names(vax_scenario_specs)[j]
  spec_i   <- vax_scenario_specs[[vax_name]]
  
  traj <- spec_i |>
    mp_simulator(time_steps = sim_horizon, outputs = c("S1", "V2", "V3", "cum_booster")) |>
    mp_trajectory() |>
    mutate(
      date   = burn_in_start + (time - 1),
      matrix = as.character(matrix),
      value  = as.numeric(value)
    )
  
  start_date <- booster_start_dates[[vax_name]]
  
  state_at_start <- traj |>
    filter(date == start_date - 1) |>
    select(matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value)
  
  final_boosters <- traj |>
    filter(matrix == "cum_booster") |>
    filter(date == max(date)) |>
    pull(value)
  
  vaccination_snapshot_rows[[j]] <- tibble(
    vax_scenario                = vax_name,
    booster_start_date          = start_date,
    `No primary series`      = 100 * state_at_start$S1 / NL_POP,
    `Primary series only`    = 100 * state_at_start$V2 / NL_POP,
    `Boosters`               = 100 * state_at_start$V3 / NL_POP,
    `Final booster coverage` = 100 * final_boosters / NL_POP
  )
}

vaccination_snapshot <- bind_rows(vaccination_snapshot_rows) |>
  mutate(across(
    c(`No primary series`, `Primary series only`,
      `Boosters`, `Final booster coverage`),
    ~ round(.x, 1)
  ))

print(vaccination_snapshot)

write.csv(
  vaccination_snapshot,
  file.path(fig_dir, "..", "output", "vaccination_snapshot.csv"),
  row.names = FALSE
)