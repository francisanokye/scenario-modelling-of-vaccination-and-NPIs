suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(shellpipes)
})

options(
  dplyr.summarise.inform = FALSE,
  scipen = 9999
)

loadEnvironments()

NL_POP <- 510550
CA_POP <- 36991981

snapshot_beta <- "base"
time_origin <- if (exists("burn_in_start")) {
  as.Date(burn_in_start)
} else {
  as.Date("2021-12-15")
}

scenario_lookup <- tibble(
  vac_code = c(
    "vax_actual_boosters",
    "vax_early_boosters",
    "vax_faster_rate",
    "vax_slower_rate",
    "vax_national_rate"
  ),
  vax_scenario = c(
    "Actual Boosters",
    "Early Boosters",
    "Faster Boosters",
    "Slower Boosters",
    "National Rate"
  ),
  total_population = c(
    NL_POP,
    NL_POP,
    NL_POP,
    NL_POP,
    CA_POP
  )
)

required_matrices <- c("S1", "V2", "V3", "cum_booster")

trajectories <- bind_rows(rdsReadList()) |>
  filter(
    beta == snapshot_beta,
    matrix %in% required_matrices
  ) |>
  transmute(
    vac_code = as.character(vac),
    matrix = as.character(matrix),
    date = time_origin + as.integer(time) - 1L,
    value = as.numeric(value)
  )

missing_matrices <- setdiff(required_matrices, unique(trajectories$matrix))

if (length(missing_matrices) > 0L) {
  stop(
    "Missing required model outputs: ",
    paste(missing_matrices, collapse = ", ")
  )
}

booster_summary <- trajectories |>
  filter(matrix == "cum_booster") |>
  group_by(vac_code) |>
  summarise(
    booster_start_date = {
      positive_dates <- date[value > 0 & is.finite(value)]
      if (length(positive_dates) > 0L) min(positive_dates) else as.Date(NA)
    },
    final_boosters = {
      valid <- is.finite(value) & !is.na(date)
      if (any(valid)) value[valid][which.max(date[valid])] else NA_real_
    },
    .groups = "drop"
  )

missing_booster_scenarios <- booster_summary |>
  filter(is.na(booster_start_date)) |>
  pull(vac_code)

if (length(missing_booster_scenarios) > 0L) {
  warning(
    "No booster administration was found for: ",
    paste(missing_booster_scenarios, collapse = ", "),
    ". These scenarios were omitted."
  )
}

valid_booster_summary <- booster_summary |>
  filter(!is.na(booster_start_date)) |>
  mutate(snapshot_date = booster_start_date - 1L)

state_at_start <- trajectories |>
  filter(matrix %in% c("S1", "V2", "V3")) |>
  inner_join(
    valid_booster_summary |>
      select(vac_code, snapshot_date),
    by = c("vac_code", "date" = "snapshot_date")
  ) |>
  group_by(vac_code, matrix) |>
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = matrix,
    values_from = value
  )

vaccination_coverage <- valid_booster_summary |>
  select(
    vac_code,
    booster_start_date,
    final_boosters
  ) |>
  left_join(state_at_start, by = "vac_code") |>
  left_join(scenario_lookup, by = "vac_code") |>
  mutate(
    vax_scenario = coalesce(vax_scenario, vac_code),
    total_population = NL_POP,
    
    `No primary series (%)` =
      100 * S1 / NL_POP,
    
    `Primary series only (%)` =
      100 * V2 / NL_POP,
    
    `Boosters at rollout start (%)` =
      100 * V3 / NL_POP,
    
    `Final booster coverage (%)` =
      100 * final_boosters / NL_POP
  ) |>
  select(
    vax_scenario,
    vac_code,
    booster_start_date,
    total_population,
    `No primary series (%)`,
    `Primary series only (%)`,
    `Boosters at rollout start (%)`,
    `Final booster coverage (%)`
  ) |>
  mutate(
    across(
      where(is.numeric) & !all_of("total_population"),
      \(x) round(x, 1)
    )
  )

print(as.data.frame(vaccination_coverage))

write.csv(
  vaccination_coverage,
  "../output/vaccination_coverage_table.csv",
  row.names = FALSE
)