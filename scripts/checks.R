library(dplyr)
library(tidyr)
library(shellpipes)

loadEnvironments()

study_start <- as.Date("2021-12-15")
study_end <- as.Date("2022-12-15")

combosims <- bind_rows(rdsReadList()) |>
  mutate(
    date = as.Date(burn_in_start) + (time - 1),
    vax_scenario = recode(
      as.character(vac),
      "vax_actual_boosters" = "Actual NL"
    ),
    policy = recode(
      as.character(beta),
      "base" = "Actual NPIs"
    )
  )

infection_contributions <- combosims |>
  filter(
    vax_scenario == "Actual NL",
    policy == "Actual NPIs",
    date >= study_start,
    date <= study_end,
    matrix %in% c("incS", "incv2", "incv3")
  ) |>
  group_by(matrix) |>
  summarise(
    infections = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    p_hosp = case_when(
      matrix == "incS"  ~ 0.05,
      matrix == "incv2" ~ 0.015,
      matrix == "incv3" ~ 0.003
    ),
    hospitalizations = infections * p_hosp
  )

print(infection_contributions)

print(
  infection_contributions |>
    summarise(
      total_infections = sum(infections),
      total_hospitalizations = sum(hospitalizations),
      effective_IHR = total_hospitalizations / total_infections
    )
)


cum_check <- combosims |>
  filter(
    vax_scenario == "Actual NL",
    policy == "Actual NPIs",
    matrix == "cum_hosp_total"
  ) |>
  arrange(date)

baseline <- cum_check |>
  filter(date < study_start) |>
  slice_max(date, n = 1) |>
  pull(value)

endpoint <- cum_check |>
  filter(date <= study_end) |>
  slice_max(date, n = 1) |>
  pull(value)

cat("Baseline cumulative hospitalizations:", baseline, "\n")
cat("Endpoint cumulative hospitalizations:", endpoint, "\n")
cat("Study-period hospitalizations:", endpoint - baseline, "\n")