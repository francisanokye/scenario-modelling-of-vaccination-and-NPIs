suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(scales)
  library(lubridate)
  library(patchwork)
})

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

scenario_labels <- c(
  "vax_actual_boosters" = "Actual Boosters",
  "vax_early_boosters" = "Early Boosters"
)

policy_labels <- c(
  "base" = "Actual NPIs",
  "strict_npis" = "Strict NPIs",
  "no_npis" = "No NPIs"
)

required_outputs <- c(
  "S1", "V2", "V3", "E1", "E2", "E3", "A1", "A2", "A3",
  "I1", "I2", "I3", "R1", "R2", "R3", "incS", "incv2", "incv3",
  "inc", "cum_inc_total", "cum_hosp_total", "cum_death_total", "cum_booster"
)

sims <- bind_rows(rdsReadList()) |>
  mutate(
    date = burn_in_start + (time - 1),
    vaccination = recode(as.character(vac), !!!scenario_labels),
    policy = recode(as.character(beta), !!!policy_labels)
  ) |>
  filter(vaccination %in% unname(scenario_labels))

missing_outputs <- setdiff(required_outputs, unique(sims$matrix))
if (length(missing_outputs) > 0) {
  stop(
    "Add these variables to model_outputs in constrain_params.R and rebuild the simulations: ",
    paste(missing_outputs, collapse = ", ")
  )
}

wide <- sims |>
  filter(matrix %in% required_outputs) |>
  select(vaccination, policy, time, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value)

diagnostics <- wide |>
  mutate(
    susceptible_leq2 = S1 + V2,
    susceptible_booster = V3,
    exposed = E1 + E2 + E3,
    infectious = A1 + A2 + A3 + I1 + I2 + I3,
    recovered = R1 + R2 + R3,
    population_total = susceptible_leq2 + susceptible_booster + exposed + infectious + recovered,
    population_error = population_total - NL_POP,
    incidence_error = inc - (incS + incv2 + incv3)
  )

baseline_date <- omicron_start - days(1)

baseline <- wide |>
  filter(date <= baseline_date) |>
  group_by(vaccination, policy) |>
  slice_max(date, n = 1, with_ties = FALSE) |>
  select(
    vaccination, policy,
    baseline_inc = cum_inc_total,
    baseline_hosp = cum_hosp_total,
    baseline_death = cum_death_total
  )

study_totals <- wide |>
  filter(date >= omicron_start, date <= analysis_end) |>
  group_by(vaccination, policy) |>
  slice_max(date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(baseline, by = c("vaccination", "policy")) |>
  transmute(
    vaccination,
    policy,
    infections = cum_inc_total - baseline_inc,
    hospitalizations = cum_hosp_total - baseline_hosp,
    deaths = cum_death_total - baseline_death,
    administered_boosters = cum_booster
  )

state_df <- diagnostics |>
  filter(date >= omicron_start) |>
  select(vaccination, policy, date, susceptible_leq2, susceptible_booster, recovered) |>
  pivot_longer(
    cols = c(susceptible_leq2, susceptible_booster, recovered),
    names_to = "state",
    values_to = "value"
  ) |>
  mutate(
    state = recode(
      state,
      susceptible_leq2 = "Susceptible: ≤2 doses",
      susceptible_booster = "Susceptible: booster",
      recovered = "Recovered"
    )
  )

incidence_df <- wide |>
  filter(date >= omicron_start) |>
  select(vaccination, policy, date, incS, incv2, incv3) |>
  pivot_longer(c(incS, incv2, incv3), names_to = "infection_group", values_to = "value") |>
  mutate(
    infection_group = recode(
      infection_group,
      incS = "Unvaccinated/one dose",
      incv2 = "Two doses",
      incv3 = "Booster"
    )
  )

p_states <- ggplot(state_df, aes(date, value, colour = vaccination)) +
  geom_line(linewidth = 1) +
  facet_grid(state ~ policy, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "People", colour = NULL, title = "Population-state decomposition") +
  theme_clean() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

p_incidence <- ggplot(incidence_df, aes(date, value, colour = vaccination)) +
  geom_line(linewidth = 1) +
  facet_grid(infection_group ~ policy, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Daily infections", colour = NULL, title = "Infection decomposition by vaccination status") +
  theme_clean() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

png("../figures/early_actual_decomposition.png", width = 4200, height = 3600, res = 300, bg = "white")
print(p_states / p_incidence)
dev.off()