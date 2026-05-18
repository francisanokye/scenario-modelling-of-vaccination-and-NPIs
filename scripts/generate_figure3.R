suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(lubridate)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

init_date <- as.Date("2021-07-24")
analysis_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sim_horizon <- 530

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

calibrator <- rdsRead("calibrate")
padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base <- rdsRead("base_beta")
padded_beta_als4 <- rdsRead("beta_als4")

stopifnot(exists("booster_inputs_df"))

generate_scenario_boosters <- function(scenario_name) {
  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name) |>
    arrange(day) |>
    transmute(
      changepoints = as.integer(day),
      vac_values = as.numeric(booster_vac_value)
    )
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster_inputs_df rows found for scenario: ", scenario_name)
  }
  
  scenario_daily <- scenario_daily |>
    group_by(changepoints) |>
    summarise(vac_values = sum(vac_values, na.rm = TRUE), .groups = "drop") |>
    arrange(changepoints)
  
  if (!any(scenario_daily$changepoints == 0)) {
    scenario_daily <- bind_rows(
      tibble(changepoints = 0L, vac_values = 0),
      scenario_daily
    ) |>
      arrange(changepoints)
  }
  
  list(
    changepoints = scenario_daily$changepoints,
    vac_values = scenario_daily$vac_values
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
  spec <- mp_optimized_spec(calibrator, "modified")
  spec$integers$booster_vac_changepoints <- as.integer(booster_days)
  spec$default$booster_vac_values <- as.numeric(booster_vals)
  spec
}

update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  
  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

spec_actual <- generate_scenario_boosters("Actual Boosters") |>
  with(set_booster_vax(changepoints, vac_values))

spec_early <- generate_scenario_boosters("Early Boosters (-35 Days)") |>
  with(set_booster_vax(changepoints, vac_values))

spec_faster <- generate_scenario_boosters("Faster Rate") |>
  with(set_booster_vax(changepoints, vac_values))

spec_slower <- generate_scenario_boosters("Slower Rate") |>
  with(set_booster_vax(changepoints, vac_values))

vax_scenario_specs <- list(
  "Slower Rate" = spec_slower,
  "Actual Boosters" = spec_actual,
  "Faster Rate" = spec_faster,
  "Early Boosters (-35 Days)" = spec_early
)

policy_list <- list(
  "No-ALS Only" = padded_beta_noals,
  "Actual ALS" = padded_beta_base,
  "ALS-4 Only" = padded_beta_als4
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "Slower rate",
  "Actual boosters",
  "Faster rate",
  "Early boosters (-35 Days)"
)

vax_x_labs <- c(
  "Slower rate" = "Slower\nboosters",
  "Actual boosters" = "Actual\nboosters",
  "Faster rate" = "Faster\nboosters",
  "Early boosters (-35 Days)" = "Earlier\nboosters"
)

dose_fill_cols <- c(
  "Booster dose" = "maroon",
  "\u22642 dose" = "gold"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  spec_i <- update_vax_spec_beta(
    vax_scenario_specs[[sim_grid$vax_name[i]]],
    policy_list[[sim_grid$policy_name[i]]]
  )
  
  all_sim_results[[i]] <- spec_i |>
    mp_simulator(
      sim_horizon,
      outputs = c(
        "hosp_leq2dose",
        "hosp_booster",
        "hosp_total",
        "death_leq2dose",
        "death_booster",
        "death_total"
      )
    ) |>
    mp_trajectory() |>
    mutate(
      date = init_date + time,
      policy = sim_grid$policy_name[i],
      vax_scenario = sim_grid$vax_name[i]
    )
}

severity_df <- bind_rows(all_sim_results) |>
  mutate(
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      vax_scenario,
      "Slower Rate" = "Slower rate",
      "Actual Boosters" = "Actual boosters",
      "Faster Rate" = "Faster rate",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  ) |>
  pivot_wider(
    names_from = matrix,
    values_from = value,
    values_fill = 0
  ) |>
  filter(
    date >= analysis_start,
    date <= analysis_end
  ) |>
  arrange(policy, vax_scenario, date) |>
  group_by(policy, vax_scenario) |>
  mutate(
    cum_hosp_leq2dose = cumsum(hosp_leq2dose),
    cum_hosp_booster = cumsum(hosp_booster),
    cum_hosp_total = cumsum(hosp_total),
    
    cum_death_leq2dose = cumsum(death_leq2dose),
    cum_death_booster = cumsum(death_booster),
    cum_death_total = cumsum(death_total)
  ) |>
  ungroup()

diagnostic_summary <- severity_df |>
  group_by(policy, vax_scenario) |>
  summarise(
    hosp_total = max(cum_hosp_total, na.rm = TRUE),
    hosp_leq2dose = max(cum_hosp_leq2dose, na.rm = TRUE),
    hosp_booster = max(cum_hosp_booster, na.rm = TRUE),
    death_total = max(cum_death_total, na.rm = TRUE),
    death_leq2dose = max(cum_death_leq2dose, na.rm = TRUE),
    death_booster = max(cum_death_booster, na.rm = TRUE),
    .groups = "drop"
  )

print(diagnostic_summary)

figure2_df <- bind_rows(
  severity_df |>
    group_by(policy, vax_scenario) |>
    summarise(
      leq2dose = max(cum_hosp_leq2dose, na.rm = TRUE),
      booster = max(cum_hosp_booster, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Hospitalizations"),
  
  severity_df |>
    group_by(policy, vax_scenario) |>
    summarise(
      leq2dose = max(cum_death_leq2dose, na.rm = TRUE),
      booster = max(cum_death_booster, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Deaths")
) |>
  pivot_longer(
    cols = c(leq2dose, booster),
    names_to = "dose_group",
    values_to = "severe_outcome"
  ) |>
  mutate(
    dose_group = recode(
      dose_group,
      "leq2dose" = "\u22642 dose",
      "booster" = "Booster dose"
    ),
    dose_group = factor(
      dose_group,
      levels = c("\u22642 dose", "Booster dose")
    ),
    policy_label = recode(
      as.character(policy),
      "No-ALS Only" = "No NPIs",
      "Actual ALS" = "Actual NPIs",
      "ALS-4 Only" = "Strict NPIs"
    ),
    policy_label = factor(
      policy_label,
      levels = c("No NPIs", "Actual NPIs", "Strict NPIs")
    ),
    outcome = factor(
      outcome,
      levels = c("Hospitalizations", "Deaths")
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

figure2 <- ggplot(
  figure2_df,
  aes(
    x = vax_scenario,
    y = severe_outcome,
    fill = dose_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_grid(
    outcome ~ policy_label,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = dose_fill_cols,
    breaks = c("\u22642 dose","Booster dose")
  ) +
  scale_x_discrete(
    labels = vax_x_labs,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Severe infections",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, colour = "black"),
    strip.text = element_text(size = 30, face = "plain"),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm"),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.background = element_blank()
  )

save_png(
  file.path(fig_dir, "figure3.png"),
  figure2,
  width = 8000,
  height = 3500
)
