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
  library(forcats)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

init_date <- as.Date("2021-07-24")
fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sim_horizon <- 530
new_V2_init <- 277707

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

gamma_pmf_discrete <- function(mean, sd, max_lag = 40) {
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  tau <- 0:max_lag
  
  w <- pgamma(tau + 1, shape = shape, scale = scale) -
    pgamma(tau, shape = shape, scale = scale)
  
  w / sum(w)
}

make_delay_state_defaults <- function(prefix, weights) {
  n_lag <- length(weights) - 1
  setNames(
    as.list(rep(0, n_lag)),
    sprintf("%s_%02d", prefix, seq_len(n_lag) - 1)
  )
}

w_hosp  <- gamma_pmf_discrete(mean = 7,  sd = 4, max_lag = 40)
w_death <- gamma_pmf_discrete(mean = 18, sd = 7, max_lag = 60)

severity_state_defaults <- c(
  make_delay_state_defaults("H_leq2dose", w_hosp),
  make_delay_state_defaults("H_booster",  w_hosp),
  make_delay_state_defaults("D_leq2dose", w_death),
  make_delay_state_defaults("D_booster",  w_death)
)

add_severity_defaults <- function(spec) {
  spec$default$p_hosp1 <- p_hosp1
  spec$default$p_hosp2 <- p_hosp2
  spec$default$p_hosp3 <- p_hosp3
  
  spec$default$p_death1 <- p_death1
  spec$default$p_death2 <- p_death2
  spec$default$p_death3 <- p_death3
  
  for (nm in names(severity_state_defaults)) {
    spec$default[[nm]] <- severity_state_defaults[[nm]]
  }
  
  spec
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
  
  spec <- add_severity_defaults(spec)
  
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

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  policy_name <- sim_grid$policy_name[i]
  vax_name <- sim_grid$vax_name[i]
  
  spec_i <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[policy_name]]
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
      policy = policy_name,
      vax_scenario = vax_name
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

figure2_df <- bind_rows(
  severity_df |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max(cum_hosp_total, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Hospitalizations"),
  
  severity_df |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max(cum_death_total, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Deaths")
) |>
  mutate(
    policy_label = recode(
      as.character(policy),
      "No-ALS Only" = "No NPIs",
      "Actual ALS" = "Actual NPIs",
      "ALS-4 Only" = "Strict NPIs"
    ),
    policy_label = factor(policy_label, levels = c("No NPIs", "Actual NPIs", "Strict NPIs")),
    outcome = factor(outcome, levels = c("Hospitalizations", "Deaths")),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

figure2 <- ggplot(
  figure2_df,
  aes(x = vax_scenario, y = severe_outcome)
) +
  geom_col(fill = "grey65", colour = "black", width = 0.72, linewidth = 0.25) +
  facet_grid(outcome ~ policy_label, scales = "free_y") +
  scale_x_discrete(labels = vax_x_labs, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Cumulative severe outcomes"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 24, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 28, colour = "black"),
    strip.text = element_text(size = 26),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(size = 36, hjust = 0.5, face = "plain"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.5, "cm"),
    plot.background = element_blank()
  )

save_png(
  file.path(fig_dir, "figure2_severe_outcomes.png"),
  figure2,
  width = 9000,
  height = 5500
)