suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
  library(forcats)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

init_date <- as.Date("2021-07-24")
analysis_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")
release_date <- as.Date("2022-03-17")

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

calibrator <- rdsRead("calibrate")
padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base <- rdsRead("base_beta")
padded_beta_als4 <- rdsRead("beta_als4")

stopifnot(exists("booster_inputs_df"))

scenario_cols <- c(
  "Slower boosters" = "#CC79A7",
  "Actual boosters" = "#4D4D4D",
  "Earlier boosters" = "#D55E00",
  "Faster boosters" = "#0072B2"
)

scenario_dose_cols <- list(
  "Actual boosters" = c(
    "≤2 dose" = "#BDBDBD",
    "Booster dose" = "#4D4D4D"
  ),
  "Slower rate" = c(
    "≤2 dose" = "#E8B6D2",
    "Booster dose" = "#CC79A7"
  ),
  "Faster rate" = c(
    "≤2 dose" = "#9ECAE1",
    "Booster dose" = "#0072B2"
  ),
  "Early boosters (-35 Days)" = c(
    "≤2 dose" = "#FDBB84",
    "Booster dose" = "#D55E00"
  )
)

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

spec_natl <- generate_scenario_boosters("National Rate") |>
  with(set_booster_vax(changepoints, vac_values)) |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))

vax_scenario_specs <- list(
  "Actual Boosters" = spec_actual,
  "Slower Rate" = spec_slower,
  "Faster Rate" = spec_faster,
  "Early Boosters (-35 Days)" = spec_early,
  "National Rate" = spec_natl
)

policy_list <- list(
  "Actual ALS" = padded_beta_base,
  "ALS-4 Only" = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
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
        "cum_booster",
        "inc",
        "inc_leq2dose",
        "inc_booster",
        "cum_inc_total",
        "cum_inc_leq2dose",
        "cum_inc_booster"
      )
    ) |>
    mp_trajectory() |>
    mutate(
      date = init_date + time,
      policy = policy_name,
      vax_scenario = vax_name
    )
}

results_df <- bind_rows(all_sim_results) |>
  mutate(
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      vax_scenario,
      "National Rate" = "National rate",
      "Slower Rate" = "Slower rate",
      "Actual Boosters" = "Actual boosters",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
      "Faster Rate" = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

incidence_df <- results_df |>
  filter(date >= analysis_start)

cumulative_df <- results_df |>
  filter(date >= analysis_start)

ab_ymax <- incidence_df |>
  filter(
    vax_scenario == "Actual boosters",
    policy %in% c("No-ALS Only", "ALS-4 Only"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(policy, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max(total_inc, na.rm = TRUE)) |>
  pull(ymax)

gh_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Early boosters (-35 Days)", "Actual boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max(total_inc, na.rm = TRUE)) |>
  pull(ymax)

de_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Slower rate", "Faster rate"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max(total_inc, na.rm = TRUE)) |>
  pull(ymax)

get_panel_cols <- function(scenario_name) {
  scenario_dose_cols[[as.character(scenario_name)]]
}

make_shaded_incidence_plot <- function(
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    legend_position = c(0.50, 0.85),
    y_limit = NULL,
    boundary_linetype = "solid"
) {
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      policy == policy_name,
      vax_scenario == scenario_name,
      matrix %in% c("inc_leq2dose", "inc_booster")
    ) |>
    select(date, matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value, values_fill = 0) |>
    mutate(total_inc = inc_leq2dose + inc_booster)
  
  y_upper <- if (is.null(y_limit)) max(plot_df$total_inc, na.rm = TRUE) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc, fill = "≤2 dose")) +
    geom_area(aes(y = inc_booster, fill = "Booster dose")) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 1.1, linetype = boundary_linetype) +
    geom_line(aes(y = inc_booster), colour = "black", linewidth = 1.1, linetype = boundary_linetype) +
    geom_vline(xintercept = release_date, linetype = "dashed", linewidth = 1.5) +
    scale_fill_manual(values = panel_cols, breaks = c("Booster dose", "≤2 dose")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y", expand = c(0.01, 0)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(title = panel_title, x = NULL, y = "Daily incidence", fill = NULL, tag = panel_tag) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = legend_position,
      legend.justification = "center",
      legend.direction = "vertical",
      legend.text = element_text(size = 35),
      legend.background = element_blank(),
      plot.title = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 35, 12, 12)
    )
}

make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      vax_scenario == scenario_name,
      policy %in% c("No-ALS Only", "ALS-4 Only"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "≤2 dose",
        "cum_inc_booster" = "Booster dose"
      ),
      policy_label = recode(
        as.character(policy),
        "No-ALS Only" = "No NPIs",
        "ALS-4 Only" = "Strict NPIs"
      ),
      dose_group = factor(dose_group, levels = c("≤2 dose", "Booster dose")),
      policy_label = factor(policy_label, levels = c("No NPIs", "Strict NPIs"))
    )
  
  ymax <- max(plot_df$value, na.rm = TRUE)
  
  ggplot(plot_df, aes(date, value, colour = dose_group, linetype = policy_label,
                      group = interaction(dose_group, policy_label))) +
    geom_line(linewidth = 3, alpha = 1) +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.25,
             label = "≤2 dose + No NPIs (Solid)", hjust = 0, size = 10,
             colour = panel_cols["≤2 dose"], fontface = "bold", alpha = 1) +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.19,
             label = "Booster dose + No NPIs (Solid)", hjust = 0, size = 10,
             colour = panel_cols["Booster dose"], fontface = "bold", alpha = 1) +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.13,
             label = "≤2 dose + Strict NPIs (Dashed)", hjust = 0, size = 10,
             colour = panel_cols["≤2 dose"], fontface = "bold", alpha = 1) +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.07,
             label = "Booster dose + Strict NPIs (Dashed)", hjust = 0, size = 10,
             colour = panel_cols["Booster dose"], fontface = "bold", alpha = 1) +
    scale_colour_manual(values = panel_cols) +
    scale_linetype_manual(values = c("No NPIs" = "solid", "Strict NPIs" = "twodash")) +
    guides(colour = "none", linetype = "none") +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y",
                 expand = expansion(mult = c(0.01, 0.08))) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Cumulative infections",
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "none",
      plot.title = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  plot_df <- df |>
    filter(
      policy == "ALS-4 Only",
      vax_scenario %in% c("Slower rate", "Faster rate"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "≤2 dose",
        "cum_inc_booster" = "Booster dose"
      ),
      rollout_speed = recode(
        as.character(vax_scenario),
        "Slower rate" = "Slower boosters",
        "Faster rate" = "Faster boosters"
      ),
      dose_group = factor(dose_group, levels = c("≤2 dose", "Booster dose")),
      rollout_speed = factor(rollout_speed, levels = c("Slower boosters", "Faster boosters"))
    )
  
  ymax <- max(plot_df$value, na.rm = TRUE)
  
  ggplot(plot_df, aes(date, value, colour = rollout_speed, linetype = dose_group,
                      group = interaction(rollout_speed, dose_group))) +
    geom_line(linewidth = 3) +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.25,
             label = "Slower boosters + ≤2 dose (Solid)", hjust = 0, size = 10,
             colour = scenario_cols["Slower boosters"], fontface = "bold") +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.19,
             label = "Slower boosters + Booster dose (Dashed)", hjust = 0, size = 10,
             colour = scenario_cols["Slower boosters"], fontface = "bold") +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.13,
             label = "Faster boosters + ≤2 dose (Solid)", hjust = 0, size = 10,
             colour = scenario_cols["Faster boosters"], fontface = "bold") +
    annotate("text", x = as.Date("2022-06-17"), y = ymax * 0.07,
             label = "Faster boosters + Booster dose (Dashed)", hjust = 0, size = 10,
             colour = scenario_cols["Faster boosters"], fontface = "bold") +
    scale_colour_manual(values = scenario_cols[c("Slower boosters", "Faster boosters")]) +
    scale_linetype_manual(values = c("≤2 dose" = "solid", "Booster dose" = "twodash")) +
    guides(colour = "none", linetype = "none") +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y",
                 expand = expansion(mult = c(0.01, 0.08))) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Cumulative infections",
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "none",
      plot.title = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

p_a <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "No-ALS Only",
  "(A)", "No NPIs", c(0.50, 0.85),
  y_limit = ab_ymax,
  boundary_linetype = "solid"
)

p_b <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "ALS-4 Only",
  "(B)", "Strict NPIs", c(0.20, 0.85),
  y_limit = ab_ymax,
  boundary_linetype = "twodash"
)

p_c <- make_cumulative_policy_plot(cumulative_df, "Actual boosters", "(C)")

p_d <- make_shaded_incidence_plot(
  incidence_df,
  "Slower rate",
  "ALS-4 Only",
  "(D)",
  "Slower booster doses",
  c(0.50, 0.85),
  y_limit = de_ymax,
  boundary_linetype = "twodash"
)

p_e <- make_shaded_incidence_plot(
  incidence_df,
  "Faster rate",
  "ALS-4 Only",
  "(E)",
  "Faster booster doses",
  c(0.20, 0.85),
  y_limit = de_ymax,
  boundary_linetype = "twodash"
)

p_f <- make_cumulative_speed_plot(cumulative_df, "(F)")

p_g <- make_shaded_incidence_plot(
  incidence_df, "Early boosters (-35 Days)", "ALS-4 Only",
  "(G)", "Early boosters", c(0.20, 0.85),
  y_limit = gh_ymax,
  boundary_linetype = "twodash"
)

p_h <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "ALS-4 Only",
  "(H)", "Actual booster doses", c(0.20, 0.85),
  y_limit = gh_ymax,
  boundary_linetype = "twodash"
)

plot_i_df <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Actual boosters", "Early boosters (-35 Days)"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(peak_daily_incidence = max(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    dose_group = recode(
      matrix,
      "inc_leq2dose" = "≤2 dose",
      "inc_booster" = "Booster dose"
    ),
    scenario_label = recode(
      as.character(vax_scenario),
      "Early boosters (-35 Days)" = "Earlier boosters",
      "Actual boosters" = "Actual boosters"
    ),
    scenario_label = factor(scenario_label, levels = c("Earlier boosters", "Actual boosters")),
    dose_group = factor(dose_group, levels = c("≤2 dose", "Booster dose"))
  )

p_i <- ggplot(
  plot_i_df,
  aes(x = dose_group, y = peak_daily_incidence, fill = scenario_label)
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.68,
    colour = "#4D4D4D",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = scenario_cols[c("Earlier boosters", "Actual boosters")]) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Earlier boosters result in higher peak incidence",
    x = NULL,
    y = "Peak daily incidence",
    fill = NULL,
    tag = "(I)"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 35, hjust = 0.5),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    legend.position = c(0.80, 0.90),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 35),
    legend.background = element_blank(),
    plot.title = element_text(size = 45, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text(size = 45, face = "bold"),
    plot.background = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )

figure1_3x3 <- (
  (p_a | p_b | p_c) /
    (p_d | p_e | p_f) /
    (p_g | p_h | p_i)
) +
  plot_layout(
    heights = c(1, 1, 1),
    widths = c(1, 1, 1)
  ) &
  theme(plot.margin = margin(8, 35, 8, 8))

save_png(
  file.path(fig_dir, "figure1.png"),
  figure1_3x3,
  width = 13500,
  height = 8000
)