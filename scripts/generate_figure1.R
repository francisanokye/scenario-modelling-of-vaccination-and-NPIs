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

analysis_start <- as.Date("2021-12-15")
offset <- 150L
init_date <- analysis_start - offset

release_date <- as.Date("2022-03-17")
sim_horizon <- 530L
analysis_end <- init_date + sim_horizon - 1L

new_V2_init <- 277707
actual_booster_start <- as.Date("2021-11-08")

stopifnot(
  init_date + offset == analysis_start,
  analysis_start + (0 - offset) == init_date
)

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

calibrator <- rdsRead("calibrate")

padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

stopifnot(
  length(padded_beta_noals) == length(padded_beta_base),
  length(padded_beta_als4) == length(padded_beta_base),
  exists("booster_inputs_df")
)

available_early <- booster_inputs_df |>
  distinct(scenario) |>
  filter(grepl("^Early Boosters", scenario, ignore.case = TRUE)) |>
  pull(scenario)

if (length(available_early) != 1) {
  stop(
    "Expected exactly one Early Boosters scenario, but found:\n",
    paste(available_early, collapse = "\n")
  )
}

early_scenario_name <- available_early

early_shift_days <- suppressWarnings(
  as.integer(
    gsub(".*\\((-?[0-9]+) Days\\).*", "\\1", early_scenario_name)
  )
)

if (is.na(early_shift_days)) {
  early_shift_days <- -35L
}

early_start_date <- actual_booster_start + days(early_shift_days)
early_plot_label <- "Earlier boosters"

stopifnot(early_start_date >= init_date)

scenario_cols <- c(
  "Slower boosters"  = "magenta",
  "Actual boosters"  = "#4D4D4D",
  "Earlier boosters" = "blue",
  "Faster boosters"  = "orange",
  "National rate"    = "springgreen3"
)

scenario_dose_cols <- list(
  "Actual boosters" = c("≤2 dose" = "grey70", "Booster dose" = "#4D4D4D"),
  "Slower rate"    = c("≤2 dose" = "#E8B6D2", "Booster dose" = "magenta"),
  "Faster rate"    = c("≤2 dose" = "#FDBB84", "Booster dose" = "orange")
)

scenario_dose_cols[[early_plot_label]] <- c(
  "≤2 dose" = "#9ECAE1",
  "Booster dose" = "blue"
)

booster_start_dates <- c(
  "Actual Boosters" = actual_booster_start,
  setNames(early_start_date, early_scenario_name),
  "Faster Rate" = actual_booster_start,
  "Slower Rate" = actual_booster_start,
  "National Rate" = actual_booster_start
)

generate_scenario_boosters <- function(scenario_name) {
  
  scenario_match <- booster_inputs_df |>
    distinct(scenario) |>
    filter(tolower(scenario) == tolower(scenario_name)) |>
    pull(scenario)
  
  if (length(scenario_match) != 1) {
    stop(
      "Scenario name not found: ", scenario_name,
      "\nAvailable scenarios are:\n",
      paste(unique(booster_inputs_df$scenario), collapse = "\n")
    )
  }
  
  scenario_name_actual <- scenario_match
  
  start_lookup <- tibble(
    scenario = names(booster_start_dates),
    start_date = as.Date(booster_start_dates),
    scenario_clean = tolower(names(booster_start_dates))
  ) |>
    filter(scenario_clean == tolower(scenario_name_actual))
  
  if (nrow(start_lookup) != 1) {
    stop("No booster start date defined for scenario: ", scenario_name_actual)
  }
  
  start_day <- as.integer(start_lookup$start_date - init_date)
  
  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name_actual) |>
    transmute(
      changepoints = as.integer(day),
      vac_values = as.numeric(booster_vac_value)
    ) |>
    filter(changepoints >= start_day) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster input rows found for scenario: ", scenario_name_actual)
  }
  
  scenario_daily <- bind_rows(
    tibble(changepoints = 0L, vac_values = 0),
    scenario_daily
  ) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
  
  active_rows <- scenario_daily |> filter(vac_values > 0)
  
  if (nrow(active_rows) == 0) {
    stop("All booster values are zero for scenario: ", scenario_name_actual)
  }
  
  print(
    active_rows |>
      summarise(
        scenario = scenario_name_actual,
        intended_start_date = start_lookup$start_date,
        first_booster_day = min(changepoints),
        first_booster_date = init_date + first_booster_day,
        last_booster_day = max(changepoints),
        last_booster_date = init_date + last_booster_day,
        total_boosters = sum(vac_values),
        max_daily_boosters = max(vac_values)
      )
  )
  
  list(
    changepoints = scenario_daily$changepoints,
    vac_values = scenario_daily$vac_values
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
  stopifnot(length(booster_days) == length(booster_vals))
  stopifnot(all(diff(booster_days) >= 0))
  
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

spec_early <- generate_scenario_boosters(early_scenario_name) |>
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
  "National Rate" = spec_natl
)

vax_scenario_specs[[early_scenario_name]] <- spec_early

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
  early_plot_label,
  "Faster rate"
)

model_outputs <- c(
  "incS", "incv2", "incv3", "inc",
  "inc_leq2dose", "inc_booster",
  "cum_inc_leq2dose", "cum_inc_booster", "cum_inc_total",
  "hosp_leq2dose", "hosp_booster", "hosp_total",
  "death_leq2dose", "death_booster", "death_total",
  "cum_hosp_leq2dose", "cum_hosp_booster", "cum_hosp_total",
  "cum_death_leq2dose", "cum_death_booster", "cum_death_total",
  "cum_booster"
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
      time_steps = sim_horizon,
      outputs = model_outputs
    ) |>
    mp_trajectory() |>
    mutate(
      date = analysis_start + (time - offset),
      policy = policy_name,
      vax_scenario = vax_name
    )
}

results_df <- bind_rows(all_sim_results) |>
  mutate(
    date = as.Date(date),
    value = as.numeric(value),
    matrix = as.character(matrix),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      as.character(vax_scenario),
      "National Rate" = "National rate",
      "Slower Rate" = "Slower rate",
      "Actual Boosters" = "Actual boosters",
      "Faster Rate" = "Faster rate",
      .default = as.character(vax_scenario)
    ),
    vax_scenario = if_else(
      tolower(vax_scenario) == tolower(early_scenario_name),
      early_plot_label,
      vax_scenario
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

saveVars(results_df)

required_outputs <- c(
  "cum_inc_leq2dose",
  "cum_inc_booster",
  "cum_hosp_total",
  "cum_death_total"
)

missing_outputs <- setdiff(required_outputs, unique(results_df$matrix))

if (length(missing_outputs) > 0) {
  stop(
    "Missing model outputs: ",
    paste(missing_outputs, collapse = ", "),
    "\nRe-run upstream after updating flows and params."
  )
}

incidence_df <- results_df |>
  filter(date >= analysis_start, date <= analysis_end)

cumulative_df <- results_df |>
  filter(date >= analysis_start, date <= analysis_end)

severity_check <- results_df |>
  filter(
    date >= analysis_start,
    date <= analysis_end,
    policy == "ALS-4 Only",
    vax_scenario %in% c("Slower rate", "Faster rate")
  ) |>
  group_by(vax_scenario) |>
  summarise(
    inc_leq2dose = max_or_na(value[matrix == "cum_inc_leq2dose"]),
    inc_booster = max_or_na(value[matrix == "cum_inc_booster"]),
    hosp_total = max_or_na(value[matrix == "cum_hosp_total"]),
    death_total = max_or_na(value[matrix == "cum_death_total"]),
    .groups = "drop"
  )

print(severity_check)

ab_ymax <- incidence_df |>
  filter(
    vax_scenario == "Actual boosters",
    policy %in% c("No-ALS Only", "ALS-4 Only"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(policy, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

gh_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c(early_plot_label, "Actual boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

de_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Slower rate", "Faster rate"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[as.character(scenario_name)]]
  if (is.null(cols)) stop("No dose colours defined for scenario: ", scenario_name)
  cols
}

make_shaded_incidence_plot <- function(
    df, scenario_name, policy_name, panel_tag, panel_title,
    legend_position = c(0.50, 0.85), y_limit = NULL,
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
  
  y_upper <- if (is.null(y_limit)) max_or_na(plot_df$total_inc) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc, fill = "≤2 dose")) +
    geom_area(aes(y = inc_booster, fill = "Booster dose")) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 2, linetype = boundary_linetype) +
    geom_line(aes(y = inc_booster), colour = "black", linewidth = 2, linetype = boundary_linetype) +
    geom_vline(xintercept = release_date, linetype = "dashed", linewidth = 2) +
    scale_fill_manual(values = panel_cols, breaks = c("Booster dose", "≤2 dose")) +
    scale_x_date(
      limits = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0)
    ) +
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
      legend.text = element_text(size = 30),
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
  
  combo_cols <- c(
    "≤2 dose + No NPIs" = unname(panel_cols["≤2 dose"]),
    "Booster dose + No NPIs" = unname(panel_cols["Booster dose"]),
    "≤2 dose + Strict NPIs" = unname(panel_cols["≤2 dose"]),
    "Booster dose + Strict NPIs" = unname(panel_cols["Booster dose"])
  )
  
  combo_types <- c(
    "≤2 dose + No NPIs" = "solid",
    "Booster dose + No NPIs" = "solid",
    "≤2 dose + Strict NPIs" = "twodash",
    "Booster dose + Strict NPIs" = "twodash"
  )
  
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
      line_group = paste(dose_group, "+", policy_label),
      line_group = factor(line_group, levels = names(combo_cols))
    ) |>
    filter(!is.na(line_group))
  
  if (nrow(plot_df) == 0) {
    stop("Panel C data is empty after creating line_group.")
  }
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = line_group,
      linetype = line_group,
      group = line_group
    )
  ) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = combo_cols, name = NULL) +
    scale_linetype_manual(values = combo_types, name = NULL) +
    scale_x_date(
      limits = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0.08)
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 4))
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 30),
      legend.key.width = unit(2, "cm"),
      legend.background = element_blank(),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  
  slow_cols <- get_panel_cols("Slower rate")
  fast_cols <- get_panel_cols("Faster rate")
  
  combo_cols <- c(
    "Slower rate + ≤2 dose" = unname(slow_cols["≤2 dose"]),
    "Slower rate + Booster dose" = unname(slow_cols["Booster dose"]),
    "Faster rate + ≤2 dose" = unname(fast_cols["≤2 dose"]),
    "Faster rate + Booster dose" = unname(fast_cols["Booster dose"])
  )
  
  combo_types <- c(
    "Slower rate + ≤2 dose" = "solid",
    "Slower rate + Booster dose" = "solid",
    "Faster rate + ≤2 dose" = "twodash",
    "Faster rate + Booster dose" = "twodash"
  )
  
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
        "Slower rate" = "Slower rate",
        "Faster rate" = "Faster rate"
      ),
      line_group = paste(rollout_speed, "+", dose_group),
      line_group = factor(line_group, levels = names(combo_cols))
    ) |>
    filter(!is.na(line_group))
  
  if (nrow(plot_df) == 0) {
    stop("Panel F data is empty after creating line_group.")
  }
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = line_group,
      linetype = line_group,
      group = line_group
    )
  ) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = combo_cols, name = NULL) +
    scale_linetype_manual(values = combo_types, name = NULL) +
    scale_x_date(
      limits = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0.08)
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 4))
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 30),
      legend.background = element_blank(),
      legend.key.width = unit(2, "cm"),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

p_a <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "No-ALS Only",
  "(A)", "No NPIs", c(0.60, 0.85), ab_ymax, "solid"
)

p_b <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "ALS-4 Only",
  "(B)", "Strict NPIs", c(0.15, 0.85), ab_ymax, "twodash"
)

p_c <- make_cumulative_policy_plot(cumulative_df, "Actual boosters", "(C)")

p_d <- make_shaded_incidence_plot(
  incidence_df, "Slower rate", "ALS-4 Only",
  "(D)", "Slower booster doses", c(0.15, 0.85), de_ymax, "solid"
)

p_e <- make_shaded_incidence_plot(
  incidence_df, "Faster rate", "ALS-4 Only",
  "(E)", "Faster booster doses", c(0.15, 0.85), de_ymax, "twodash"
)

p_f <- make_cumulative_speed_plot(cumulative_df, "(F)")

p_g <- make_shaded_incidence_plot(
  incidence_df, early_plot_label, "ALS-4 Only",
  "(G)", "Early boosters", c(0.15, 0.85), gh_ymax, "twodash"
)

p_h <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "ALS-4 Only",
  "(H)", "Actual booster doses", c(0.15, 0.85), gh_ymax, "twodash"
)

plot_i_df <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Actual boosters", early_plot_label),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(
    peak_daily_incidence = max_or_na(value),
    .groups = "drop"
  ) |>
  mutate(
    dose_group = recode(
      matrix,
      "inc_leq2dose" = "≤2 dose",
      "inc_booster" = "Booster dose"
    ),
    scenario_label = if_else(
      vax_scenario == early_plot_label,
      "Earlier boosters",
      "Actual boosters"
    ),
    scenario_label = factor(
      scenario_label,
      levels = c("Earlier boosters", "Actual boosters")
    ),
    dose_group = factor(
      dose_group,
      levels = c("≤2 dose", "Booster dose")
    )
  )

p_i <- ggplot(plot_i_df, aes(x = dose_group, y = peak_daily_incidence, fill = scenario_label)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.68,
    colour = "#4D4D4D",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = scenario_cols[c("Earlier boosters", "Actual boosters")]) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Peak daily incidence", fill = NULL, tag = "(I)") +
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
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1, 1)) &
  theme(plot.margin = margin(8, 35, 8, 8))

save_png(
  file.path(fig_dir, "figure2.png"),
  figure1_3x3,
  width = 13500,
  height = 8000
)

severity_fill_cols <- c(
  "National rate" = "springgreen3",
  "Slower rate" = "magenta",
  "Actual boosters" = "#4D4D4D",
  "Earlier boosters" = "blue",
  "Faster rate" = "orange"
)

vax_short_labs <- c(
  "National rate" = "National\nrate",
  "Slower rate" = "Slower\nboosters",
  "Actual boosters" = "Actual\nboosters",
  "Earlier boosters" = "Earlier\nboosters",
  "Faster rate" = "Faster\nboosters"
)

figure2_df <- bind_rows(
  results_df |>
    filter(
      date >= analysis_start,
      date <= analysis_end,
      matrix == "cum_hosp_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max_or_na(value),
      .groups = "drop"
    ) |>
    mutate(outcome = "Hospitalizations"),
  
  results_df |>
    filter(
      date >= analysis_start,
      date <= analysis_end,
      matrix == "cum_death_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max_or_na(value),
      .groups = "drop"
    ) |>
    mutate(outcome = "Deaths")
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
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

if (nrow(figure2_df) == 0) {
  stop("figure2_df is empty. Check that cum_hosp_total and cum_death_total exist in results_df.")
}

x_order_df <- figure2_df |>
  filter(outcome == "Hospitalizations") |>
  group_by(policy_label) |>
  arrange(desc(severe_outcome), .by_group = TRUE) |>
  mutate(x_pos = row_number()) |>
  ungroup() |>
  select(policy_label, vax_scenario, x_pos)

figure2_df <- figure2_df |>
  left_join(x_order_df, by = c("policy_label", "vax_scenario")) |>
  mutate(
    x_lab = vax_short_labs[as.character(vax_scenario)],
    value_label = comma(round(severe_outcome, 0))
  ) |>
  filter(!is.na(x_pos))

figure2 <- ggplot(
  figure2_df,
  aes(x = x_pos, y = severe_outcome, fill = vax_scenario)
) +
  geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
  geom_text(aes(label = value_label), vjust = -0.35, size = 10) +
  facet_grid(outcome ~ policy_label, scales = "free_y") +
  scale_fill_manual(
    values = severity_fill_cols,
    breaks = vax_order,
    labels = c(
      "National rate" = "National rate",
      "Slower rate" = "Slower boosters",
      "Actual boosters" = "Actual boosters",
      "Earlier boosters" = "Earlier boosters",
      "Faster rate" = "Faster boosters"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:5,
    labels = NULL,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Severe outcomes",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35, colour = "black"),
    strip.text = element_text(size = 35),
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.45, "cm"),
    plot.title = element_text(size = 40, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm"),
    plot.margin = margin(t = 10, r = 20, b = 55, l = 20),
    plot.background = element_blank()
  ) +
  geom_text(
    data = figure2_df |> filter(outcome == "Deaths"),
    aes(y = -Inf, label = x_lab),
    vjust = 1.15,
    size = 8,
    lineheight = 0.85,
    inherit.aes = TRUE
  ) +
  coord_cartesian(clip = "off")

save_png(
  file.path(fig_dir, "figure3_same_init.png"),
  figure2,
  width = 10000,
  height = 5000
)