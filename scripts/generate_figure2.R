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
offset         <- 150L
init_date      <- analysis_start - offset

release_date   <- as.Date("2022-03-17")
sim_horizon    <- 530L
analysis_end   <- init_date + sim_horizon - 1L

actual_booster_start <- as.Date("2021-11-07")
early_start_date     <- as.Date("2021-10-04")

new_V2_init <- 277707

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
  init_date + offset == analysis_start,
  length(padded_beta_noals) == length(padded_beta_base),
  length(padded_beta_als4)  == length(padded_beta_base),
  exists("booster_inputs_df")
)

scenario_cols <- c(
  "National Rate"   = "springgreen3",
  "Slower Boosters" = "magenta",
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "blue",
  "Faster Boosters" = "orange"
)

scenario_dose_cols <- list(
  "Actual Boosters" = c("≤2 Dose" = "grey70",   "Booster Dose" = "#4D4D4D"),
  "Slower Boosters" = c("≤2 Dose" = "#E8B6D2",  "Booster Dose" = "magenta"),
  "Faster Boosters" = c("≤2 Dose" = "#FDBB84",  "Booster Dose" = "orange"),
  "Early Boosters"  = c("≤2 Dose" = "#9ECAE1",  "Booster Dose" = "blue")
)

generate_scenario_boosters <- function(scenario_name, start_date) {
  
  scenario_daily <- booster_inputs_df |>
    filter(
      scenario == scenario_name,
      date >= start_date
    ) |>
    arrange(day)
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster input rows found for scenario: ", scenario_name)
  }
  
  changepoints <- scenario_daily$day
  vac_values   <- scenario_daily$booster_vac_value
  
  if (!any(changepoints == 0)) {
    changepoints <- c(0L, changepoints)
    vac_values   <- c(0, vac_values)
  }
  
  list(
    changepoints = as.integer(changepoints),
    vac_values   = as.numeric(vac_values)
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
  stopifnot(length(booster_days) == length(booster_vals))
  stopifnot(all(diff(booster_days) >= 0))
  
  spec <- mp_optimized_spec(calibrator, "modified")
  spec$integers$booster_vac_changepoints <- as.integer(booster_days)
  spec$default$booster_vac_values        <- as.numeric(booster_vals)
  spec
}

update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  
  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

spec_actual <- generate_scenario_boosters("Actual Boosters", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values))

spec_early <- generate_scenario_boosters("Early Boosters", early_start_date) |>
  with(set_booster_vax(changepoints, vac_values))

spec_faster <- generate_scenario_boosters("Faster Boosters", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values))

spec_slower <- generate_scenario_boosters("Slower Boosters", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values))

spec_natl <- generate_scenario_boosters("National Rate", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values)) |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))

vax_scenario_specs <- list(
  "National Rate"   = spec_natl,
  "Slower Boosters" = spec_slower,
  "Actual Boosters" = spec_actual,
  "Early Boosters"  = spec_early,
  "Faster Boosters" = spec_faster
)

policy_list <- list(
  "Actual ALS"  = padded_beta_base,
  "ALS-4 Only"  = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "National Rate",
  "Slower Boosters",
  "Actual Boosters",
  "Early Boosters",
  "Faster Boosters"
)

model_outputs <- c(
  "inc_symp_leq2dose",
  "inc_symp_booster",
  "inc_symp_total",
  "inc_asymp_leq2dose",
  "inc_asymp_booster",
  "inc_asymp_total",
  "cum_booster"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name    = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  spec_i <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[policy_name]]
  )
  
  all_sim_results[[i]] <- spec_i |>
    mp_simulator(
      time_steps = sim_horizon,
      outputs    = model_outputs
    ) |>
    mp_trajectory() |>
    mutate(
      date         = analysis_start + (time - offset),
      policy       = policy_name,
      vax_scenario = vax_name
    )
}

results_df <- bind_rows(all_sim_results) |>
  mutate(
    date         = as.Date(date),
    value        = as.numeric(value),
    matrix       = as.character(matrix),
    policy       = factor(policy, levels = policy_order),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order)
  )

saveVars(results_df)

required_outputs <- c(
  "inc_symp_leq2dose",
  "inc_symp_booster",
  "inc_asymp_leq2dose",
  "inc_asymp_booster"
)

missing_outputs <- setdiff(required_outputs, unique(results_df$matrix))

if (length(missing_outputs) > 0) {
  stop(
    "Missing model outputs: ",
    paste(missing_outputs, collapse = ", "),
    "\nRe-run upstream after updating flows."
  )
}

incidence_df <- results_df |>
  filter(date >= analysis_start, date <= analysis_end)

make_cumulative_outputs <- function(df, inc_leq2_name, inc_booster_name,
                                    cum_leq2_name, cum_booster_name) {
  
  df |>
    filter(matrix %in% c(inc_leq2_name, inc_booster_name)) |>
    group_by(policy, vax_scenario, matrix) |>
    arrange(date, .by_group = TRUE) |>
    mutate(
      value = cumsum(value),
      matrix = recode(
        matrix,
        !!inc_leq2_name := cum_leq2_name,
        !!inc_booster_name := cum_booster_name
      )
    ) |>
    ungroup()
}

symp_cumulative_df <- make_cumulative_outputs(
  incidence_df,
  "inc_symp_leq2dose",
  "inc_symp_booster",
  "cum_symp_leq2dose",
  "cum_symp_booster"
)

asymp_cumulative_df <- make_cumulative_outputs(
  incidence_df,
  "inc_asymp_leq2dose",
  "inc_asymp_booster",
  "cum_asymp_leq2dose",
  "cum_asymp_booster"
)

get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[as.character(scenario_name)]]
  if (is.null(cols)) stop("No dose colours defined for scenario: ", scenario_name)
  cols
}

make_ymax <- function(df, scenario_names, policy_names, leq2_output, booster_output) {
  df |>
    filter(
      vax_scenario %in% scenario_names,
      policy %in% policy_names,
      matrix %in% c(leq2_output, booster_output)
    ) |>
    group_by(policy, vax_scenario, date) |>
    summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
    summarise(ymax = max_or_na(total_inc)) |>
    pull(ymax)
}

make_shaded_incidence_plot <- function(
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    leq2_output,
    booster_output,
    y_label,
    legend_position = c(0.50, 0.85),
    y_limit = NULL,
    boundary_linetype = "solid"
) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      policy == policy_name,
      vax_scenario == scenario_name,
      matrix %in% c(leq2_output, booster_output)
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        !!leq2_output := "leq2",
        !!booster_output := "booster"
      )
    ) |>
    select(date, dose_group, value) |>
    pivot_wider(names_from = dose_group, values_from = value, values_fill = 0) |>
    mutate(total_inc = leq2 + booster)
  
  y_upper <- if (is.null(y_limit)) max_or_na(plot_df$total_inc) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc, fill = "≤2 Dose")) +
    geom_area(aes(y = booster, fill = "Booster Dose")) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 2,
              linetype = boundary_linetype) +
    geom_line(aes(y = booster), colour = "black", linewidth = 2,
              linetype = boundary_linetype) +
    geom_vline(xintercept = release_date, linetype = "dashed", linewidth = 2) +
    scale_fill_manual(
      values = panel_cols,
      breaks = c("Booster Dose", "≤2 Dose")
    ) +
    scale_x_date(
      limits = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(title = panel_title, x = NULL, y = y_label, fill = NULL, tag = panel_tag) +
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

make_cumulative_policy_plot <- function(
    df,
    scenario_name,
    panel_tag,
    cum_leq2_output,
    cum_booster_output,
    y_label
) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  combo_cols <- c(
    "≤2 Dose + No NPIs" = unname(panel_cols["≤2 Dose"]),
    "Booster Dose + No NPIs" = unname(panel_cols["Booster Dose"]),
    "≤2 Dose + Strict NPIs" = unname(panel_cols["≤2 Dose"]),
    "Booster Dose + Strict NPIs" = unname(panel_cols["Booster Dose"])
  )
  
  combo_types <- c(
    "≤2 Dose + No NPIs" = "solid",
    "Booster Dose + No NPIs" = "solid",
    "≤2 Dose + Strict NPIs" = "twodash",
    "Booster Dose + Strict NPIs" = "twodash"
  )
  
  plot_df <- df |>
    filter(
      vax_scenario == scenario_name,
      policy %in% c("No-ALS Only", "ALS-4 Only"),
      matrix %in% c(cum_leq2_output, cum_booster_output)
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        !!cum_leq2_output := "≤2 Dose",
        !!cum_booster_output := "Booster Dose"
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
  
  ggplot(
    plot_df,
    aes(x = date, y = value, colour = line_group,
        linetype = line_group, group = line_group)
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
    labs(x = NULL, y = y_label, tag = panel_tag) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE,
                            override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4))
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
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(
    df,
    panel_tag,
    cum_leq2_output,
    cum_booster_output,
    y_label
) {
  
  slow_cols <- get_panel_cols("Slower Boosters")
  fast_cols <- get_panel_cols("Faster Boosters")
  
  combo_cols <- c(
    "Slower Boosters + ≤2 Dose" = unname(slow_cols["≤2 Dose"]),
    "Slower Boosters + Booster Dose" = unname(slow_cols["Booster Dose"]),
    "Faster Boosters + ≤2 Dose" = unname(fast_cols["≤2 Dose"]),
    "Faster Boosters + Booster Dose" = unname(fast_cols["Booster Dose"])
  )
  
  combo_types <- c(
    "Slower Boosters + ≤2 Dose" = "solid",
    "Slower Boosters + Booster Dose" = "solid",
    "Faster Boosters + ≤2 Dose" = "twodash",
    "Faster Boosters + Booster Dose" = "twodash"
  )
  
  plot_df <- df |>
    filter(
      policy == "ALS-4 Only",
      vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
      matrix %in% c(cum_leq2_output, cum_booster_output)
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        !!cum_leq2_output := "≤2 Dose",
        !!cum_booster_output := "Booster Dose"
      ),
      line_group = paste(vax_scenario, "+", dose_group),
      line_group = factor(line_group, levels = names(combo_cols))
    ) |>
    filter(!is.na(line_group))
  
  ggplot(
    plot_df,
    aes(x = date, y = value, colour = line_group,
        linetype = line_group, group = line_group)
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
    labs(x = NULL, y = y_label, tag = panel_tag) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE,
                            override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4))
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
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

make_peak_bar_plot <- function(
    df,
    panel_tag,
    leq2_output,
    booster_output,
    y_label
) {
  
  plot_df <- df |>
    filter(
      policy == "ALS-4 Only",
      vax_scenario %in% c("Early Boosters", "Actual Boosters"),
      matrix %in% c(leq2_output, booster_output)
    ) |>
    group_by(vax_scenario, matrix) |>
    summarise(peak_daily_incidence = max_or_na(value), .groups = "drop") |>
    mutate(
      dose_group = recode(
        matrix,
        !!leq2_output := "≤2 Dose",
        !!booster_output := "Booster Dose"
      ),
      scenario_label = if_else(
        vax_scenario == "Early Boosters",
        "Early Boosters",
        "Actual Boosters"
      ),
      scenario_label = factor(
        scenario_label,
        levels = c("Early Boosters", "Actual Boosters")
      ),
      dose_group = factor(dose_group, levels = c("≤2 Dose", "Booster Dose"))
    )
  
  ggplot(
    plot_df,
    aes(x = dose_group, y = peak_daily_incidence, fill = scenario_label)
  ) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.68,
      colour = "#4D4D4D",
      linewidth = 0.25
    ) +
    scale_fill_manual(values = scenario_cols[c("Early Boosters", "Actual Boosters")]) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = y_label, fill = NULL, tag = panel_tag) +
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
}

make_figure_3x3 <- function(
    incidence_data,
    cumulative_data,
    file_name,
    leq2_output,
    booster_output,
    cum_leq2_output,
    cum_booster_output,
    inc_y_label,
    cum_y_label,
    peak_y_label,
    plot_title
) {
  
  ab_ymax <- make_ymax(
    incidence_data,
    "Actual Boosters",
    c("No-ALS Only", "ALS-4 Only"),
    leq2_output,
    booster_output
  )
  
  de_ymax <- make_ymax(
    incidence_data,
    c("Slower Boosters", "Faster Boosters"),
    "ALS-4 Only",
    leq2_output,
    booster_output
  )
  
  gh_ymax <- make_ymax(
    incidence_data,
    c("Early Boosters", "Actual Boosters"),
    "ALS-4 Only",
    leq2_output,
    booster_output
  )
  
  p_a <- make_shaded_incidence_plot(
    incidence_data, "Actual Boosters", "No-ALS Only",
    "(A)", "No NPIs", leq2_output, booster_output,
    inc_y_label, c(0.60, 0.85), ab_ymax, "solid"
  )
  
  p_b <- make_shaded_incidence_plot(
    incidence_data, "Actual Boosters", "ALS-4 Only",
    "(B)", "Strict NPIs", leq2_output, booster_output,
    inc_y_label, c(0.15, 0.85), ab_ymax, "solid"
  )
  
  p_c <- make_cumulative_policy_plot(
    cumulative_data, "Actual Boosters", "(C)",
    cum_leq2_output, cum_booster_output, cum_y_label
  )
  
  p_d <- make_shaded_incidence_plot(
    incidence_data, "Slower Boosters", "ALS-4 Only",
    "(D)", "Slower Boosters", leq2_output, booster_output,
    inc_y_label, c(0.15, 0.85), de_ymax, "solid"
  )
  
  p_e <- make_shaded_incidence_plot(
    incidence_data, "Faster Boosters", "ALS-4 Only",
    "(E)", "Faster Boosters", leq2_output, booster_output,
    inc_y_label, c(0.15, 0.85), de_ymax, "twodash"
  )
  
  p_f <- make_cumulative_speed_plot(
    cumulative_data, "(F)",
    cum_leq2_output, cum_booster_output, cum_y_label
  )
  
  p_g <- make_shaded_incidence_plot(
    incidence_data, "Early Boosters", "ALS-4 Only",
    "(G)", "Early Boosters", leq2_output, booster_output,
    inc_y_label, c(0.15, 0.85), gh_ymax, "solid"
  )
  
  p_h <- make_shaded_incidence_plot(
    incidence_data, "Actual Boosters", "ALS-4 Only",
    "(H)", "Actual Boosters", leq2_output, booster_output,
    inc_y_label, c(0.15, 0.85), gh_ymax, "twodash"
  )
  
  p_i <- make_peak_bar_plot(
    incidence_data, "(I)",
    leq2_output, booster_output, peak_y_label
  )
  
  fig <- (
    (p_a | p_b | p_c) /
      (p_d | p_e | p_f) /
      (p_g | p_h | p_i)
  ) +
    plot_layout(heights = c(1, 1, 1), widths = c(1, 1, 1)) +
    plot_annotation(title = plot_title) &
    theme(
      plot.title = element_text(
        size = 55,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20)
      ),
      plot.margin = margin(8, 35, 8, 8)
    )
  
  save_png(
    file.path(fig_dir, file_name),
    fig,
    width = 13500,
    height = 8000
  )
  
  fig
}

symptomatic_figure <- make_figure_3x3(
  incidence_data = incidence_df,
  cumulative_data = symp_cumulative_df,
  file_name = "figure_symptomatic_infections.png",
  leq2_output = "inc_symp_leq2dose",
  booster_output = "inc_symp_booster",
  cum_leq2_output = "cum_symp_leq2dose",
  cum_booster_output = "cum_symp_booster",
  inc_y_label = "Daily symptomatic infections",
  cum_y_label = "Cumulative symptomatic infections",
  peak_y_label = "Peak daily symptomatic infections",
  plot_title = "Symptomatic infections"
)

asymptomatic_figure <- make_figure_3x3(
  incidence_data = incidence_df,
  cumulative_data = asymp_cumulative_df,
  file_name = "figure_asymptomatic_infections.png",
  leq2_output = "inc_asymp_leq2dose",
  booster_output = "inc_asymp_booster",
  cum_leq2_output = "cum_asymp_leq2dose",
  cum_booster_output = "cum_asymp_booster",
  inc_y_label = "Daily asymptomatic infections",
  cum_y_label = "Cumulative asymptomatic infections",
  peak_y_label = "Peak daily asymptomatic infections",
  plot_title = "Asymptomatic infections"
)