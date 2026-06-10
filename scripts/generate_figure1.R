
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

actual_booster_start <- as.Date("2021-11-08")
early_start_date     <- as.Date("2021-10-04")

NL_POP       <- 510550
CA_PRIMARY_ONLY_PROP_START <- 0.750
# new_V2_init <- round(0.750 * NL_POP)
new_V2_init <- NA_real_

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

booster_start_dates <- c(
  "Actual Boosters" = actual_booster_start,
  "Early Boosters"  = early_start_date,
  "Faster Boosters" = actual_booster_start,
  "Slower Boosters" = actual_booster_start,
  "National Rate"   = actual_booster_start
)

generate_scenario_boosters <- function(scenario_name, start_date) {
  
  scenario_daily <- booster_inputs_df |>
    filter(
      scenario == scenario_name,
      date     >= start_date
    ) |>
    arrange(day)
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster input rows found for scenario: ", scenario_name)
  }
  
  changepoints <- scenario_daily$day
  vac_values   <- scenario_daily$booster_vac_value
  
  if (!any(changepoints == 0)) {
    changepoints <- c(0L, changepoints)
    vac_values   <- c(0,  vac_values)
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

spec_early  <- generate_scenario_boosters("Early Boosters", early_start_date) |>
  with(set_booster_vax(changepoints, vac_values))

spec_faster <- generate_scenario_boosters("Faster Boosters", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values))

spec_slower <- generate_scenario_boosters("Slower Boosters", actual_booster_start) |>
  with(set_booster_vax(changepoints, vac_values))

calibrate_national_V2_init <- function(target_prop = CA_PRIMARY_ONLY_PROP_START) {
  
  natl_boosters <- generate_scenario_boosters("National Rate", actual_booster_start)
  
  get_prop <- function(v2_init_candidate) {
    spec_i <- set_booster_vax(
      natl_boosters$changepoints,
      natl_boosters$vac_values
    ) |>
      mp_tmb_insert(default = list(V2_init = as.numeric(v2_init_candidate)))
    
    traj <- spec_i |>
      mp_simulator(time_steps = sim_horizon, outputs = c("V2")) |>
      mp_trajectory() |>
      mutate(date = analysis_start + (time - offset))
    
    traj |>
      filter(date == actual_booster_start, matrix == "V2") |>
      summarise(value = first(as.numeric(value))) |>
      pull(value) / NL_POP
  }
  
  objective <- function(x) get_prop(x) - target_prop
  round(uniroot(objective, lower = 0, upper = NL_POP)$root)
}

new_V2_init <- calibrate_national_V2_init()

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
    policy       = factor(policy,       levels = policy_order),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order)
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

incidence_df  <- results_df |> filter(date >= analysis_start, date <= analysis_end)
cumulative_df <- results_df |> filter(date >= analysis_start, date <= analysis_end)

print(
  results_df |>
    filter(
      date   == analysis_start,
      matrix == "cum_booster",
      policy == "Actual ALS"
    ) |>
    select(vax_scenario, value) |>
    arrange(vax_scenario)
)

ab_ymax <- incidence_df |>
  filter(
    vax_scenario == "Actual Boosters",
    policy %in% c("No-ALS Only", "ALS-4 Only"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(policy, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

gh_ymax <- incidence_df |>
  filter(
    policy       == "ALS-4 Only",
    vax_scenario %in% c("Early Boosters", "Actual Boosters"),
    matrix       %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

de_ymax <- incidence_df |>
  filter(
    policy       == "ALS-4 Only",
    vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
    matrix       %in% c("inc_leq2dose", "inc_booster")
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
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    legend_position  = c(0.50, 0.85),
    y_limit          = NULL,
    boundary_linetype = "solid"
) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      policy       == policy_name,
      vax_scenario == scenario_name,
      matrix       %in% c("inc_leq2dose", "inc_booster")
    ) |>
    select(date, matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value, values_fill = 0) |>
    mutate(total_inc = inc_leq2dose + inc_booster)
  
  y_upper <- if (is.null(y_limit)) max_or_na(plot_df$total_inc) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc,  fill = "≤2 Dose")) +
    geom_area(aes(y = inc_booster, fill = "Booster Dose")) +
    geom_line(aes(y = total_inc),  colour = "black", linewidth = 2,
              linetype = boundary_linetype) +
    geom_line(aes(y = inc_booster), colour = "black", linewidth = 2,
              linetype = boundary_linetype) +
    geom_vline(xintercept = release_date, linetype = "dashed", linewidth = 2) +
    scale_fill_manual(values = panel_cols,
                      breaks = c("Booster Dose", "≤2 Dose")) +
    scale_x_date(
      limits       = c(analysis_start, analysis_end),
      date_breaks  = "1 month",
      date_labels  = "%d %b %Y",
      expand       = c(0.01, 0)
    ) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(title = panel_title, x = NULL, y = "Daily incidence",
         fill = NULL, tag = panel_tag) +
    theme_clean() +
    theme(
      axis.text.x      = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y      = element_text(size = 35),
      axis.title.y     = element_text(size = 35),
      legend.position  = legend_position,
      legend.justification = "center",
      legend.direction = "vertical",
      legend.text      = element_text(size = 30),
      legend.background = element_blank(),
      plot.title       = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.tag         = element_text(size = 45, face = "bold"),
      plot.background  = element_blank(),
      plot.margin      = margin(12, 35, 12, 12)
    )
}

make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  combo_cols <- c(
    "≤2 Dose + No NPIs"        = unname(panel_cols["≤2 Dose"]),
    "Booster Dose + No NPIs"   = unname(panel_cols["Booster Dose"]),
    "≤2 Dose + Strict NPIs"    = unname(panel_cols["≤2 Dose"]),
    "Booster Dose + Strict NPIs" = unname(panel_cols["Booster Dose"])
  )
  
  combo_types <- c(
    "≤2 Dose + No NPIs"        = "solid",
    "Booster Dose + No NPIs"   = "solid",
    "≤2 Dose + Strict NPIs"    = "twodash",
    "Booster Dose + Strict NPIs" = "twodash"
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
        "cum_inc_leq2dose" = "≤2 Dose",
        "cum_inc_booster"  = "Booster Dose"
      ),
      policy_label = recode(
        as.character(policy),
        "No-ALS Only" = "No NPIs",
        "ALS-4 Only"  = "Strict NPIs"
      ),
      line_group = paste(dose_group, "+", policy_label),
      line_group = factor(line_group, levels = names(combo_cols))
    ) |>
    filter(!is.na(line_group))
  
  if (nrow(plot_df) == 0) stop("Panel C data is empty after creating line_group.")
  
  ggplot(plot_df,
         aes(x = date, y = value,
             colour   = line_group,
             linetype = line_group,
             group    = line_group)) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = combo_cols, name = NULL) +
    scale_linetype_manual(values = combo_types, name = NULL) +
    scale_x_date(
      limits      = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand      = c(0.01, 0.08)
    ) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour   = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4))
    ) +
    theme_clean() +
    theme(
      axis.text.x      = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y      = element_text(size = 35),
      axis.title.y     = element_text(size = 35),
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.text      = element_text(size = 30),
      legend.key.width = unit(2, "cm"),
      legend.background = element_blank(),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin    = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag         = element_text(size = 45, face = "bold"),
      plot.background  = element_blank(),
      plot.margin      = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  
  slow_cols <- get_panel_cols("Slower Boosters")
  fast_cols <- get_panel_cols("Faster Boosters")
  
  combo_cols <- c(
    "Slower Boosters + ≤2 Dose"     = unname(slow_cols["≤2 Dose"]),
    "Slower Boosters + Booster Dose" = unname(slow_cols["Booster Dose"]),
    "Faster Boosters + ≤2 Dose"     = unname(fast_cols["≤2 Dose"]),
    "Faster Boosters + Booster Dose" = unname(fast_cols["Booster Dose"])
  )
  
  combo_types <- c(
    "Slower Boosters + ≤2 Dose"     = "solid",
    "Slower Boosters + Booster Dose" = "solid",
    "Faster Boosters + ≤2 Dose"     = "twodash",
    "Faster Boosters + Booster Dose" = "twodash"
  )
  
  plot_df <- df |>
    filter(
      policy       == "ALS-4 Only",
      vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
      matrix       %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "≤2 Dose",
        "cum_inc_booster"  = "Booster Dose"
      ),
      line_group = paste(vax_scenario, "+", dose_group),
      line_group = factor(line_group, levels = names(combo_cols))
    ) |>
    filter(!is.na(line_group))
  
  if (nrow(plot_df) == 0) stop("Panel F data is empty after creating line_group.")
  
  ggplot(plot_df,
         aes(x = date, y = value,
             colour   = line_group,
             linetype = line_group,
             group    = line_group)) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = combo_cols, name = NULL) +
    scale_linetype_manual(values = combo_types, name = NULL) +
    scale_x_date(
      limits      = c(analysis_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand      = c(0.01, 0.08)
    ) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour   = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4)),
      linetype = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 4))
    ) +
    theme_clean() +
    theme(
      axis.text.x      = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y      = element_text(size = 35),
      axis.title.y     = element_text(size = 35),
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.text      = element_text(size = 30),
      legend.background = element_blank(),
      legend.key.width = unit(2, "cm"),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin    = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag         = element_text(size = 45, face = "bold"),
      plot.background  = element_blank(),
      plot.margin      = margin(12, 55, 12, 12)
    )
}

p_a <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "No-ALS Only",
  "(A)", "No NPIs", c(0.60, 0.85), ab_ymax, "solid"
)

p_b <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "ALS-4 Only",
  "(B)", "Strict NPIs", c(0.15, 0.85), ab_ymax, "solid"
)

p_c <- make_cumulative_policy_plot(cumulative_df, "Actual Boosters", "(C)")

p_d <- make_shaded_incidence_plot(
  incidence_df, "Slower Boosters", "ALS-4 Only",
  "(D)", "Slower Boosters Doses", c(0.15, 0.85), de_ymax, "solid"
)

p_e <- make_shaded_incidence_plot(
  incidence_df, "Faster Boosters", "ALS-4 Only",
  "(E)", "Faster Boosters Doses", c(0.15, 0.85), de_ymax, "twodash"
)

p_f <- make_cumulative_speed_plot(cumulative_df, "(F)")

p_g <- make_shaded_incidence_plot(
  incidence_df, "Early Boosters", "ALS-4 Only",
  "(G)", "Early Boosters", c(0.15, 0.85), gh_ymax, "solid"
)

p_h <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "ALS-4 Only",
  "(H)", "Actual Boosters Doses", c(0.15, 0.85), gh_ymax, "twodash"
)

plot_i_df <- incidence_df |>
  filter(
    policy       == "ALS-4 Only",
    vax_scenario %in% c("Early Boosters", "Actual Boosters"),
    matrix       %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(
    peak_daily_incidence = max_or_na(value),
    .groups = "drop"
  ) |>
  mutate(
    dose_group = recode(
      matrix,
      "inc_leq2dose" = "≤2 Dose",
      "inc_booster"  = "Booster Dose"
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
    dose_group = factor(
      dose_group,
      levels = c("≤2 Dose", "Booster Dose")
    )
  )

p_i <- ggplot(plot_i_df,
              aes(x = dose_group, y = peak_daily_incidence,
                  fill = scenario_label)) +
  geom_col(
    position  = position_dodge(width = 0.72),
    width     = 0.68,
    colour    = "#4D4D4D",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = scenario_cols[c("Early Boosters", "Actual Boosters")]
  ) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Peak daily incidence", fill = NULL, tag = "(I)") +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 35, hjust = 0.5),
    axis.text.y      = element_text(size = 35),
    axis.title.y     = element_text(size = 35),
    legend.position  = c(0.80, 0.90),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text      = element_text(size = 35),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.tag         = element_text(size = 45, face = "bold"),
    plot.background  = element_blank(),
    plot.margin      = margin(12, 35, 12, 12)
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
  width  = 13500,
  height = 8000
)

severity_fill_cols <- scenario_cols

vax_short_labs <- c(
  "National Rate"   = "National\nRate",
  "Slower Boosters" = "Slower\nBoosters",
  "Actual Boosters" = "Actual\nBoosters",
  "Early Boosters"  = "Early\nBoosters",
  "Faster Boosters" = "Faster\nBoosters"
)

figure2_df <- bind_rows(
  results_df |>
    filter(
      date   >= analysis_start,
      date   <= analysis_end,
      matrix == "cum_hosp_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Hospitalizations"),
  
  results_df |>
    filter(
      date   >= analysis_start,
      date   <= analysis_end,
      matrix == "cum_death_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Deaths")
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
    policy_label = recode(
      as.character(policy),
      "No-ALS Only" = "No NPIs",
      "Actual ALS"  = "Actual NPIs",
      "ALS-4 Only"  = "Strict NPIs"
    ),
    policy_label = factor(
      policy_label,
      levels = c("No NPIs", "Actual NPIs", "Strict NPIs")
    ),
    outcome = factor(outcome, levels = c("Hospitalizations", "Deaths")),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order)
  )

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
    x_lab       = vax_short_labs[as.character(vax_scenario)],
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
    drop   = FALSE
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
    x     = NULL,
    y     = "Severe outcomes",
    fill  = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 30),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 35),
    axis.title.y     = element_text(size = 35, colour = "black"),
    strip.text       = element_text(size = 35),
    strip.background = element_blank(),
    legend.position  = "top",
    legend.text      = element_text(size = 25),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.45, "cm"),
    plot.title       = element_text(size = 40, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing    = unit(2, "cm"),
    plot.margin      = margin(t = 10, r = 20, b = 55, l = 20),
    plot.background  = element_blank()
  ) +
  geom_text(
    data = figure2_df |> filter(outcome == "Deaths"),
    aes(y = -Inf, label = x_lab),
    vjust      = 1.15,
    size       = 8,
    lineheight = 0.85,
    inherit.aes = TRUE
  ) +
  coord_cartesian(clip = "off")

save_png(
  file.path(fig_dir, "figure3.png"),
  figure2,
  width  = 10000,
  height = 5000
)
