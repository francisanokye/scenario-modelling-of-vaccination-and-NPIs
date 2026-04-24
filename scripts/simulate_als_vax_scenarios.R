
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

# ============================================================
# CONFIG
# ============================================================
init_date      <- as.Date("2021-07-24")
booster_start  <- as.Date("2021-11-08")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
cut_date       <- as.Date("2022-03-14")
fig_dir        <- "../figures"
diag_dir       <- file.path(fig_dir, "diagnostics_")

eligible_pop   <- 426198
NL_POP         <- 510550

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

ensure_dir(fig_dir)
ensure_dir(diag_dir)

# ============================================================
# INPUTS
# ============================================================
calibrator        <- rdsRead("calibrate")
fitserodata       <- rdsRead("fitsero.rds")
padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

stopifnot(exists("weekly_scenarios_df"))
stopifnot(exists("booster_inputs_df"))

sim_horizon <- 530
new_V2_init <- 277707

# ============================================================
# Use daily booster inputs
# ============================================================
generate_scenario_boosters <- function(scenario_name) {

  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name) |>
    arrange(day) |>
    transmute(changepoints = as.integer(day), vac_values   = as.numeric(booster_vac_value))

  if (nrow(scenario_daily) == 0) {
    stop(paste("No booster_inputs_df rows found for scenario:", scenario_name))
  }

  scenario_daily <- scenario_daily |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)

  if (!any(scenario_daily$changepoints == 0)) {
    scenario_daily <- bind_rows(
      tibble(changepoints = 0L, vac_values = 0),
      scenario_daily
    ) |>
      arrange(changepoints)
  }

  list(changepoints = scenario_daily$changepoints, vac_values   = scenario_daily$vac_values)
  
  }

# ============================================================
# replace booster schedule in calibrated spec
# ============================================================
set_booster_vax <- function(booster_days, booster_vals) {
  stopifnot(length(booster_days) == length(booster_vals))
  stopifnot(all(diff(booster_days) >= 0))

  spec <- mp_optimized_spec(calibrator, "modified")
  spec$integers$booster_vac_changepoints <- as.integer(booster_days)
  spec$default$booster_vac_values        <- as.numeric(booster_vals)
  spec
}

# ============================================================
# vaccination scenario specs
# ============================================================
actual_boosters_results <- generate_scenario_boosters("Actual Boosters")
spec_actual <- set_booster_vax(
  actual_boosters_results$changepoints,
  actual_boosters_results$vac_values
)

early_boosters_results <- generate_scenario_boosters("Early Boosters (-35 Days)")
spec_early <- set_booster_vax(
  early_boosters_results$changepoints,
  early_boosters_results$vac_values
)

faster_rate_results <- generate_scenario_boosters("Faster Rate")
spec_faster <- set_booster_vax(
  faster_rate_results$changepoints,
  faster_rate_results$vac_values
)

slower_rate_results <- generate_scenario_boosters("Slower Rate")
spec_slower <- set_booster_vax(
  slower_rate_results$changepoints,
  slower_rate_results$vac_values
)

natl_rate_results <- generate_scenario_boosters("National Rate")
spec_natl <- set_booster_vax(
  natl_rate_results$changepoints,
  natl_rate_results$vac_values
)
spec_natl <- spec_natl |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))

# ============================================================
# baseline calibrated model
# ============================================================
spec_fixed <- mp_optimized_spec(calibrator, "modified")
time_steps <- spec_fixed$default$time_steps

# ============================================================
# list of model specs with applied vaccination scenarios
# ============================================================
vax_scenario_specs <- list(
  "Actual Boosters"           = spec_actual,
  "Slower Rate"               = spec_slower,
  "Faster Rate"               = spec_faster,
  "Early Boosters (-35 Days)" = spec_early,
  "National Rate"             = spec_natl
)

# ============================================================
# apply ALS policies on the vaccination scenario specs
# ============================================================
update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)

  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

policy_list <- list(
  "Actual ALS"  = padded_beta_base,
  "ALS-4 Only"  = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

# ============================================================
# COMMON VAX ORDER
# ============================================================
vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
)

vax_cols <- c(
  "National rate"             = "springgreen3",
  "Slower rate"               = "magenta",
  "Actual boosters"           = "black",
  "Early boosters (-35 Days)" = "blue",
  "Faster rate"               = "orange"
)

scenario_labs <- c(
  "National Rate"             = "National rate",
  "Slower Rate"               = "Slower rate",
  "Actual Boosters"           = "Actual boosters",
  "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
  "Faster Rate"               = "Faster rate"
)

scenario_levels <- c(
  "National Rate",
  "Slower Rate",
  "Actual Boosters",
  "Early Boosters (-35 Days)",
  "Faster Rate"
)

# ============================================================
# Simulate each policy x vaccination scenario
# ============================================================
policy_names <- names(policy_list)
vax_names    <- names(vax_scenario_specs)

sim_grid <- expand.grid(
  policy_name = policy_names,
  vax_name    = vax_names,
  stringsAsFactors = FALSE
) |>
  mutate(
    combined_label = paste(policy_name, "+", vax_name)
  )

n_sims <- nrow(sim_grid)

spec_list       <- vector("list", n_sims)
all_sim_results <- vector("list", n_sims)
diag_plot_list  <- vector("list", n_sims)

for (i in seq_len(n_sims)) {
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]

  spec_list[[i]] <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[policy_name]]
  )
}

for (i in seq_len(n_sims)) {

  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]

  sim_dat <- spec_list[[i]] |>
    mp_simulator(
      sim_horizon,
      outputs = c("cum_booster", "total_inc", "inc") # mp_state_vars(spec_list[[i]]), 
    ) |>
    mp_trajectory() |>
    mutate(
      date = init_date + time,
      policy = policy_name,
      vax_scenario = vax_name,
      combined_label = sim_grid$combined_label[i]
    )

  all_sim_results[[i]] <- sim_dat

  p_diag <- sim_dat |>
    ggplot(aes(x = date, y = value)) +
    geom_line() +
    labs(
      title = paste("Diagnostic compartment trajectories:", policy_name, "-", vax_name),
      x = NULL,
      y = NULL
    ) +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "month",
      date_labels = "%b %d %Y"
    ) +
    facet_wrap(~matrix, scales = "free") +
    theme_clean() +
    theme(
      axis.text.x       = element_text(size = 10, hjust = 1, angle = 45, margin = margin(t = 8)),
      axis.text.y       = element_text(size = 10),
      axis.title.x      = element_text(size = 10),
      axis.title.y      = element_text(size = 10),
      plot.title        = element_text(size = 12, hjust = 0.5),
      strip.text        = element_text(size = 12),
      strip.background  = element_blank(),
      panel.spacing     = unit(1.5, "cm")
    )

  diag_plot_list[[i]] <- p_diag

  file_stub <- paste0(
    "diag_",
    gsub("[^A-Za-z0-9]+", "_", policy_name),
    "__",
    gsub("[^A-Za-z0-9]+", "_", vax_name)
  )

  save_png(
    file.path(diag_dir, paste0(file_stub, ".png")),
    p_diag,
    width = 3500,
    height = 800
  )
}

# ============================================================
# Bind all results together
# ============================================================
final_results_df <- bind_rows(all_sim_results) |>
  mutate(
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order)
  ) |>
  filter(date >= booster_start) |>
  mutate(
    vax_scenario = recode(
      vax_scenario,
      "National Rate"             = "National rate",
      "Slower Rate"               = "Slower rate",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
      "Actual Boosters"           = "Actual boosters",
      "Faster Rate"               = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

# ============================================================
# VACCINATION-SCENARIO COMPARISON
# Use booster_inputs_df and same init_date clock
# ============================================================
input_daily_df <- booster_inputs_df |>
  mutate(
    date = init_date + day
  ) |>
  filter(date >= booster_start, date <= analysis_end) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  arrange(scenario, date) |>
  group_by(scenario) |>
  mutate(
    cum_val = cumsum(replace_na(booster_vac_value, 0))
  ) |>
  ungroup()

input_cum_booster_df <- input_daily_df |>
  group_by(scenario, week_end) |>
  summarise(
    cum_val = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cum_prop = cum_val / eligible_pop,
    series = "Input curve"
  )

model_daily_df <- final_results_df |>
  filter(
    matrix == "cum_booster",
    policy == "Actual ALS",
    date >= booster_start,
    date <= analysis_end
  ) |>
  mutate(
    scenario = recode(
      as.character(vax_scenario),
      "National rate"             = "National Rate",
      "Slower rate"               = "Slower Rate",
      "Actual boosters"           = "Actual Boosters",
      "Early boosters (-35 Days)" = "Early Boosters (-35 Days)",
      "Faster rate"               = "Faster Rate"
    ),
    scenario = factor(scenario, levels = scenario_levels),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(scenario, week_end) |>
  summarise(
    cum_val = max(value, na.rm = TRUE),
    .groups = "drop"
  )

model_cum_booster_df <- input_cum_booster_df |>
  select(scenario, week_end) |>
  left_join(model_daily_df, by = c("scenario", "week_end")) |>
  group_by(scenario) |>
  arrange(week_end, .by_group = TRUE) |>
  tidyr::fill(cum_val, .direction = "down") |>
  ungroup() |>
  mutate(
    cum_prop = cum_val / eligible_pop,
    series = "Model output"
  )

booster_compare_df <- bind_rows(
  input_cum_booster_df |>
    select(scenario, week_end, cum_val, cum_prop, series),
  model_cum_booster_df |>
    select(scenario, week_end, cum_val, cum_prop, series)
) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    series = factor(series, levels = c("Input curve", "Model output"))
  )

vax_cols_panel <- c(
  "Slower Rate"               = "magenta",
  "Early Boosters (-35 Days)" = "blue",
  "Actual Boosters"           = "black",
  "Faster Rate"               = "orange",
  "National Rate"             = "springgreen3"
)

booster_shot_facets <- ggplot(
  booster_compare_df,
  aes(
    x = week_end,
    y = cum_prop,
    linetype = series,
    group = series,
    colour = scenario
  )
) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ scenario,
    ncol = 3,
    scales = "fixed",
    labeller = labeller(scenario = scenario_labs)
  ) +
  geom_vline(
    xintercept = booster_start,
    colour = "brown",
    linetype = 2,
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = vax_cols_panel,
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c(
      "Input curve"  = "solid",
      "Model output" = "dashed"
    )
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0.02, 0)
  ) +
  labs(
    title = "Cumulative booster uptake: input curve vs model output",
    x = NULL,
    y = "Proportion of NL eligible population",
    linetype = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 25, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    strip.text = element_text(size = 35),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 40, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.spacing = unit(1.0, "cm")
  ) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.99)
  )

save_png(
  "../figures/booster_shots_policy.png",
  booster_shot_facets,
  width = 9000,
  height = 4500
)

# ============================================================
# Incidence ribbon summary
# ============================================================
inc_ribbon_df <- final_results_df |>
  filter(matrix == "inc") |>
  group_by(date, policy, vax_scenario) |>
  summarise(
    mean  = mean(value, na.rm = TRUE),
    se    = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) |>
  mutate(
    lower = ifelse(is.na(lower), mean, lower),
    upper = ifelse(is.na(upper), mean, upper)
  ) |>
  filter(is.finite(lower), is.finite(upper), upper >= lower)

# ============================================================
# Incidence plot for all policies
# ============================================================
p_inc_policy_facets <- final_results_df |>
  filter(matrix == "inc") |>
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_ribbon(
    data = inc_ribbon_df,
    aes(x = date, ymin = lower, ymax = upper, fill = vax_scenario, group = vax_scenario),
    inherit.aes = FALSE,
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~ policy, ncol = 3, scales = "fixed") +
  scale_color_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_fill_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0.01, 0)
  ) +
  geom_vline(
    xintercept = cut_date,
    colour = "brown",
    linetype = 4,
    linewidth = 1.5
  ) +
  geom_vline(
    xintercept = analysis_end,
    colour = "purple",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Daily incidence by ALS policies and vaccination scenario",
    x = NULL,
    y = "Daily incidence",
    colour = "Vaccination scenario"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 15, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 28),
    axis.title.y = element_text(size = 28, color = "black"),
    strip.text = element_text(size = 28),
    strip.background = element_rect(fill = "gray"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 40, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

save_png(
  "../figures/incidence_all_policy.png",
  p_inc_policy_facets,
  width = 7500,
  height = 2500
)

# ============================================================
# Incidence plot for ALS-4 Only
# ============================================================
inc_ribbon_als4 <- inc_ribbon_df |>
  filter(policy == "ALS-4 Only")

inc_als4_only <- final_results_df |>
  filter(matrix == "inc", policy == "ALS-4 Only") |>
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_ribbon(
    data = inc_ribbon_als4,
    aes(x = date, ymin = lower, ymax = upper, fill = vax_scenario, group = vax_scenario),
    inherit.aes = FALSE,
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_fill_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  geom_vline(
    xintercept = cut_date,
    colour = "brown",
    linetype = 4,
    linewidth = 1.5
  ) +
  geom_vline(
    xintercept = analysis_end,
    colour = "purple",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Daily incidence under ALS-4 Only across vaccination scenario",
    x = NULL,
    y = "Daily incidence",
    colour = "Vaccination scenario"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 14, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18, color = "black"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

save_png(
  "../figures/incidence_als4_policy.png",
  inc_als4_only,
  width = 4500,
  height = 1800
)

# ============================================================
# Total infections split by period
# ============================================================
policy_x_labs <- c(
  "No-ALS Only" = "No-ALS\nOnly",
  "Actual ALS"  = "Actual\nALS",
  "ALS-4 Only"  = "ALS-4\nOnly"
)

vax_order_fixed <- c(
  "Slower rate",
  "Actual boosters",
  "National rate",
  "Early boosters (-35 Days)",
  "Faster rate"
)

vax_x_labs <- c(
  "Slower rate"               = "Slower\nrate",
  "National rate"             = "National\nrate",
  "Actual boosters"           = "Actual\nboosters",
  "Early boosters (-35 Days)" = "Early\nboosters\n(-35 Days)",
  "Faster rate"               = "Faster\nrate"
)

period_cols <- c(
  "ALS period (≤ Mar 14, 2022)"    = "blue",
  "No-ALS period (> Mar 14, 2022)" = "goldenrod"
)

totals_split_df <- final_results_df |>
  filter(matrix == "total_inc", date >= analysis_start) |>
  group_by(policy, vax_scenario) |>
  summarise(
    total_cut = max(value[date <= cut_date], na.rm = TRUE),
    total_end = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    total_cut = ifelse(is.finite(total_cut), total_cut, 0),
    total_end = ifelse(is.finite(total_end), total_end, 0),
    infections_before_cut = pmax(total_cut, 0),
    infections_after_cut  = pmax(total_end - total_cut, 0),
    policy = factor(as.character(policy), levels = policy_order),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order_fixed)
  )

totals_split_long <- totals_split_df |>
  pivot_longer(
    cols = c(infections_before_cut, infections_after_cut),
    names_to = "period",
    values_to = "infections"
  ) |>
  mutate(
    period = factor(
      period,
      levels = c("infections_after_cut", "infections_before_cut"),
      labels = c("No-ALS period (> Mar 14, 2022)", "ALS period (≤ Mar 14, 2022)")
    ),
    policy = factor(as.character(policy), levels = policy_order),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order_fixed)
  )

facet_labs <- function(x) {
  out <- unname(vax_x_labs[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

p_totals_split <- ggplot(
  totals_split_long,
  aes(x = vax_scenario, y = infections, fill = period)
) +
  geom_col(width = 0.75) +
  facet_wrap(
    ~ policy,
    nrow = 1,
    drop = FALSE,
    labeller = labeller(vax_scenario = facet_labs)
  ) +
  scale_fill_manual(
    values = period_cols,
    breaks = c("ALS period (≤ Mar 14, 2022)", "No-ALS period (> Mar 14, 2022)"),
    name = NULL
  ) +
  scale_x_discrete(labels = vax_x_labs, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Total infections split by ALS policy period and after ALS removal",
    x = NULL,
    y = "Total infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 15, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20, color = "black"),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 24, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

save_png(
  "../figures/vax_stackbar_als.png",
  p_totals_split,
  width = 6500,
  height = 2500
)

# ============================================================
# Peak figures
# ============================================================
peak_df <- final_results_df |>
  filter(matrix == "inc", date <= analysis_end) |>
  group_by(vax_scenario, policy) |>
  summarise(
    peak_inc = max(value, na.rm = TRUE),
    peak_date = date[which.max(value)],
    .groups = "drop"
  ) |>
  mutate(
    peak_lab = comma(round(peak_inc, 0)),
    peak_date_lab = format(peak_date, "%b %d, %Y")
  ) |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

vax_peak_labs <- c(
  "Slower rate"               = "Slower\nboosters",
  "Actual boosters"           = "Actual\nboosters",
  "National rate"             = "National\nrate",
  "Faster rate"               = "Faster\nboosters",
  "Early boosters (-35 Days)" = "Earlier\nboosters\n(-35 days)"
)

policy_cols <- c(
  "No-ALS Only" = "red",
  "Actual ALS"  = "orange4",
  "ALS-4 Only"  = "darkgreen"
)

p_peak_point <- ggplot(
  peak_df,
  aes(x = peak_inc, y = vax_scenario, colour = policy)
) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  geom_text(
    aes(label = peak_lab),
    hjust = -0.20,
    size = 9,
    position = position_dodge(width = 0.5),
    show.legend = FALSE
  ) +
  scale_colour_manual(values = policy_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.1, 0.18))) +
  scale_y_discrete(labels = rev(vax_peak_labs)) +
  labs(
    title = "Peak daily infections by vaccination and policy scenario",
    x = "Peak daily infections",
    y = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30, color = "black"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.background = element_blank(),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 45, hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

peak_df2 <- final_results_df |>
  filter(date >= booster_start, date <= analysis_end) |>
  group_by(vax_scenario, policy) |>
  summarise(
    peak_inc = max(value[matrix == "inc"], na.rm = TRUE),
    peak_date = date[matrix == "inc"][which.max(value[matrix == "inc"])],
    .groups = "drop"
  ) |>
  mutate(
    peak_date_lab = format(peak_date, "%b %d %Y")
  ) |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

heat_df <- peak_df2 |>
  mutate(vax_scenario = factor(vax_scenario, levels = rev(vax_order)))

ref_peak <- heat_df |>
  filter(policy == "Actual ALS", vax_scenario == "Actual boosters") |>
  pull(peak_inc)

heat_df_rel <- heat_df |>
  mutate(
    pct_change_peak = 100 * (peak_inc - ref_peak) / ref_peak,
    pct_lab = sprintf("%.0f%%", pct_change_peak)
  )

vax_y_labs <- c(
  "Slower rate"               = "Slower rate",
  "Actual boosters"           = "Actual boosters",
  "National rate"             = "National rate",
  "Faster rate"               = "Faster rate",
  "Early boosters (-35 Days)" = "Earlier\nboosters\n(-35 days)"
)

p_heat_relative_enhanced <- ggplot(
  heat_df_rel,
  aes(x = policy, y = vax_scenario, fill = pct_change_peak)
) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = pct_lab), size = 10, fontface = "bold", vjust = -0.25) +
  geom_text(aes(label = peak_date_lab), size = 10, vjust = 1.2) +
  geom_tile(
    data = heat_df_rel |> filter(policy == "Actual ALS", vax_scenario == "Actual boosters"),
    fill = NA,
    colour = "black",
    linewidth = 1.5
  ) +
  scale_fill_gradient2(
    low = "#00FFFF",
    mid = "lightyellow",
    high = "red",
    midpoint = 0,
    labels = function(x) paste0(round(x), "%"),
    name = "% change in peak\nvs observed reference"
  ) +
  scale_y_discrete(labels = vax_y_labs, drop = FALSE) +
  labs(
    title = "Relative change in peak daily infections and timing across scenarios",
    x = NULL,
    y = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 45, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.background = element_blank(),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    panel.grid = element_blank()
  ) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

p_peak_combo <- p_peak_point / p_heat_relative_enhanced +
  plot_layout(heights = c(1.15, 1))

save_png(
  "../figures/peak_combo_enhanced.png",
  p_peak_combo,
  width = 7500,
  height = 5200
)

# ============================================================
# Difference in cumulative infections
# ============================================================
ref_vax <- "Actual boosters"
als_levels <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")
vax_compare_levels <- c(
  "National rate",
  "Slower rate",
  "Early boosters (-35 Days)",
  "Faster rate"
)

period_palette <- c(
  "ALS period (≤ Mar 14, 2022)" = "navy",
  "No-ALS period (> Mar 14, 2022)" = "darkgreen",
  "Total" = "red"
)

als_vax_raw <- final_results_df |>
  filter(matrix == "total_inc", date >= booster_start) |>
  group_by(policy, vax_scenario) |>
  summarise(
    total_cut = max(value[date <= cut_date], na.rm = TRUE),
    total_end = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    total_cut = ifelse(is.finite(total_cut), total_cut, 0),
    total_end = ifelse(is.finite(total_end), total_end, 0),
    infections_dec15_mar14 = pmax(total_cut, 0),
    infections_mar15_end   = pmax(total_end - total_cut, 0)
  ) |>
  select(policy, vax_scenario, infections_dec15_mar14, infections_mar15_end) |>
  pivot_longer(
    cols = c(infections_dec15_mar14, infections_mar15_end),
    names_to = "period",
    values_to = "total_infections"
  ) |>
  mutate(
    period = recode(
      period,
      "infections_dec15_mar14" = "ALS period (≤ Mar 14, 2022)",
      "infections_mar15_end"   = "No-ALS period (> Mar 14, 2022)"
    )
  )

als_vax_clean <- als_vax_raw |>
  filter(vax_scenario %in% c(ref_vax, vax_compare_levels))

ref_df <- als_vax_clean |>
  filter(vax_scenario == ref_vax) |>
  select(policy, period, ref_total = total_infections)

comp_df <- als_vax_clean |>
  filter(vax_scenario != ref_vax) |>
  rename(comparison_vax = vax_scenario)

period_diffs <- comp_df |>
  left_join(ref_df, by = c("policy", "period")) |>
  mutate(diff = ref_total - total_infections) |>
  select(policy, comparison_vax, period, diff)

total_diffs <- period_diffs |>
  group_by(policy, comparison_vax) |>
  summarise(diff = sum(diff, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "Total")

plot_df <- bind_rows(period_diffs, total_diffs) |>
  mutate(
    period = factor(period, levels = c("ALS period (≤ Mar 14, 2022)", "No-ALS period (> Mar 14, 2022)", "Total")),
    comparison_vax = factor(comparison_vax, levels = vax_compare_levels),
    policy = factor(policy, levels = als_levels)
  )

vax_x_labs2 <- c(
  "National rate"             = "National\nrate",
  "Slower rate"               = "Slower\nrate",
  "Early boosters (-35 Days)" = "Early\nboosters\n(-35 Days)",
  "Faster rate"               = "Faster\nrate"
)

p_actual_vs_all_vax <- ggplot(
  plot_df,
  aes(x = comparison_vax, y = diff, fill = period)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.6) +
  facet_wrap(~ policy, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = period_palette, name = "Period") +
  scale_x_discrete(labels = vax_x_labs2) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Vaccination scenarios",
    y = expression(Delta * " cumulative infections"),
    title = "Change in cumulative infections: Actual boosters vs all vaccination scenarios"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 28, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(size = 28, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    strip.text = element_text(size = 30),
    strip.background = element_rect(fill = "gray"),
    legend.position = "bottom",
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
    panel.spacing = unit(2, "cm")
  )

save_png(
  "../figures/actual_boosters_vs_all_vax_diffs.png",
  p_actual_vs_all_vax,
  width = 9000,
  height = 3500
)


