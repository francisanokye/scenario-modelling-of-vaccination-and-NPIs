# ============================================================
# Nature-style vaccination scenario figures with Rc(t)
# Single script: simulation + Rc + publication-style plots
# ============================================================

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
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

# ============================================================
# 1. CONFIG
# ============================================================
init_date      <- as.Date("2021-07-24")
booster_start  <- as.Date("2021-11-08")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
cut_date       <- as.Date("2022-03-14")
als4_start     <- as.Date("2022-01-08")
als4_end       <- as.Date("2022-02-07")
school_close   <- as.Date("2021-12-20")
school_open    <- as.Date("2022-01-26")

eligible_pop <- 426198
NL_POP       <- 510550
new_V2_init  <- 277707

fig_dir  <- "../figures"
diag_dir <- file.path(fig_dir, "diagnostics_vax_rc")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_png <- function(path, plot_obj, width = 5000, height = 3200,
                     res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

ensure_dir(fig_dir)
ensure_dir(diag_dir)

# ============================================================
# 2. INPUTS
# ============================================================
calibrator        <- rdsRead("calibrate")
fitserodata       <- rdsRead("fitsero.rds")
padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

stopifnot(exists("weekly_scenarios_df"))
stopifnot(exists("booster_inputs_df"))

# Use the same simulation window for every scenario.
sim_horizon <- as.integer(analysis_end - init_date)

# ============================================================
# 3. GLOBAL LABELS, ORDERS, AND COLOURS
# ============================================================
policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

policy_labs <- c(
  "No-ALS Only" = "No ALS",
  "Actual ALS"  = "Actual ALS",
  "ALS-4 Only"  = "ALS-4 only"
)

policy_cols <- c(
  "No-ALS Only" = "#D55E00",
  "Actual ALS"  = "#000000",
  "ALS-4 Only"  = "#0072B2"
)

scenario_levels <- c(
  "National Rate",
  "Slower Rate",
  "Actual Boosters",
  "Early Boosters (-35 Days)",
  "Faster Rate"
)

vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
)

scenario_labs <- c(
  "National Rate"             = "National rate",
  "Slower Rate"               = "Slower rate",
  "Actual Boosters"           = "Actual boosters",
  "Early Boosters (-35 Days)" = "Earlier boosters (-35 days)",
  "Faster Rate"               = "Faster rate"
)

vax_labs <- c(
  "National rate"             = "National rate",
  "Slower rate"               = "Slower rate",
  "Actual boosters"           = "Actual boosters",
  "Early boosters (-35 Days)" = "Earlier boosters\n(-35 days)",
  "Faster rate"               = "Faster rate"
)

vax_cols <- c(
  "National rate"             = "#009E73",
  "Slower rate"               = "#CC79A7",
  "Actual boosters"           = "#000000",
  "Early boosters (-35 Days)" = "#0072B2",
  "Faster rate"               = "#E69F00"
)

# ============================================================
# 4. THEME
# ============================================================
theme_nat <- theme_clean() +
  theme(
    axis.text.x       = element_text(size = 20, colour = "black"),
    axis.text.y       = element_text(size = 20, colour = "black"),
    axis.title.x      = element_text(size = 23, colour = "black"),
    axis.title.y      = element_text(size = 23, colour = "black"),
    plot.title        = element_text(size = 28, hjust = 0.5, colour = "black"),
    strip.text        = element_text(size = 22, colour = "black"),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 19),
    legend.background = element_blank(),
    plot.background   = element_blank()
  )

# ============================================================
# 5. BUILD VACCINATION SCENARIO SPECS
# ============================================================
generate_scenario_boosters <- function(scenario_name) {
  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name) |>
    arrange(day) |>
    transmute(
      changepoints = as.integer(day),
      vac_values   = as.numeric(booster_vac_value)
    )
  
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
  
  stopifnot(all(diff(scenario_daily$changepoints) >= 0))
  stopifnot(all(scenario_daily$vac_values >= 0))
  
  list(
    changepoints = scenario_daily$changepoints,
    vac_values   = scenario_daily$vac_values
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
  stopifnot(length(booster_days) == length(booster_vals))
  stopifnot(all(diff(booster_days) >= 0))
  stopifnot(all(booster_vals >= 0))
  
  spec <- mp_optimized_spec(calibrator, "modified")
  spec$integers$booster_vac_changepoints <- as.integer(booster_days)
  spec$default$booster_vac_values        <- as.numeric(booster_vals)
  spec
}

make_spec_for_scenario <- function(scenario_name) {
  x <- generate_scenario_boosters(scenario_name)
  set_booster_vax(x$changepoints, x$vac_values)
}

spec_actual <- make_spec_for_scenario("Actual Boosters")
spec_early  <- make_spec_for_scenario("Early Boosters (-35 Days)")
spec_faster <- make_spec_for_scenario("Faster Rate")
spec_slower <- make_spec_for_scenario("Slower Rate")
spec_natl   <- make_spec_for_scenario("National Rate") |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))

vax_scenario_specs <- list(
  "Actual Boosters"           = spec_actual,
  "Slower Rate"               = spec_slower,
  "Faster Rate"               = spec_faster,
  "Early Boosters (-35 Days)" = spec_early,
  "National Rate"             = spec_natl
)

# ============================================================
# 6. APPLY POLICY-SPECIFIC BETA(t)
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

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name    = names(vax_scenario_specs),
  stringsAsFactors = FALSE
) |>
  mutate(combined_label = paste(policy_name, "+", vax_name))

spec_list <- vector("list", nrow(sim_grid))
all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  spec_list[[i]] <- update_vax_spec_beta(
    vax_scenario_specs[[sim_grid$vax_name[i]]],
    policy_list[[sim_grid$policy_name[i]]]
  )
}

# ============================================================
# 7. SIMULATE ALL POLICY x VACCINATION SCENARIOS
# ============================================================
for (i in seq_len(nrow(sim_grid))) {
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  sim_dat <- spec_list[[i]] |>
    mp_simulator(
      sim_horizon,
      outputs = c("cum_booster", "total_inc", "inc")
    ) |>
    mp_trajectory() |>
    mutate(
      date           = init_date + time,
      policy         = policy_name,
      vax_scenario   = vax_name,
      combined_label = sim_grid$combined_label[i]
    )
  
  all_sim_results[[i]] <- sim_dat
}

final_results_all <- bind_rows(all_sim_results) |>
  mutate(
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      vax_scenario,
      "National Rate"             = "National rate",
      "Slower Rate"               = "Slower rate",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
      "Actual Boosters"           = "Actual boosters",
      "Faster Rate"               = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  ) |>
  filter(date >= booster_start, date <= analysis_end)

# Epidemiological outcomes should use the analysis window.
epi_results_df <- final_results_all |>
  filter(date >= analysis_start, date <= analysis_end)

# Check comparable windows across scenarios.
window_check <- epi_results_df |>
  group_by(policy, vax_scenario) |>
  summarise(
    start_date = min(date),
    end_date   = max(date),
    n_days     = n_distinct(date),
    .groups = "drop"
  )

print(window_check)
stopifnot(length(unique(window_check$start_date)) == 1)
stopifnot(length(unique(window_check$end_date)) == 1)
stopifnot(length(unique(window_check$n_days)) == 1)

# ============================================================
# 8. BOOSTER COVERAGE INPUTS AND MODEL OUTPUTS
# ============================================================
input_daily_df <- booster_inputs_df |>
  mutate(date = init_date + day) |>
  filter(date >= booster_start, date <= analysis_end) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  arrange(scenario, date) |>
  group_by(scenario) |>
  mutate(cum_val = cumsum(replace_na(booster_vac_value, 0))) |>
  ungroup()

input_cum_booster_df <- input_daily_df |>
  group_by(scenario, week_end) |>
  summarise(cum_val = max(cum_val, na.rm = TRUE), .groups = "drop") |>
  mutate(
    cum_prop = clip01(cum_val / eligible_pop),
    series = "Input curve"
  )

model_cum_booster_df <- final_results_all |>
  filter(matrix == "cum_booster", policy == "Actual ALS") |>
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
  summarise(cum_val = max(value, na.rm = TRUE), .groups = "drop") |>
  right_join(input_cum_booster_df |> select(scenario, week_end),
             by = c("scenario", "week_end")) |>
  group_by(scenario) |>
  arrange(week_end, .by_group = TRUE) |>
  fill(cum_val, .direction = "down") |>
  ungroup() |>
  mutate(
    cum_prop = clip01(cum_val / eligible_pop),
    series = "Model output"
  )

booster_compare_df <- bind_rows(
  input_cum_booster_df |> select(scenario, week_end, cum_val, cum_prop, series),
  model_cum_booster_df |> select(scenario, week_end, cum_val, cum_prop, series)
) |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    series = factor(series, levels = c("Input curve", "Model output"))
  )

# Daily booster coverage from model output for phase-plane plots.
booster_model_daily <- final_results_all |>
  filter(matrix == "cum_booster") |>
  transmute(
    date,
    policy,
    vax_scenario,
    booster_prop = clip01(value / eligible_pop)
  )

# ============================================================
# 9. Rc(t) CALCULATION
# ============================================================
# Rc(t) = beta(t) * [(mu/gamma_i) + ((1-mu)*zeta/gamma_a)] *
#         [p1*kappa1 + p2*kappa2 + p3*kappa3]
kappa1  <- 1
kappa2  <- 0.91
kappa3  <- 0.3
gamma_i <- 1 / 7
gamma_a <- 1 / 10
mu      <- 0.678
zeta    <- 0.75
p1      <- 0.15
p2      <- 0.85
p3      <- 0

bracket_term   <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)
susceptibility <- p1 * kappa1 + p2 * kappa2 + p3 * kappa3
mult_const     <- bracket_term * susceptibility

offset0 <- as.integer(analysis_start - init_date)

beta_vec_to_Rc_df <- function(beta_vec, policy_name) {
  tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number() - 1L,
           date = init_date + time,
           Rc = beta * mult_const,
           policy = policy_name) |>
    filter(date >= analysis_start, date <= analysis_end) |>
    transmute(date, policy, value = Rc)
}

Rc_policy_df <- bind_rows(
  beta_vec_to_Rc_df(padded_beta_base,  "Actual ALS"),
  beta_vec_to_Rc_df(padded_beta_noals, "No-ALS Only"),
  beta_vec_to_Rc_df(padded_beta_als4,  "ALS-4 Only")
) |>
  mutate(policy = factor(policy, levels = policy_order))

# Calibrated actual Rc(t) with uncertainty from beta_thing.
Rc_actual_ci_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  filter(matrix == "beta_thing") |>
  mutate(
    date = init_date + time,
    policy = "Actual ALS",
    value = value * mult_const,
    conf.low = conf.low * mult_const,
    conf.high = conf.high * mult_const
  ) |>
  filter(date >= analysis_start, date <= analysis_end) |>
  select(date, policy, value, conf.low, conf.high)

# Join Rc and booster coverage.
phase_df <- booster_model_daily |>
  filter(date >= analysis_start, date <= analysis_end) |>
  left_join(Rc_policy_df, by = c("date", "policy")) |>
  rename(Rc = value) |>
  filter(!is.na(Rc))

# ============================================================
# 10. SUMMARY METRICS
# ============================================================
inc_daily_df <- epi_results_df |>
  filter(matrix == "inc") |>
  transmute(date, policy, vax_scenario, inc = value)

total_daily_df <- epi_results_df |>
  filter(matrix == "total_inc") |>
  transmute(date, policy, vax_scenario, total_inc = value)

summary_df <- inc_daily_df |>
  group_by(policy, vax_scenario) |>
  summarise(
    peak_inc  = max(inc, na.rm = TRUE),
    peak_date = date[which.max(inc)],
    cum_inc   = sum(inc, na.rm = TRUE),
    attack_rate = cum_inc / NL_POP,
    .groups = "drop"
  ) |>
  left_join(
    booster_model_daily |>
      filter(date == analysis_end) |>
      select(policy, vax_scenario, booster_prop_end = booster_prop),
    by = c("policy", "vax_scenario")
  ) |>
  left_join(
    Rc_policy_df |>
      group_by(policy) |>
      summarise(
        first_Rc_below_1 = if (any(value < 1, na.rm = TRUE)) min(date[value < 1], na.rm = TRUE) else as.Date(NA),
        days_Rc_above_1  = sum(value > 1, na.rm = TRUE),
        days_Rc_below_1  = sum(value < 1, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "policy"
  ) |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

# Reference: Actual boosters + Actual ALS.
ref_summary <- summary_df |>
  filter(policy == "Actual ALS", vax_scenario == "Actual boosters")

ref_peak <- ref_summary$peak_inc[1]
ref_cum  <- ref_summary$cum_inc[1]
ref_date <- ref_summary$peak_date[1]

summary_df <- summary_df |>
  mutate(
    pct_change_peak = 100 * (peak_inc - ref_peak) / ref_peak,
    pct_change_cum  = 100 * (cum_inc - ref_cum) / ref_cum,
    days_peak_shift = as.numeric(peak_date - ref_date),
    peak_lab        = comma(round(peak_inc, 0)),
    cum_lab         = comma(round(cum_inc, 0)),
    attack_lab      = percent(attack_rate, accuracy = 0.1),
    peak_date_lab   = format(peak_date, "%b %d")
  )

# Relative reduction against No-ALS within each vaccination scenario.
noals_ref <- summary_df |>
  filter(policy == "No-ALS Only") |>
  select(vax_scenario, noals_cum = cum_inc, noals_peak = peak_inc)

relative_reduction_df <- summary_df |>
  left_join(noals_ref, by = "vax_scenario") |>
  mutate(
    cum_reduction_vs_noals  = 100 * (noals_cum - cum_inc) / noals_cum,
    peak_reduction_vs_noals = 100 * (noals_peak - peak_inc) / noals_peak
  )

# Difference curves relative to Actual boosters within each policy.
ref_curve <- inc_daily_df |>
  filter(vax_scenario == "Actual boosters") |>
  select(date, policy, ref_inc = inc)

diff_curve_df <- inc_daily_df |>
  left_join(ref_curve, by = c("date", "policy")) |>
  mutate(delta_inc = inc - ref_inc) |>
  filter(vax_scenario != "Actual boosters")

# ============================================================
# 11. FIGURE 1: MECHANISM PANEL
# ALS timeline + booster uptake + Rc(t)
# ============================================================
als_timeline_df <- tibble(
  xmin = c(analysis_start, als4_start, school_close, school_open),
  xmax = c(analysis_end, als4_end, school_open, analysis_end),
  y = c(4, 3, 2, 1),
  label = c("Analysis window", "ALS-4 period", "K-12 closed", "K-12 open")
)

p_timeline <- ggplot(als_timeline_df) +
  geom_segment(
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 2.3,
    lineend = "round"
  ) +
  geom_text(aes(x = xmin, y = y + 0.25, label = label), hjust = 0, size = 6) +
  scale_x_date(limits = c(analysis_start, analysis_end), date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(limits = c(0.5, 4.7)) +
  labs(x = NULL, y = NULL, title = "Intervention timeline") +
  theme_void() +
  theme(plot.title = element_text(size = 25, hjust = 0.5))

p_booster <- ggplot(input_cum_booster_df,
                    aes(x = week_end, y = cum_prop, colour = scenario)) +
  geom_line(linewidth = 1.7) +
  scale_colour_manual(values = setNames(vax_cols[scenario_labs], names(scenario_labs)),
                      labels = scenario_labs, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_date(limits = c(analysis_start, analysis_end), date_breaks = "1 month", date_labels = "%b") +
  labs(title = "Booster rollout scenarios", x = NULL, y = "Booster coverage") +
  theme_nat

p_rc <- ggplot() +
  annotate("rect", xmin = als4_start, xmax = als4_end, ymin = -Inf, ymax = Inf,
           alpha = 0.18, fill = "grey70") +
  geom_ribbon(
    data = Rc_actual_ci_df,
    aes(x = date, ymin = conf.low, ymax = conf.high),
    alpha = 0.25,
    fill = "black"
  ) +
  geom_line(data = Rc_policy_df,
            aes(x = date, y = value, colour = policy), linewidth = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  scale_colour_manual(values = policy_cols, labels = policy_labs, drop = FALSE) +
  scale_x_date(limits = c(analysis_start, analysis_end), date_breaks = "1 month", date_labels = "%b") +
  labs(title = expression("Control reproduction number, "*R[c]*"(t)"),
       x = NULL, y = expression(R[c])) +
  theme_nat

fig_mechanism <- p_timeline / p_booster / p_rc +
  plot_layout(heights = c(0.7, 1.2, 1.4))

save_png(file.path(fig_dir, "fig1_mechanism_timeline_boosters_Rc.png"),
         fig_mechanism, width = 5600, height = 5200)

# ============================================================
# 12. FIGURE 2: EPIDEMIC PHASE PLANE
# Booster coverage vs Rc(t)
# ============================================================
p_phase_rc <- ggplot(
  phase_df,
  aes(x = booster_prop, y = Rc, colour = policy, group = interaction(policy, vax_scenario))
) +
  geom_path(alpha = 0.75, linewidth = 1.2) +
  geom_point(data = phase_df |> filter(date %in% c(analysis_start, cut_date, analysis_end)),
             size = 2.2, alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  facet_wrap(~ vax_scenario, labeller = labeller(vax_scenario = vax_labs), ncol = 3) +
  scale_colour_manual(values = policy_cols, labels = policy_labs, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = expression("Control landscape: booster coverage versus "*R[c]*"(t)"),
    x = "Booster coverage",
    y = expression(R[c])
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig2_phase_plane_booster_coverage_Rc.png"),
         p_phase_rc, width = 6200, height = 4200)

# ============================================================
# 13. FIGURE 3: PEAK TIMING VS PEAK MAGNITUDE
# ============================================================
p_peak_tradeoff <- ggplot(
  summary_df,
  aes(x = peak_date, y = peak_inc, colour = policy, shape = vax_scenario)
) +
  geom_point(size = 5.5) +
  geom_text(aes(label = peak_date_lab), vjust = -1.0, size = 5.2, show.legend = FALSE) +
  scale_colour_manual(values = policy_cols, labels = policy_labs, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    title = "Peak timing versus peak daily infections",
    x = "Peak date",
    y = "Peak daily infections"
  ) +
  theme_nat +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_png(file.path(fig_dir, "fig3_peak_timing_vs_peak_magnitude.png"),
         p_peak_tradeoff, width = 5600, height = 3600)

# ============================================================
# 14. FIGURE 4: RELATIVE REDUCTION AGAINST NO-ALS
# ============================================================
p_rel_reduction <- ggplot(
  relative_reduction_df |> filter(policy != "No-ALS Only"),
  aes(x = vax_scenario, y = cum_reduction_vs_noals, colour = policy, group = policy)
) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_line(linewidth = 1.2, position = position_dodge(width = 0.25)) +
  geom_point(size = 4.8, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = policy_cols, labels = policy_labs, drop = FALSE) +
  scale_x_discrete(labels = vax_labs) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Relative reduction in cumulative infections compared with No ALS",
    x = NULL,
    y = "Reduction in cumulative infections"
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig4_relative_reduction_vs_noals.png"),
         p_rel_reduction, width = 5600, height = 3600)

# ============================================================
# 15. FIGURE 5: HEATMAP OF VACCINATION x ALS INTERACTION
# ============================================================
p_heat_cum <- ggplot(
  summary_df,
  aes(x = policy, y = vax_scenario, fill = attack_rate)
) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = attack_lab), size = 6.5, fontface = "bold") +
  scale_x_discrete(labels = policy_labs, drop = FALSE) +
  scale_y_discrete(labels = vax_labs, drop = FALSE) +
  scale_fill_gradient(low = "white", high = "firebrick", labels = percent_format(accuracy = 1)) +
  labs(
    title = "Vaccination–policy interaction: cumulative infection burden",
    x = NULL,
    y = NULL,
    fill = "Attack rate"
  ) +
  theme_nat +
  theme(panel.grid = element_blank())

save_png(file.path(fig_dir, "fig5_heatmap_vaccination_policy_attack_rate.png"),
         p_heat_cum, width = 5600, height = 3600)

# ============================================================
# 16. FIGURE 6: TIME BOUGHT BY POLICY AND VACCINATION
# Days peak shifted compared with Actual ALS + Actual boosters
# ============================================================
p_time_bought <- ggplot(
  summary_df,
  aes(x = booster_prop_end, y = days_peak_shift, colour = policy, label = vax_scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(group = policy), linewidth = 1.2) +
  geom_point(size = 5.2) +
  scale_colour_manual(values = policy_cols, labels = policy_labs, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Time bought by vaccination and policy scenarios",
    x = "Booster coverage at end of analysis window",
    y = "Peak shift relative to actual scenario, days"
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig6_time_bought_peak_shift.png"),
         p_time_bought, width = 5600, height = 3600)

# ============================================================
# 17. FIGURE 7: COUNTERFACTUAL DIFFERENCE CURVES
# Difference in daily infections vs Actual boosters, within each policy
# ============================================================
p_diff_curve <- ggplot(
  diff_curve_df,
  aes(x = date, y = delta_inc, colour = vax_scenario)
) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ policy, labeller = labeller(policy = policy_labs), ncol = 1, scales = "free_y") +
  scale_colour_manual(values = vax_cols, labels = vax_labs, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "Daily infections gained or averted compared with actual boosters",
    x = NULL,
    y = expression(Delta*" daily infections")
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig7_counterfactual_difference_curves.png"),
         p_diff_curve, width = 5600, height = 5200)

# ============================================================
# 18. FIGURE 8: CUMULATIVE INFECTION FRACTION OVER TIME
# ============================================================
cum_frac_df <- inc_daily_df |>
  group_by(policy, vax_scenario) |>
  arrange(date, .by_group = TRUE) |>
  mutate(cum_frac = cumsum(inc) / NL_POP) |>
  ungroup()

p_cum_frac <- ggplot(
  cum_frac_df,
  aes(x = date, y = cum_frac, colour = vax_scenario)
) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ policy, labeller = labeller(policy = policy_labs), ncol = 1) +
  scale_colour_manual(values = vax_cols, labels = vax_labs, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "Cumulative infection fraction across vaccination scenarios",
    x = NULL,
    y = "Cumulative infections / population"
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig8_cumulative_infection_fraction.png"),
         p_cum_frac, width = 5600, height = 5200)

# ============================================================
# 19. FIGURE 9: Rc THRESHOLD CROSSING SUMMARY
# ============================================================
threshold_df <- Rc_policy_df |>
  group_by(policy) |>
  summarise(
    first_Rc_below_1 = if (any(value < 1, na.rm = TRUE)) min(date[value < 1], na.rm = TRUE) else as.Date(NA),
    days_Rc_above_1  = sum(value > 1, na.rm = TRUE),
    days_Rc_below_1  = sum(value < 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(policy = factor(policy, levels = policy_order))

p_threshold <- ggplot(threshold_df, aes(x = policy, y = days_Rc_above_1, fill = policy)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = days_Rc_above_1), vjust = -0.4, size = 7) +
  scale_x_discrete(labels = policy_labs, drop = FALSE) +
  scale_fill_manual(values = policy_cols, guide = "none") +
  labs(
    title = expression("Days with "*R[c]*"(t) > 1"),
    x = NULL,
    y = "Number of days"
  ) +
  theme_nat

save_png(file.path(fig_dir, "fig9_Rc_threshold_crossing_days.png"),
         p_threshold, width = 4200, height = 3400)

# ============================================================
# 20. FIGURE 10: MAIN MULTI-PANEL SUMMARY
# Recommended for manuscript main text
# ============================================================
fig_main_summary <- (p_booster / p_rc) | (p_peak_tradeoff / p_heat_cum) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))

save_png(file.path(fig_dir, "fig10_main_summary_vax_policy_story.png"),
         fig_main_summary, width = 8200, height = 6200)

# ============================================================
# 21. SAVE DATA OBJECTS FOR TABLES OR SUPPLEMENT
# ============================================================
saveRDS(final_results_all, file.path(fig_dir, "vax_policy_simulation_results_all.rds"))
saveRDS(epi_results_df,   file.path(fig_dir, "vax_policy_epi_results_analysis_window.rds"))
saveRDS(Rc_policy_df,    file.path(fig_dir, "Rc_policy_df.rds"))
saveRDS(summary_df,      file.path(fig_dir, "summary_metrics_vax_policy.rds"))

write.csv(summary_df, file.path(fig_dir, "summary_metrics_vax_policy.csv"), row.names = FALSE)
write.csv(window_check, file.path(fig_dir, "simulation_window_check.csv"), row.names = FALSE)

message("Done. All Nature-style vaccination and Rc figures saved in: ", normalizePath(fig_dir))
