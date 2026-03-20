
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
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
cut_date       <- as.Date("2022-03-14")
fig_dir        <- "../figures"

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

# ============================================================
# INPUTS
# ============================================================
calibrator         <- rdsRead("calibrate")
fitserodata        <- rdsRead("fitsero.rds")
padded_beta_noals  <- rdsRead("policy_neutral")
padded_beta_base   <- rdsRead("base_beta")
padded_beta_als4   <- rdsRead("beta_als4")
booster_offset <- 107
sim_horizon <- 500

# V2_init <-  234,600 ## nl
# V2_init <-  277,707 ## natl

new_V2_init <- 277707 # use this nationa rate

day0 <- as.Date("2021-11-08") - lubridate::days(booster_offset)
date0 <- as.Date("2021-12-15")

# load weekly_scenarios_df
weekly_scenarios_df <- weekly_scenarios_df

generate_scenario_boosters <- function(scenario_name) {
  booster_inputs_df <- weekly_scenarios_df |>
    dplyr::filter(scenario == scenario_name) |>
    dplyr::select(week_end, scenario, weekly_rate) |>
    mutate(
      daydiff = as.numeric(week_end - lag(week_end, default = week_end[1] - 7)),
      daily_rate = weekly_rate / daydiff,
      days = as.numeric(week_end - day0)
    ) |>
    dplyr::filter(between(days, 0, 202 + booster_offset))
  
  changepoints <- booster_inputs_df$days - (booster_offset %% 7) - 1
  values <- booster_inputs_df$daily_rate
  
  if (!any(changepoints == 0)) {
    changepoints <- c(0, changepoints)
    values <- c(0, values)
  }
  
  return(list(changepoints = changepoints, vac_values = values))
}

# ============================================================
# Replace booster schedule in calibrated spec
# the model already uses:
#   booster_vac_values
#   booster_vac_changepoints
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
# Vaccination scenario specs
# ============================================================
early_boosters_results <- generate_scenario_boosters("Early Boosters (-35 Days)")
early_changepoints <- early_boosters_results$changepoints
early_vac_values   <- early_boosters_results$vac_values
spec_early <- set_booster_vax(early_boosters_results$changepoints, early_boosters_results$vac_values)

faster_rate_results <- generate_scenario_boosters("Faster Rate")
faster_changepoints <- faster_rate_results$changepoints
faster_vac_values   <- faster_rate_results$vac_values
spec_faster <- set_booster_vax(faster_rate_results$changepoints, faster_rate_results$vac_values)

slower_rate_results <- generate_scenario_boosters("Slower Rate")
slower_changepoints <- slower_rate_results$changepoints
slower_vac_values   <- slower_rate_results$vac_values
spec_slower <- set_booster_vax(slower_rate_results$changepoints, slower_rate_results$vac_values)

natl_rate_results <- generate_scenario_boosters("National Rate")
natl_changepoints <- natl_rate_results$changepoints
natl_vac_values   <- natl_rate_results$vac_values
spec_natl <- set_booster_vax(natl_rate_results$changepoints, natl_rate_results$vac_values)
spec_natl <- spec_natl |> 
  mp_tmb_insert(default = list(V2_init = new_V2_init))

# ============================================================
# Baseline calibrated model
# ============================================================
spec_fixed <- mp_optimized_spec(calibrator, "modified")
time_steps <- spec_fixed$default$time_steps

# ------------------------------------------------------------
# list of model specs with applied vaccination scenarios
# ------------------------------------------------------------
vax_scenario_specs <- list(
  "Actual Boosters"           = spec_fixed,
  "Slower Rate"          = spec_slower,
  "Faster Rate"          = spec_faster,
  "Early Boosters (-35 Days)" = spec_early,
  "National Rate"             = spec_natl
)

# ============================================================
# Apply ALS policies on the vaccination scenario specs
# ============================================================
update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

policy_list <- list(
  "Actual ALS"  = padded_beta_base,
  "ALS-4 Only"  = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)
policy_order <- names(policy_list)

# ============================================================
# COMMON VAX ORDER FIX
# ============================================================
vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
)

# ============================================================
# Simulate one scenario
# ============================================================
policy_names <- names(policy_list)
vax_names    <- names(vax_scenario_specs)

sim_grid <- expand.grid(
  policy_name = policy_names,
  vax_name    = vax_names,
  stringsAsFactors = FALSE
) |>
  mutate(combined_label = paste(policy_name, "+", vax_name))

n_sims <- nrow(sim_grid)
spec_list <- vector("list", n_sims)
all_sim_results <- vector("list", n_sims)

# prebuild specs
for (i in seq_len(n_sims)) {
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  spec_list[[i]] <-
    if (policy_name == "Actual ALS" && vax_name == "Actual Boosters") {
      spec_fixed
    } else {
      update_vax_spec_beta(
        vax_scenario_specs[[vax_name]],
        policy_list[[policy_name]]
      )
    }
}

# Simulate
for (i in seq_len(n_sims)) {
  all_sim_results[[i]] <- spec_list[[i]] |>
    mp_simulator(time_steps = sim_horizon, outputs = c("inc", "total_inc")) |>
    mp_trajectory() |>
    mutate(
      date = date0 + (time - booster_offset),
      policy = sim_grid$policy_name[i],
      vax_scenario = sim_grid$vax_name[i],
      combined_label = sim_grid$combined_label[i]
    )
}

final_results_df <- bind_rows(all_sim_results) |>
  mutate(policy = factor(policy, levels = c("No-ALS Only", "Actual ALS", "ALS-4 Only"))) |>
  filter(date >= analysis_start) |>
  mutate(
    vax_scenario = recode(
      vax_scenario,
      "National Rate" =  "National rate",
      "Slower Rate" = "Slower rate",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
      "Actual Boosters" = "Actual boosters",
      "Faster Rate" = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

# ============================================================
# Plotting
# ============================================================
vax_cols <- c(
  "National rate"             = "springgreen3",
  "Slower rate"               = "magenta",
  "Actual boosters"           = "black",
  "Early boosters (-35 Days)" = "blue",
  "Faster rate"               = "orange"
)

# ------------------------------------------------------------
# incidence plot for all policies
# ------------------------------------------------------------
p_inc_policy_facets <- final_results_df |>
  filter(matrix == "inc") |>
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ policy, ncol = 3, scales = "fixed") +
  scale_color_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  geom_vline(
    xintercept = as.Date("2022-03-14"),
    colour = "brown",
    linetype = 4,
    linewidth = 1.5
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Daily incidence by ALS policies and vaccination scenario",
    x = NULL,
    y = "Daily incidence",
    colour = "Vaccination scenario"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 15, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 25, color = "black"),
    strip.text = element_text(size = 25),
    strip.background = element_rect(fill = "gray"),
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
    plot.title = element_text(size = 35, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

print(p_inc_policy_facets)

png("../figures/incidence_all_policy.png", width = 6000, height = 2500, res = 300, bg = "white")
print(p_inc_policy_facets)
dev.off()

# ------------------------------------------------------------
# incidence plot for ALS-4 Only
# ------------------------------------------------------------
inc_als4_only <- final_results_df |>
  filter(matrix == "inc" & policy == "ALS-4 Only") |>
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  geom_vline(
    xintercept = as.Date("2022-03-14"),
    colour = "brown",
    linetype = 4,
    linewidth = 1.5
  ) +
  scale_y_continuous(labels = scales::comma) +
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
    strip.text = element_blank(),
    strip.background = element_blank(),
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
    plot.title = element_text(size = 20, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

print(inc_als4_only)

png("../figures/incidence_als4_only.png", width = 4500, height = 1800, res = 300, bg = "white")
print(inc_als4_only)
dev.off()

# ============================================================
# Plot stackbar plot
# ============================================================
start_date <- as.Date("2021-12-15")
cut_date   <- as.Date("2022-03-14")

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

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
  "Slower rate"               = "Slower rate",
  "National rate"             = "National rate",
  "Actual boosters"           = "Actual boosters",
  "Faster rate"               = "Faster rate",
  "Early boosters (-35 Days)" = "Early boosters (-35 Days)"
  
)

period_cols <- c(
  "ALS period (≤ Mar 14, 2022)"    = "blue",
  "No-ALS period (> Mar 14, 2022)" = "pink4"
)

totals_split_df <- final_results_df |>
  filter(
    matrix == "total_inc",
    date >= start_date
  ) |>
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
  aes(x = policy, y = infections, fill = period)
) +
  geom_col(width = 0.75) +
  facet_wrap(
    ~ vax_scenario,
    nrow = 1,
    drop = FALSE,
    labeller = labeller(vax_scenario = facet_labs)
  ) +
  scale_fill_manual(
    values = period_cols,
    breaks = c(
      "ALS period (≤ Mar 14, 2022)",
      "No-ALS period (> Mar 14, 2022)"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = policy_x_labs, drop = FALSE) +
  scale_y_continuous(labels = scales::comma) +
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
    plot.title = element_text(size = 24, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

print(p_totals_split)

png("../figures/vax_stackbar_als.png", width = 5000, height = 2500, res = 300, bg = "white")
print(p_totals_split)
dev.off()

# ============================================================
# Peak summary from daily incidence
# ============================================================
peak_df <- final_results_df |>
  filter(
    matrix == "inc",
    date >= analysis_start,
    date <= as.Date("2022-05-22")
  ) |>
  group_by(vax_scenario, policy) |>
  summarise(
    peak_inc = max(value, na.rm = TRUE),
    peak_date = date[which.max(value)],
    .groups = "drop"
  ) |>
  mutate(
    peak_lab = comma(round(peak_inc, 0)),
    peak_day = as.integer(peak_date - min(peak_date, na.rm = TRUE)),
    peak_date_lab = format(peak_date, "%b %d, %Y")
  )

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")
peak_df <- peak_df |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

vax_x_labs <- c(
  "Slower rate"               = "Slower\nboosters",
  "National rate"             = "National\nrate",
  "Actual boosters"           = "Actual\nboosters",
  "Early boosters (-35 Days)" = "Earlier boosters\n(-35 days)",
  "Faster rate"               = "Faster\nboosters"
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
    size = 8,
    position = position_dodge(width = 0.5),
    show.legend = FALSE
  ) +
  scale_colour_manual(values = policy_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.1, 0.18))) +
  scale_y_discrete(labels = rev(vax_x_labs)) +
  labs(
    title = "Peak daily infections by vaccination and policy scenario",
    x = "Peak daily infections",
    y = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    strip.text = element_text(size = 22),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 22),
    legend.background = element_blank(),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 25, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  ) +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_peak_point)

# ============================================================
# USER OPTION: optional threshold for policy communication
# set to NULL if not needed
# ============================================================
peak_threshold <- NULL
# peak_threshold <- 5000

# ============================================================
# Peak summary from daily incidence
# ============================================================
peak_df <- final_results_df |>
  filter(
    date >= analysis_start,
    date <= as.Date("2022-05-22")
  ) |>
  group_by(vax_scenario, policy) |>
  summarise(
    peak_inc = max(value[matrix == "inc"], na.rm = TRUE),
    peak_date = date[matrix == "inc"][which.max(value[matrix == "inc"])],
    final_total_inc = max(value[matrix == "total_inc"], na.rm = TRUE),
    cut_total_inc = max(value[matrix == "total_inc" & date <= as.Date("2022-03-14")], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    peak_lab = comma(round(peak_inc, 0)),
    peak_date_lab = format(peak_date, "%b %d %Y"),
    peak_date_num = as.numeric(peak_date),
    is_observed = policy == "Actual ALS" & vax_scenario == "Actual boosters"
  )

peak_df <- peak_df |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

vax_y_labs <- c(
  "Slower rate"               = "Slower rate",
  "National rate"             = "National rate",
  "Actual boosters"           = "Actual boosters",
  "Early boosters (-35 Days)" = "Earlier boosters\n(-35 days)",
  "Faster rate"               = "Faster rate"
)

heat_df <- peak_df |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = rev(vax_order))
  )

ref_peak <- heat_df |>
  filter(
    policy == "Actual ALS",
    vax_scenario == "Actual boosters"
  ) |>
  pull(peak_inc)

heat_df_rel <- heat_df |>
  mutate(
    pct_change_peak = 100 * (peak_inc - ref_peak) / ref_peak,
    pct_lab = sprintf("%.0f%%", pct_change_peak)
  )

p_heat_relative_enhanced <- ggplot(
  heat_df_rel,
  aes(x = policy, y = vax_scenario, fill = pct_change_peak)
) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = pct_lab),
    size = 8,
    fontface = "bold",
    vjust = -0.25
  ) +
  geom_text(
    aes(label = peak_date_lab),
    size = 8,
    vjust = 1.2
  ) +
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
    axis.text.x = element_text(size = 20, hjust = 0.5),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "plain"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.background = element_blank(),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.background = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_heat_relative_enhanced)

# ============================================================
# OPTIONAL combined figure
# ============================================================
p_peak_combo <- p_peak_point / p_heat_relative_enhanced +
  plot_layout(heights = c(1.15, 1))

png("../figures/peak_combo_enhanced.png", width = 5000, height = 5200, res = 300, bg = "white")
print(p_peak_combo)
dev.off()