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
booster_start  <- as.Date("2021-11-08")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
cut_date       <- as.Date("2022-03-14")
fig_dir        <- "../figures"
diag_dir  <- file.path(fig_dir, "diagnostics_")

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
calibrator         <- rdsRead("calibrate")
fitserodata        <- rdsRead("fitsero.rds")
padded_beta_noals  <- rdsRead("policy_neutral")
padded_beta_base   <- rdsRead("base_beta")
padded_beta_als4   <- rdsRead("beta_als4")
booster_offset <- 150#107
sim_horizon <- 530

# V2_init <-  234,600 ## nl   
new_V2_init <- 277707 # use this national rate

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
    dplyr::filter(between(days, 0, 159 + booster_offset))
  
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
  "Slower Rate"               = spec_slower,
  "Faster Rate"               = spec_faster,
  "Early Boosters"            = spec_early,
  "National Rate"             = spec_natl
)

# ============================================================
# Apply ALS policies on the vaccination scenario specs
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
# Simulate each policy x vaccination scenario
# Save diagnostic plot for each spec_list[[i]]
# Bind all results together
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

spec_list        <- vector("list", n_sims)
all_sim_results  <- vector("list", n_sims)
diag_plot_list   <- vector("list", n_sims)

# ------------------------------------------------------------
# build scenario-specific specs
# ------------------------------------------------------------
for (i in seq_len(n_sims)) {
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  spec_list[[i]] <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[policy_name]]
  )
}

# ------------------------------------------------------------
# simulate each spec, save diagnostic plot, store results
# ------------------------------------------------------------
for (i in seq_len(n_sims)) {
  
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  sim_dat <- spec_list[[i]] |>
    mp_simulator(
      sim_horizon,
      outputs = c("cum_booster", "total_inc", "inc")
    ) |>
    mp_trajectory() |>
    mutate(
      date = date0 + (time - booster_offset),
      policy = policy_name,
      vax_scenario = vax_name,
      combined_label = sim_grid$combined_label[i]
    )
  
  all_sim_results[[i]] <- sim_dat
  
  p_diag <- sim_dat |>
    ggplot(aes(x = date, y = value)) +
    geom_rect(
      data = data.frame(
        xmin = date0 - booster_offset,
        xmax = date0
      ),
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey30",
      alpha = 0.25
    ) +
    geom_vline(
      xintercept = date0 + (231 - booster_offset),
      color = "black",
      linetype = 4,
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = date0 + (300 - booster_offset),
      colour = "purple",
      linetype = 4,
      linewidth = 1.5
    ) +
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
      axis.text.x       = element_text(size = 10,  hjust = 1, angle = 45, margin = margin(t = 8)),
      axis.text.y       = element_text(size = 10),
      axis.title.x      = element_text(size = 10),
      axis.title.y      = element_text(size = 10),
      plot.title        = element_text(size = 12, hjust = 0.5, face = "plain"),
      strip.text        = element_text(size = 12, face = "plain", hjust = 0.5),
      strip.background  = element_blank(),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.spacing.x  = unit(1.2, "cm"),
      legend.spacing.y  = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.text       = element_text(size = 12, face = "plain"),
      panel.spacing     = unit(1.5, "cm"),
      plot.background   = element_blank()
    )
  
  diag_plot_list[[i]] <- p_diag
  
  file_stub <- paste0(
    "diag_",
    gsub("[^A-Za-z0-9]+", "_", policy_name),
    "__",
    gsub("[^A-Za-z0-9]+", "_", vax_name)
  )
  
  fn <- file.path(diag_dir, paste0(file_stub, ".png"))
  
  png(fn, width = 3500, height = 800, res = 300, bg = "white", type = "cairo")
  print(p_diag)
  dev.off()
}

# ------------------------------------------------------------
# bind all results together
# ------------------------------------------------------------
final_results_df <- bind_rows(all_sim_results) |>
  mutate(
    value = as.integer(value),
    policy = factor(policy, levels = c("No-ALS Only", "Actual ALS", "ALS-4 Only"))
  ) |>
  filter(date >= booster_start) |>
  mutate(
    vax_scenario = recode(
      vax_scenario,
      "National Rate"   = "National rate",
      "Slower Rate"     = "Slower rate",
      "Early Boosters"  = "Early boosters (-35 Days)",
      "Actual Boosters" = "Actual boosters",
      "Faster Rate"     = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

vax_cols <- c(
  "National rate"             = "springgreen3",
  "Slower rate"               = "magenta",
  "Actual boosters"           = "black",
  "Early boosters (-35 Days)" = "blue",
  "Faster rate"               = "orange"
)

booster_shot_facets <- final_results_df |>
  filter(matrix == "cum_booster" & policy == "Actual ALS") |>  
  ggplot(aes(x = date, y = 100 * (value / 510550), colour = vax_scenario, group = vax_scenario)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.1) +
  facet_wrap(~ policy, ncol = 3, scales = "fixed") +
  scale_color_manual(values = vax_cols, breaks = vax_order, drop = FALSE) +
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  geom_vline(
    xintercept = as.Date("2022-03-14"),
    colour = "brown",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = as.Date("2022-05-22"),
    colour = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  # Fixed the scales error here:
  scale_y_continuous(labels = scales::label_number(suffix = "%")) + 
  labs(
    title = "Cumulative Booster Uptake by Policy and Scenario",
    x = NULL,
    y = "Cumulative Coverage (%)",
    colour = NULL
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

print(booster_shot_facets)

png("../figures/booster_shots_policy.png", width = 6000, height = 2500, res = 300, bg = "white")
print(booster_shot_facets)
dev.off()


final_results_df <- final_results_df |>
  filter(date >= analysis_start)

# ------------------------------------------------------------
# incidence plot for all policies
# ------------------------------------------------------------
p_inc_policy_facets <- final_results_df |>
  filter(matrix == "inc") |>
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_line(linewidth = 1.5) +
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
  geom_vline(
    xintercept = as.Date("2022-05-22"),
    colour = "purple",
    linetype = "dashed",
    linewidth = 1.2
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
    plot.title = element_text(size = 40, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

print(p_inc_policy_facets)

png("../figures/incidence_all_policy.png", width = 5000, height = 3500, res = 300, bg = "white")
print(p_inc_policy_facets)
dev.off()

# ------------------------------------------------------------
# incidence plot for ALS-4 Only
# ------------------------------------------------------------
inc_als4_only <- final_results_df |>
  filter(matrix == "inc" & policy == "ALS-4 Only") |> 
  ggplot(aes(x = date, y = value, colour = vax_scenario, group = vax_scenario)) +
  geom_line(linewidth = 1.5) +
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
  geom_vline(
    xintercept = as.Date("2022-05-22"),
    colour = "purple",
    linetype = "dashed",
    linewidth = 1.2
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

png("../figures/incidence_als4_policy.png", width = 4500, height = 1800, res = 300, bg = "white")
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
  "Early boosters (-35 Days)" = "Earlier\nboosters\n(-35 days)",
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
    size = 10,
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
    axis.text.x = element_text(size = 30, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    strip.text = element_text(size = 30),
    strip.background = element_rect(fill = "grey90"),
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
    plot.title = element_text(size = 45, face = "plain", hjust = 0.5),
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
  "Early boosters (-35 Days)" = "Earlier\nboosters\n(-35 days)",
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
    size = 10,
    fontface = "bold",
    vjust = -0.25
  ) +
  geom_text(
    aes(label = peak_date_lab),
    size = 10,
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
    axis.text.x = element_text(size = 30, hjust = 0.5),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 45, hjust = 0.5, face = "plain"),
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
    plot.background = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# OPTIONAL combined figure
# ============================================================
p_peak_combo <- p_peak_point / p_heat_relative_enhanced +
  plot_layout(heights = c(1.15, 1))

png("../figures/peak_combo_enhanced.png", width = 7500, height = 5200, res = 300, bg = "white")
print(p_peak_combo)
dev.off()


# ============================================================
# Difference in cumulative infections:
# Actual boosters vs all other vaccination scenarios
# using final_results_df
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
  "No-ALS period (> Mar 14, 2022)"    = "darkgreen",
  "Total"         = "red"
)

# ------------------------------------------------------------
# Build cumulative totals by period from total_inc
# ------------------------------------------------------------
als_vax_raw <- final_results_df |>
  filter(
    matrix == "total_inc",
    date >= as.Date("2021-12-15")
  ) |>
  group_by(policy, vax_scenario) |>
  summarise(
    total_cut = max(value[date <= as.Date("2022-03-14")], na.rm = TRUE),
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

# ------------------------------------------------------------
# Keep Actual boosters + all comparison scenarios
# ------------------------------------------------------------
als_vax_clean <- als_vax_raw |>
  filter(vax_scenario %in% c(ref_vax, vax_compare_levels))

# ------------------------------------------------------------
# Reference values: Actual boosters
# ------------------------------------------------------------
ref_df <- als_vax_clean |>
  filter(vax_scenario == ref_vax) |>
  select(policy, period, ref_total = total_infections)

# ------------------------------------------------------------
# Comparison values: all non-reference scenarios
# ------------------------------------------------------------
comp_df <- als_vax_clean |>
  filter(vax_scenario != ref_vax) |>
  rename(comparison_vax = vax_scenario)

# ------------------------------------------------------------
# Compute differences by period
# diff = Actual boosters - comparison scenario
# ------------------------------------------------------------
period_diffs <- comp_df |>
  left_join(ref_df, by = c("policy", "period")) |>
  mutate(diff = ref_total - total_infections) |>
  select(policy, comparison_vax, period, diff)

# ------------------------------------------------------------
# Compute total difference across both periods
# ------------------------------------------------------------
total_diffs <- period_diffs |>
  group_by(policy, comparison_vax) |>
  summarise(diff = sum(diff, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "Total")

# ------------------------------------------------------------
# Final plotting data
# ------------------------------------------------------------
plot_df <- bind_rows(period_diffs, total_diffs) |>
  mutate(
    period = factor(
      period,
      levels = c("ALS period (≤ Mar 14, 2022)", "No-ALS period (> Mar 14, 2022)", "Total")
    ),
    comparison_vax = factor(
      comparison_vax,
      levels = vax_compare_levels
    ),
    policy = factor(
      policy,
      levels = als_levels
    )
  )

# ------------------------------------------------------------
# Optional nicer x labels
# ------------------------------------------------------------
vax_x_labs <- c(
  "National rate"             = "National\nrate",
  "Slower rate"               = "Slower\nrate",
  "Early boosters (-35 Days)" = "Early\nboosters\n(-35 Days)",
  "Faster rate"               = "Faster\nrate"
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p_actual_vs_all_vax <- ggplot(
  plot_df,
  aes(x = comparison_vax, y = diff, fill = period)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.6
  ) +
  facet_wrap(~ policy, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = period_palette, name = "Period") +
  scale_x_discrete(labels = vax_x_labs) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Vaccination scenarios",
    y = expression(Delta * " cumulative infections (Actual boosters - alternative scenario)"),
    title = "Difference in cumulative infections: Actual boosters vs all vaccination scenarios"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 28, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(size = 28, color = "black"),
    axis.title.y = element_text(size = 28, color = "black"),
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
    plot.title = element_text(size = 40, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm")
  )

print(p_actual_vs_all_vax)

png("../figures/actual_boosters_vs_all_vax_diffs.png",width = 9000,height = 4500,res = 300,bg = "white",type = "cairo")
print(p_actual_vs_all_vax)
dev.off()

