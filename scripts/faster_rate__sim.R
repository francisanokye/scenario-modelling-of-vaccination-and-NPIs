
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

faster_rate_results <- generate_scenario_boosters("Faster Rate")
spec_faster <- set_booster_vax(
  faster_rate_results$changepoints,
  faster_rate_results$vac_values
)

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
  "Faster Rate"               = spec_faster
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
  "Actual boosters",
  "Faster rate"
)

vax_cols <- c(
  "Actual boosters"           = "black",
  "Faster rate"               = "orange"
)

scenario_labs <- c(
  "Actual Boosters"           = "Actual boosters",
  "Faster Rate"               = "Faster rate"
)

scenario_levels <- c(
  "Actual Boosters",
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
      "Actual boosters"           = "Actual Boosters",
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
  "Actual Boosters"           = "black",
  "Faster Rate"               = "orange"
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
