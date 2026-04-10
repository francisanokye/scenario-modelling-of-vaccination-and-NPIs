
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
calibrator    <- rdsRead("calibrate")
fitserodata   <- rdsRead("fitsero.rds")

# ALS-4 duration beta series
padded_beta_ALS4_14days  <- rdsRead("beta_ALS4_14days")
padded_beta_ALS4_28days  <- rdsRead("beta_ALS4_28days")
padded_beta_ALS4_40days  <- rdsRead("beta_ALS4_40days")
padded_beta_ALS4_60days  <- rdsRead("beta_ALS4_60days")
padded_beta_ALS4_90days  <- rdsRead("beta_ALS4_90days")
padded_beta_ALS4_120days <- rdsRead("beta_ALS4_120days")
padded_beta_ALS4_150days <- rdsRead("beta_ALS4_150days")
padded_beta_ALS4_180days <- rdsRead("beta_ALS4_180days")

booster_offset <- 150
sim_horizon    <- 500
day0           <- as.Date("2021-11-08") - lubridate::days(booster_offset)
date0          <- as.Date("2021-12-15")
new_V2_init <- 277707

# weekly_scenarios_df should already exist from your vaccination-scenario script
weekly_scenarios_df <- weekly_scenarios_df

# ============================================================
# HELPERS
# ============================================================
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

# ============================================================
# BASE SPEC + VACCINATION SCENARIO SPECS
# IMPORTANT:
# use exact source names from weekly_scenarios_df here
# ============================================================
spec_fixed <- mp_optimized_spec(calibrator, "modified")
time_steps <- spec_fixed$default$time_steps

# source scenario names from weekly_scenarios_df
slower_res <- generate_scenario_boosters("Slower Rate")
faster_res <- generate_scenario_boosters("Faster Rate")
early_res  <- generate_scenario_boosters("Early Boosters (-35 Days)")
natl_res   <- generate_scenario_boosters("National Rate")

spec_slower <- set_booster_vax(slower_res$changepoints, slower_res$vac_values)
spec_faster <- set_booster_vax(faster_res$changepoints, faster_res$vac_values)
spec_early  <- set_booster_vax(early_res$changepoints,  early_res$vac_values)
spec_natl   <- set_booster_vax(natl_res$changepoints,   natl_res$vac_values)
spec_natl <- spec_natl |> 
  mp_tmb_insert(default = list(V2_init = new_V2_init))

# display order used everywhere
vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
)

# display-name spec list
vax_scenario_specs <- list(
  "Actual boosters"           = spec_fixed,
  "Slower rate"               = spec_slower,
  "Faster rate"               = spec_faster,
  "Early boosters (-35 Days)" = spec_early,
  "National rate"             = spec_natl
)

# ============================================================
# POLICY SCENARIOS: ALS-4 DURATION
# ============================================================
policy_list <- list(
  "14"  = padded_beta_ALS4_14days,
  "28"  = padded_beta_ALS4_28days,
  "40"  = padded_beta_ALS4_40days,
  "60"  = padded_beta_ALS4_60days,
  "90"  = padded_beta_ALS4_90days,
  "120" = padded_beta_ALS4_120days,
  "150" = padded_beta_ALS4_150days,
  "180" = padded_beta_ALS4_180days
)

policy_days_order <- c("14", "28", "40", "60", "90", "120", "150", "180")

# ============================================================
# SIMULATE ALL COMBINATIONS
# ============================================================
policy_names <- names(policy_list)
vax_names    <- names(vax_scenario_specs)

sim_grid <- expand.grid(
  als4_days = policy_names,
  vax_name  = vax_names,
  stringsAsFactors = FALSE
)

n_sims <- nrow(sim_grid)
spec_list <- vector("list", n_sims)
all_sim_results <- vector("list", n_sims)

for (i in seq_len(n_sims)) {
  als4_days <- sim_grid$als4_days[i]
  vax_name  <- sim_grid$vax_name[i]
  
  spec_list[[i]] <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[als4_days]]
  )
}

for (i in seq_len(n_sims)) {
  all_sim_results[[i]] <- spec_list[[i]] |>
    mp_simulator(time_steps = sim_horizon, outputs = c("inc", "total_inc")) |>
    mp_trajectory() |>
    pivot_wider(names_from = matrix, values_from = value) |>
    mutate(
      date = date0 + (time - booster_offset),
      als4_days = sim_grid$als4_days[i],
      vax_scenario = sim_grid$vax_name[i]
    )
}

final_results_df <- bind_rows(all_sim_results) |>
  filter(date >= analysis_start, date <= analysis_end) |>
  mutate(
    als4_days = factor(als4_days, levels = policy_days_order),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

# ============================================================
# MAXIMUM PEAK DAILY INFECTIONS
# one value per vaccination scenario × ALS-4 duration
# ============================================================
peak_df <- final_results_df |>
  group_by(vax_scenario, als4_days) |>
  summarise(
    peak_inc = max(inc, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    als4_days = factor(als4_days, levels = policy_days_order),
    als4_days_num = as.numeric(as.character(als4_days)),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

# ============================================================
# COLOURS
# ============================================================
vax_cols <- c(
  "National rate"             = "springgreen3",
  "Slower rate"               = "magenta",
  "Actual boosters"           = "black",
  "Early boosters (-35 Days)" = "blue",
  "Faster rate"               = "orange"
)

# ============================================================
# THRESHOLD + CROSSING DAYS
# label only crossing days
# ============================================================
threshold <- 1500

cross_df <- peak_df |>
  arrange(vax_scenario, als4_days_num) |>
  group_by(vax_scenario) |>
  group_modify(~{
    df <- .x
    x <- df$als4_days_num
    y <- df$peak_inc
    
    idx <- which(y[-length(y)] > threshold & y[-1] <= threshold)[1]
    
    if (is.na(idx)) {
      return(tibble(
        cross_day = NA_real_,
        cross_peak = NA_real_
      ))
    }
    
    x1 <- x[idx]
    x2 <- x[idx + 1]
    y1 <- y[idx]
    y2 <- y[idx + 1]
    
    cross_day <- x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
    
    tibble(
      cross_day = cross_day,
      cross_peak = threshold
    )
  }) |>
  ungroup() |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = vax_order),
    cross_lab = ifelse(is.na(cross_day), NA_character_, paste0("~", round(cross_day)))
  )

# ============================================================
# PLOT
# ============================================================
p_peak_lines <- ggplot(
  peak_df,
  aes(
    x = als4_days_num,
    y = peak_inc,
    color = vax_scenario,
    group = vax_scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 5) +
  geom_hline(
    yintercept = threshold,
    color = "red",
    linetype = 2.5,
    linewidth = 3
  ) +
  geom_point(
    data = cross_df |> filter(!is.na(cross_day)),
    aes(x = cross_day, y = cross_peak, color = vax_scenario),
    size = 6,
    shape = 8,
    fill = "gold",
    stroke = 2,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = cross_df |> filter(!is.na(cross_day)),
    aes(x = cross_day, xend = cross_day, y = 0, yend = cross_peak, color = vax_scenario),
    linetype = 4,
    linewidth = 1.5,
    alpha = 1,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = cross_df |> filter(!is.na(cross_day)),
    aes(x = cross_day, y = cross_peak, label = cross_lab, color = vax_scenario),
    vjust = 5,
    hjust = 0.5,
    size = 9,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = vax_cols,
    breaks = vax_order,
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(14, 28, 40, 60, 90, 120, 150, 180),
    labels = c("14", "28", "40", "60", "90", "120", "150", "180")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Maximum daily infections across ALS-4 durations and vaccination scenarios",
    x = "ALS-4 duration (days)",
    y = "Maximum daily infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x        = element_text(size = 27, angle = 0, hjust = 0.5),
    axis.text.y        = element_text(size = 30),
    axis.title.x       = element_text(size = 30),
    axis.title.y       = element_text(size = 30),
    plot.title         = element_text(size = 35, hjust = 0.5, face = "plain"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.box.just    = "left",
    legend.background  = element_blank(),
    legend.text        = element_text(size = 27),
    plot.background    = element_blank(),
    plot.margin        = margin(20, 20, 20, 40),
    panel.spacing.x    = unit(1.2, "cm"),
    panel.grid.major.x = element_blank()
  )

print(p_peak_lines)

save_png(
  file.path(fig_dir, "als4_duration.png"),
  p_peak_lines,
  width = 6000, height = 3200, res = 300
)