# 
# suppressPackageStartupMessages({
#   library(macpan2)
#   library(tidyverse)
#   library(shellpipes)
#   library(ggthemes)
#   library(scales)
#   library(grid)
# })
# 
# options(dplyr.summarise.inform = FALSE)
# options(scipen = 9999)
# 
# loadEnvironments()
# 
# calib     <- rdsRead("calibrate")
# base_beta <- rdsRead("base")
# 
# strict_npi_constant_value <- 0.245
# 
# study_start  <- as.Date("2021-12-15")
# primary_end  <- as.Date("2022-05-22")
# extended_end <- as.Date("2023-01-01")
# 
# primary_threshold  <- 1612
# extended_threshold <- 3055
# 
# primary_threshold_label  <- "Threshold"
# extended_threshold_label <- "Threshold"
# 
# max_duration <- length(base_beta) - as.integer(step_start) + 1L
# 
# durations <- seq(0L, 180L, by = 7L)
# durations <- durations[durations <= max_duration]
# 
# vax_scenarios <- c(
#   "Actual Boosters" = "vax_actual_boosters",
#   "Early Boosters"  = "vax_early_boosters",
#   "Faster Boosters" = "vax_faster_rate",
#   "Slower Boosters" = "vax_slower_rate",
#   "National Rate"   = "vax_national_rate"
# )
# 
# vax_cols <- c(
#   "Actual Boosters" = "#4D4D4D",
#   "Early Boosters"  = "#1f78b4",
#   "Faster Boosters" = "#ff7f00",
#   "Slower Boosters" = "magenta",
#   "National Rate"   = "#33a02c"
# )
# 
# load_booster <- function(stem) {
#   env <- new.env(parent = emptyenv())
#   load(paste0(stem, ".rda"), envir = env)
#   
#   if (!exists("new_booster", envir = env, inherits = FALSE)) {
#     stop("Object 'new_booster' was not found in ", paste0(stem, ".rda"))
#   }
#   
#   get("new_booster", envir = env, inherits = FALSE)
# }
# 
# booster_list <- lapply(vax_scenarios, load_booster)
# 
# base_spec <- mp_optimized_spec(calib, "modified") |>
#   mp_tmb_delete(phase = "before", at = 26L)
# 
# make_beta <- function(duration) {
#   duration <- as.integer(duration)
#   
#   beta <- base_beta
#   
#   if (duration > 0L) {
#     start_index <- as.integer(step_start)
#     end_index <- min(start_index + duration - 1L, length(beta))
#     beta[start_index:end_index] <- strict_npi_constant_value
#   }
#   
#   beta
# }
# 
# prepare_booster_inputs <- function(new_booster) {
#   required_columns <- c("day", "daily_rate")
#   missing_columns <- setdiff(required_columns, names(new_booster))
#   
#   booster_inputs <- new_booster |>
#     transmute(
#       day = as.integer(day),
#       daily_rate = as.numeric(daily_rate)
#     ) |>
#     filter(is.finite(day), is.finite(daily_rate)) |>
#     arrange(day) |>
#     group_by(day) |>
#     summarise(daily_rate = sum(daily_rate), .groups = "drop")
#   
#   if (!any(booster_inputs$day == 0L)) {
#     booster_inputs <- bind_rows(
#       tibble(day = 0L, daily_rate = 0),
#       booster_inputs
#     )
#   }
#   
#   booster_inputs
# }
# 
# total_hosp_one <- function(new_booster, duration, study_end) {
#   study_end <- as.Date(study_end)
#   
#   new_beta <- make_beta(duration)
#   booster_inputs <- prepare_booster_inputs(new_booster)
#   
#   spec <- base_spec |>
#     mp_tmb_insert(default = list(timeseries_beta = new_beta))
#   
#   spec$integers$booster_vac_changepoints <- booster_inputs$day
#   spec$default$booster_vac_values <- booster_inputs$daily_rate
#   
#   trajectory <- spec |>
#     mp_simulator(
#       time_steps = sim_horizon,
#       outputs = "cum_hosp_total"
#     ) |>
#     mp_trajectory() |>
#     mutate(date = as.Date(burn_in_start) + (time - 1)) |>
#     filter(matrix == "cum_hosp_total", is.finite(value)) |>
#     arrange(date)
#   
#   baseline_value <- trajectory |>
#     filter(date < study_start) |>
#     slice_max(order_by = date, n = 1, with_ties = FALSE) |>
#     pull(value)
#   
#   endpoint_value <- trajectory |>
#     filter(date <= study_end) |>
#     slice_max(order_by = date, n = 1, with_ties = FALSE) |>
#     pull(value)
#   
#   if (length(baseline_value) == 0L) {
#     baseline_value <- 0
#   }
#   
#   if (length(endpoint_value) == 0L) {
#     return(NA_real_)
#   }
#   
#   total_hospitalizations <- endpoint_value - baseline_value
#   
#   if (
#     is.finite(total_hospitalizations) &&
#     total_hospitalizations < 0 &&
#     abs(total_hospitalizations) < 1e-8
#   ) {
#     total_hospitalizations <- 0
#   }
#   
#   total_hospitalizations
# }
# 
# run_duration_sweep <- function(study_end, window_name) {
#   map_dfr(names(vax_scenarios), function(scenario_name) {
#     
#     booster <- booster_list[[scenario_name]]
#     
#     map_dfr(durations, function(days) {
#       tibble(
#         vax_scenario = scenario_name,
#         duration = days,
#         total_hosp = total_hosp_one(
#           new_booster = booster,
#           duration = days,
#           study_end = study_end
#         )
#       )
#     })
#   }) |>
#     mutate(
#       study_end = as.Date(study_end),
#       window = window_name,
#       vax_scenario = factor(vax_scenario, levels = names(vax_scenarios))
#     )
# }
# 
# duration_to_meet <- function(duration, total_hosp, threshold) {
#   keep <- is.finite(duration) & is.finite(total_hosp)
#   duration <- duration[keep]
#   total_hosp <- total_hosp[keep]
#   
#   order_index <- order(duration)
#   duration <- duration[order_index]
#   total_hosp <- total_hosp[order_index]
#   
#   if (all(total_hosp <= threshold)) {
#     return(min(duration))
#   }
#   
#   if (all(total_hosp > threshold)) {
#     return(NA_real_)
#   }
#   
#   crossing_index <- which(total_hosp <= threshold)[1]
#   
#   if (crossing_index == 1L) {
#     return(duration[1])
#   }
#   
#   x1 <- duration[crossing_index - 1L]
#   x2 <- duration[crossing_index]
#   y1 <- total_hosp[crossing_index - 1L]
#   y2 <- total_hosp[crossing_index]
#   
#   if (isTRUE(all.equal(y1, y2))) {
#     return(x2)
#   }
#   
#   x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
# }
# 
# summarise_threshold <- function(sweep_df, threshold) {
#   sweep_df |>
#     group_by(vax_scenario) |>
#     summarise(
#       hosp_days_needed = duration_to_meet(
#         duration = duration,
#         total_hosp = total_hosp,
#         threshold = threshold
#       ),
#       .groups = "drop"
#     ) |>
#     mutate(hosp_days_needed = round(hosp_days_needed))
# }
# 
# make_duration_plot <- function(
#     sweep_df,
#     threshold,
#     threshold_label,
#     study_end,
#     output_file,
#     plot_title,
#     legend_position = c(0.88, 0.78),
#     legend_direction = "vertical"
# ) {
#   required_tbl <- summarise_threshold(
#     sweep_df = sweep_df,
#     threshold = threshold
#   )
#   
#   cross_df <- required_tbl |>
#     filter(is.finite(hosp_days_needed)) |>
#     transmute(
#       vax_scenario,
#       strict_npi_days = hosp_days_needed,
#       cross_lab = paste0("~", hosp_days_needed, " days"),
#       label_y = threshold * 0.18
#     )
#   
#   plot_object <- ggplot(
#     sweep_df,
#     aes(x = duration, y = total_hosp, colour = vax_scenario)
#   ) +
#     geom_line(linewidth = 2) +
#     geom_point(size = 4) +
#     geom_hline(
#       yintercept = threshold,
#       colour = "red",
#       linetype = "dashed",
#       linewidth = 2
#     ) +
#     annotate(
#       "text",
#       x = max(sweep_df$duration, na.rm = TRUE),
#       y = threshold,
#       label = paste0(threshold_label, ": ", comma(round(threshold))),
#       hjust = 1,
#       vjust = -0.6,
#       size = 8,
#       colour = "red"
#     ) +
#     geom_segment(
#       data = cross_df,
#       aes(
#         x = strict_npi_days,
#         xend = strict_npi_days,
#         y = 0,
#         yend = threshold,
#         colour = vax_scenario
#       ),
#       linetype = "dashed",
#       linewidth = 0.8,
#       inherit.aes = FALSE,
#       show.legend = FALSE
#     ) +
#     geom_point(
#       data = cross_df,
#       aes(x = strict_npi_days, y = threshold, colour = vax_scenario),
#       shape = 8,
#       size = 8,
#       stroke = 1.8,
#       inherit.aes = FALSE,
#       show.legend = FALSE
#     ) +
#     geom_text(
#       data = cross_df,
#       aes(
#         x = strict_npi_days,
#         y = label_y,
#         label = cross_lab,
#         colour = vax_scenario
#       ),
#       angle = 90,
#       vjust = -0.5,
#       size = 8,
#       inherit.aes = FALSE,
#       show.legend = FALSE
#     ) +
#     scale_colour_manual(
#       values = vax_cols,
#       breaks = names(vax_scenarios),
#       name = NULL,
#       drop = FALSE
#     ) +
#     scale_x_continuous(
#       breaks = seq(0, max(durations), by = 14),
#       expand = expansion(mult = c(0.01, 0.03))
#     ) +
#     scale_y_continuous(
#       labels = comma,
#       expand = expansion(mult = c(0, 0.12))
#     ) +
#     labs(
#       title = plot_title,
#       x = "Strict NPI duration (days)",
#       y = "Total hospitalizations"
#     ) +
#     guides(
#       colour = guide_legend(
#         override.aes = list(shape = 16, linetype = 1, size = 4)
#       )
#     ) +
#     theme_clean() +
#     theme(
#       axis.text.x = element_text(size = 30, hjust = 0.5, angle = 90),
#       axis.text.y = element_text(size = 30),
#       axis.title.x = element_text(size = 30),
#       axis.title.y = element_text(size = 30),
#       plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
#       plot.subtitle = element_text(size = 25, hjust = 0.5),
#       legend.position = legend_position,
#       legend.justification = "center",
#       legend.direction = legend_direction,
#       legend.text = element_text(size = 27),
#       legend.spacing.y = unit(0.8, "cm"),
#       legend.key.spacing.y = unit(0.5, "cm"),
#       legend.background = element_blank(),
#       panel.grid.major.x = element_blank(),
#       plot.background = element_blank(),
#       plot.margin = margin(12, 35, 12, 12)
#     )
#   
#   png(
#     output_file,
#     width = 6000,
#     height = 3600,
#     res = 300,
#     bg = "white",
#     type = "cairo"
#   )
#   
#   print(plot_object)
#   dev.off()
#   
#   list(
#     plot = plot_object,
#     required_duration = required_tbl
#   )
# }
# 
# primary_sweep <- run_duration_sweep(
#   study_end = primary_end,
#   window_name = "Primary: 15 Dec 2021 to 22 May 2022"
# )
# 
# extended_sweep <- run_duration_sweep(
#   study_end = extended_end,
#   window_name = "Sensitivity: 15 Dec 2021 to 1 Jan 2023"
# )
# 
# primary_result <- make_duration_plot(
#   sweep_df = primary_sweep,
#   threshold = primary_threshold,
#   threshold_label = primary_threshold_label,
#   study_end = primary_end,
#   output_file = "../figures/figure4a.png",
#   plot_title = "Total hospitalizations across strict NPI duration and booster scenarios",
#   legend_position = c(0.88, 0.78),
#   legend_direction = "vertical"
# )
# 
# extended_result <- make_duration_plot(
#   sweep_df = extended_sweep,
#   threshold = extended_threshold,
#   threshold_label = extended_threshold_label,
#   study_end = extended_end,
#   output_file = "../figures/figure4b.png",
#   plot_title = "Total hospitalizations across strict NPI duration and booster scenarios",
#   legend_position = "top",
#   legend_direction = "horizontal"
# )
# 





suppressPackageStartupMessages({
  library(macpan2)
  library(tidyverse)
  library(shellpipes)
  library(ggthemes)
  library(scales)
  library(grid)
})

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

calib     <- rdsRead("calibrate")
base_beta <- rdsRead("base")

strict_npi_constant_value <- 0.245

study_start  <- as.Date("2021-12-15")
primary_end  <- as.Date("2022-05-22")
extended_end <- as.Date("2023-01-01")

primary_threshold  <- 1612
extended_threshold <- 3055

primary_threshold_label  <- "Threshold"
extended_threshold_label <- "Threshold"

max_duration <- length(base_beta) - as.integer(step_start) + 1L

durations <- seq(0L, 180L, by = 7L)
durations <- durations[durations <= max_duration]

vax_scenarios <- c(
  "Actual Boosters" = "vax_actual_boosters",
  "Early Boosters"  = "vax_early_boosters",
  "Faster Boosters" = "vax_faster_rate",
  "Slower Boosters" = "vax_slower_rate",
  "National Rate"   = "vax_national_rate"
)

vax_cols <- c(
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "#1f78b4",
  "Faster Boosters" = "#ff7f00",
  "Slower Boosters" = "magenta",
  "National Rate"   = "#33a02c"
)

load_booster <- function(stem) {
  env <- new.env(parent = emptyenv())
  load(paste0(stem, ".rda"), envir = env)
  
  if (!exists("new_booster", envir = env, inherits = FALSE)) {
    stop("Object 'new_booster' was not found in ", paste0(stem, ".rda"))
  }
  
  get("new_booster", envir = env, inherits = FALSE)
}

booster_list <- lapply(vax_scenarios, load_booster)

base_spec <- mp_optimized_spec(calib, "modified") |>
  mp_tmb_delete(phase = "before", at = 26L)

make_beta <- function(duration) {
  duration <- as.integer(duration)
  
  beta <- base_beta
  
  if (duration > 0L) {
    start_index <- as.integer(step_start)
    end_index <- min(start_index + duration - 1L, length(beta))
    beta[start_index:end_index] <- strict_npi_constant_value
  }
  
  beta
}

prepare_booster_inputs <- function(new_booster) {
  required_columns <- c("day", "daily_rate")
  missing_columns <- setdiff(required_columns, names(new_booster))
  
  booster_inputs <- new_booster |>
    transmute(
      day = as.integer(day),
      daily_rate = as.numeric(daily_rate)
    ) |>
    filter(is.finite(day), is.finite(daily_rate)) |>
    arrange(day) |>
    group_by(day) |>
    summarise(daily_rate = sum(daily_rate), .groups = "drop")
  
  if (!any(booster_inputs$day == 0L)) {
    booster_inputs <- bind_rows(
      tibble(day = 0L, daily_rate = 0),
      booster_inputs
    )
  }
  
  booster_inputs
}

booster_inputs_list <- lapply(booster_list, prepare_booster_inputs)

simulate_trajectory <- function(booster_inputs, duration) {
  new_beta <- make_beta(duration)
  
  spec <- base_spec |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
  
  spec$integers$booster_vac_changepoints <- booster_inputs$day
  spec$default$booster_vac_values <- booster_inputs$daily_rate
  
  spec |>
    mp_simulator(
      time_steps = sim_horizon,
      outputs = "cum_hosp_total"
    ) |>
    mp_trajectory() |>
    mutate(date = as.Date(burn_in_start) + (time - 1)) |>
    filter(matrix == "cum_hosp_total", is.finite(value)) |>
    arrange(date)
}

total_hosp_from_trajectory <- function(trajectory, study_end) {
  study_end <- as.Date(study_end)
  
  baseline_value <- trajectory |>
    filter(date < study_start) |>
    slice_max(order_by = date, n = 1, with_ties = FALSE) |>
    pull(value)
  
  endpoint_value <- trajectory |>
    filter(date <= study_end) |>
    slice_max(order_by = date, n = 1, with_ties = FALSE) |>
    pull(value)
  
  if (length(baseline_value) == 0L) {
    baseline_value <- 0
  }
  
  if (length(endpoint_value) == 0L) {
    return(NA_real_)
  }
  
  total_hospitalizations <- endpoint_value - baseline_value
  
  if (
    is.finite(total_hospitalizations) &&
    total_hospitalizations < 0 &&
    abs(total_hospitalizations) < 1e-8
  ) {
    total_hospitalizations <- 0
  }
  
  total_hospitalizations
}

run_all_sweeps <- function() {
  map_dfr(names(vax_scenarios), function(scenario_name) {
    booster_inputs <- booster_inputs_list[[scenario_name]]
    
    map_dfr(durations, function(days) {
      trajectory <- simulate_trajectory(booster_inputs, days)
      
      tibble(
        vax_scenario = scenario_name,
        duration = days,
        primary_hosp  = total_hosp_from_trajectory(trajectory, primary_end),
        extended_hosp = total_hosp_from_trajectory(trajectory, extended_end)
      )
    })
  }) |>
    mutate(vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)))
}

all_sweeps <- run_all_sweeps()

primary_sweep <- all_sweeps |>
  transmute(vax_scenario, duration, total_hosp = primary_hosp,
            study_end = as.Date(primary_end),
            window = "Primary: 15 Dec 2021 to 22 May 2022")

extended_sweep <- all_sweeps |>
  transmute(vax_scenario, duration, total_hosp = extended_hosp,
            study_end = as.Date(extended_end),
            window = "Sensitivity: 15 Dec 2021 to 1 Jan 2023")

duration_to_meet <- function(duration, total_hosp, threshold) {
  keep <- is.finite(duration) & is.finite(total_hosp)
  duration <- duration[keep]
  total_hosp <- total_hosp[keep]
  
  order_index <- order(duration)
  duration <- duration[order_index]
  total_hosp <- total_hosp[order_index]
  
  if (all(total_hosp <= threshold)) {
    return(min(duration))
  }
  
  if (all(total_hosp > threshold)) {
    return(NA_real_)
  }
  
  crossing_index <- which(total_hosp <= threshold)[1]
  
  if (crossing_index == 1L) {
    return(duration[1])
  }
  
  x1 <- duration[crossing_index - 1L]
  x2 <- duration[crossing_index]
  y1 <- total_hosp[crossing_index - 1L]
  y2 <- total_hosp[crossing_index]
  
  if (isTRUE(all.equal(y1, y2))) {
    return(x2)
  }
  
  x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
}

summarise_threshold <- function(sweep_df, threshold) {
  sweep_df |>
    group_by(vax_scenario) |>
    summarise(
      hosp_days_needed = duration_to_meet(
        duration = duration,
        total_hosp = total_hosp,
        threshold = threshold
      ),
      .groups = "drop"
    ) |>
    mutate(hosp_days_needed = round(hosp_days_needed))
}

make_duration_plot <- function(
    sweep_df,
    threshold,
    threshold_label,
    study_end,
    output_file,
    plot_title,
    legend_position = c(0.88, 0.78),
    legend_direction = "vertical"
) {
  required_tbl <- summarise_threshold(
    sweep_df = sweep_df,
    threshold = threshold
  )
  
  cross_df <- required_tbl |>
    filter(is.finite(hosp_days_needed)) |>
    transmute(
      vax_scenario,
      strict_npi_days = hosp_days_needed,
      cross_lab = paste0("~", hosp_days_needed, " days"),
      label_y = threshold * 0.18
    )
  
  plot_object <- ggplot(
    sweep_df,
    aes(x = duration, y = total_hosp, colour = vax_scenario)
  ) +
    geom_line(linewidth = 2) +
    geom_point(size = 4) +
    geom_hline(
      yintercept = threshold,
      colour = "red",
      linetype = "dashed",
      linewidth = 2
    ) +
    annotate(
      "text",
      x = max(sweep_df$duration, na.rm = TRUE),
      y = threshold,
      label = paste0(threshold_label, ": ", comma(round(threshold))),
      hjust = 1,
      vjust = -0.6,
      size = 8,
      colour = "red"
    ) +
    geom_segment(
      data = cross_df,
      aes(
        x = strict_npi_days,
        xend = strict_npi_days,
        y = 0,
        yend = threshold,
        colour = vax_scenario
      ),
      linetype = "dashed",
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = cross_df,
      aes(x = strict_npi_days, y = threshold, colour = vax_scenario),
      shape = 8,
      size = 8,
      stroke = 1.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = cross_df,
      aes(
        x = strict_npi_days,
        y = label_y,
        label = cross_lab,
        colour = vax_scenario
      ),
      angle = 90,
      vjust = -0.5,
      size = 8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = vax_cols,
      breaks = names(vax_scenarios),
      name = NULL,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq(0, max(durations), by = 14),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = plot_title,
      x = "Strict NPI duration (days)",
      y = "Total hospitalizations"
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(shape = 16, linetype = 1, size = 4)
      )
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 0.5, angle = 90),
      axis.text.y = element_text(size = 30),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 25, hjust = 0.5),
      legend.position = legend_position,
      legend.justification = "center",
      legend.direction = legend_direction,
      legend.text = element_text(size = 27),
      legend.spacing.y = unit(0.8, "cm"),
      legend.key.spacing.y = unit(0.5, "cm"),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(12, 35, 12, 12)
    )
  
  png(
    output_file,
    width = 6000,
    height = 3600,
    res = 300,
    bg = "white",
    type = "cairo"
  )
  
  print(plot_object)
  dev.off()
  
  list(
    plot = plot_object,
    required_duration = required_tbl
  )
}


primary_result <- make_duration_plot(
  sweep_df = primary_sweep,
  threshold = primary_threshold,
  threshold_label = primary_threshold_label,
  study_end = primary_end,
  output_file = "../figures/figure4a.png",
  plot_title = "Total hospitalizations across strict NPI duration and booster scenarios",
  legend_position = c(0.88, 0.78),
  legend_direction = "vertical"
)

extended_result <- make_duration_plot(
  sweep_df = extended_sweep,
  threshold = extended_threshold,
  threshold_label = extended_threshold_label,
  study_end = extended_end,
  output_file = "../figures/figure4b.png",
  plot_title = "Total hospitalizations across strict NPI duration and booster scenarios",
  legend_position = "top",
  legend_direction = "horizontal"
)