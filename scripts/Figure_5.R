
suppressPackageStartupMessages({
  library(macpan2)
  library(tidyverse)
  library(shellpipes)
  library(ggthemes)
  library(scales)
  library(future)
  library(furrr)
  library(patchwork)
})

options(dplyr.summarise.inform = FALSE, scipen = 9999)
options(future.rng.onMisuse = "ignore")

loadEnvironments()

study_start <- as.Date("2021-12-15")
study_end   <- as.Date("2023-01-01")

strict_beta <- 0.245
coarse_step <- 7L

calibrator <- rdsRead("calibrate")
actual_beta <- as.numeric(rdsRead("base"))
no_npi_beta <- as.numeric(rdsRead("no_npis"))

vax_scenarios <- c(
  "actual" = "vax_actual_boosters",
  "earlier"  = "vax_early_boosters",
  "faster rate" = "vax_faster_rate",
  "slower rate" = "vax_slower_rate",
  "national rate"   = "vax_national_rate"
)

workers <- max(1L, min(length(vax_scenarios), future::availableCores()))
plan(multisession, workers = workers)

vax_colours <- c(
  "actual" = "#4D4D4D",
  "earlier"  = "#1F78B4",
  "faster rate" = "#FF7F00",
  "slower rate" = "#CC79A7",
  "national rate"   = "#33A02C"
)

vax_scen_legends <- c(
  "actual" = "Actual NL",
  "earlier"  = "Earlier",
  "faster rate" = "Faster rate",
  "slower rate" = "Slower rate",
  "national rate"   = "Canadian rate"
)

burn_in_start <- as.Date(burn_in_start)
study_start_index <- as.integer(study_start - burn_in_start) + 1L
study_end_index   <- as.integer(study_end - burn_in_start) + 1L
simulation_steps  <- study_end_index

load_booster <- function(stem) {
  file <- paste0(stem, ".rda")
  env <- new.env(parent = emptyenv())
  
  load(file, envir = env)
  
  if (!exists("new_booster", envir = env, inherits = FALSE)) {
    stop("Object 'new_booster' was not found in ", file, ".")
  }
  
  get("new_booster", envir = env, inherits = FALSE)
}

prepare_booster <- function(new_booster) {
  required_columns <- c("day", "daily_rate")
  missing_columns <- setdiff(required_columns, names(new_booster))
  
  booster <- new_booster |>
    transmute(
      day = as.integer(day),
      daily_rate = as.numeric(daily_rate)
    ) |>
    filter(
      is.finite(day),
      is.finite(daily_rate),
      day >= 0L,
      day < simulation_steps
    ) |>
    arrange(day) |>
    group_by(day) |>
    summarise(daily_rate = sum(daily_rate), .groups = "drop")
  
  if (!any(booster$day == 0L)) {
    booster <- bind_rows(
      tibble(day = 0L, daily_rate = 0),
      booster
    )
  }
  
  booster
}

booster_list <- lapply(vax_scenarios, load_booster) |>
  lapply(prepare_booster)

base_spec <- mp_optimized_spec(calibrator, "modified") |>
  mp_tmb_delete(phase = "before", at = 26L)

build_counterfactual_beta <- function(strict_days) {
  strict_days <- as.integer(strict_days)
  
  if (length(strict_days) != 1L || is.na(strict_days) || strict_days < 0L) {
    stop("strict_days must be one non-negative integer.")
  }
  
  beta <- actual_beta
  
  beta[study_start_index:study_end_index] <-
    no_npi_beta[study_start_index:study_end_index]
  
  if (strict_days > 0L) {
    strict_end_index <- min(
      study_start_index + strict_days - 1L,
      study_end_index
    )
    
    beta[study_start_index:strict_end_index] <- strict_beta
  }
  
  beta
}

simulate_hospitalizations <- function(booster, beta) {
  spec <- base_spec |>
    mp_tmb_insert(default = list(timeseries_beta = beta))
  
  spec$integers$booster_vac_changepoints <- booster$day
  spec$default$booster_vac_values <- booster$daily_rate
  
  trajectory <- spec |>
    mp_simulator(
      time_steps = simulation_steps,
      outputs = "cum_hosp_total"
    ) |>
    mp_trajectory() |>
    filter(matrix == "cum_hosp_total", is.finite(value)) |>
    transmute(
      date = burn_in_start + (time - 1L),
      value = as.numeric(value)
    ) |>
    arrange(date) |>
    mutate(daily_hospitalizations = value - lag(value))
  
  baseline <- trajectory |>
    filter(date < study_start) |>
    slice_max(date, n = 1L, with_ties = FALSE) |>
    pull(value)
  
  endpoint <- trajectory |>
    filter(date <= study_end) |>
    slice_max(date, n = 1L, with_ties = FALSE) |>
    pull(value)
  
  peak <- trajectory |>
    filter(date >= study_start, date <= study_end) |>
    summarise(value = max(daily_hospitalizations, na.rm = TRUE)) |>
    pull(value)
  
  if (length(baseline) == 0L) {
    baseline <- 0
  }
  
  total <- endpoint - baseline
  
  if (is.finite(total) && total < 0 && abs(total) < 1e-8) {
    total <- 0
  }
  
  list(
    total_hospitalizations = as.numeric(total),
    peak_hospitalizations = as.numeric(peak)
  )
}

actual_booster <- booster_list[["actual"]]

hospitalization_thresholds <- simulate_hospitalizations(
  booster = actual_booster,
  beta = actual_beta
)

hospitalization_threshold <- hospitalization_thresholds$total_hospitalizations
peak_hospitalization_threshold <- hospitalization_thresholds$peak_hospitalizations

# print(hospitalization_thresholds)
# print(peak_hospitalization_threshold)

max_strict_days <- 150
coarse_durations <- unique(c(
  seq.int(0L, max_strict_days, by = coarse_step),
  max_strict_days
))

analyse_scenario <- function(scenario) {
  simulation_cache <- new.env(parent = emptyenv())
  
  simulate_duration <- function(strict_days) {
    key <- as.character(strict_days)
    
    if (exists(key, envir = simulation_cache, inherits = FALSE)) {
      return(get(key, envir = simulation_cache, inherits = FALSE))
    }
    
    outcomes <- simulate_hospitalizations(
      booster = booster_list[[scenario]],
      beta = build_counterfactual_beta(strict_days)
    )
    
    assign(key, outcomes, envir = simulation_cache)
    outcomes
  }
  
  message("Simulating ", scenario, "...")
  
  coarse_result <- map_dfr(coarse_durations, function(strict_days) {
    outcomes <- simulate_duration(strict_days)
    
    tibble(
      vax_scenario = scenario,
      strict_npi_days = strict_days,
      total_hospitalizations = outcomes$total_hospitalizations,
      peak_hospitalizations = outcomes$peak_hospitalizations
    )
  })
  
  first_total_hit <- coarse_result |>
    filter(total_hospitalizations <= hospitalization_threshold) |>
    arrange(strict_npi_days) |>
    slice_head(n = 1L)
  
  first_peak_hit <- coarse_result |>
    filter(peak_hospitalizations <= peak_hospitalization_threshold) |>
    arrange(strict_npi_days) |>
    slice_head(n = 1L)
  
  total_candidate_days <- integer()
  peak_candidate_days <- integer()
  
  if (nrow(first_total_hit) > 0L) {
    total_upper_day <- first_total_hit$strict_npi_days[[1]]
    total_lower_day <- max(0L, total_upper_day - coarse_step + 1L)
    total_candidate_days <- seq.int(total_lower_day, total_upper_day)
  }
  
  if (nrow(first_peak_hit) > 0L) {
    peak_upper_day <- first_peak_hit$strict_npi_days[[1]]
    peak_lower_day <- max(0L, peak_upper_day - coarse_step + 1L)
    peak_candidate_days <- seq.int(peak_lower_day, peak_upper_day)
  }
  
  candidate_days <- sort(unique(c(
    total_candidate_days,
    peak_candidate_days
  )))
  
  candidate_results <- map_dfr(candidate_days, function(strict_days) {
    outcomes <- simulate_duration(strict_days)
    
    tibble(
      strict_npi_days = strict_days,
      total_hospitalizations = outcomes$total_hospitalizations,
      peak_hospitalizations = outcomes$peak_hospitalizations
    )
  })
  
  minimum_total <- candidate_results |>
    filter(total_hospitalizations <= hospitalization_threshold) |>
    arrange(strict_npi_days) |>
    slice_head(n = 1L)
  
  if (nrow(minimum_total) == 0L) {
    minimum_total <- tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      total_hospitalizations = NA_real_,
      threshold_met = FALSE
    )
  } else {
    minimum_total <- minimum_total |>
      transmute(
        vax_scenario = scenario,
        minimum_strict_npi_days = strict_npi_days,
        total_hospitalizations,
        threshold_met = TRUE
      )
  }
  
  minimum_peak <- candidate_results |>
    filter(peak_hospitalizations <= peak_hospitalization_threshold) |>
    arrange(strict_npi_days) |>
    slice_head(n = 1L)
  
  if (nrow(minimum_peak) == 0L) {
    minimum_peak <- tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      peak_hospitalizations = NA_real_,
      threshold_met = FALSE
    )
  } else {
    minimum_peak <- minimum_peak |>
      transmute(
        vax_scenario = scenario,
        minimum_strict_npi_days = strict_npi_days,
        peak_hospitalizations,
        threshold_met = TRUE
      )
  }
  
  list(
    coarse_result = coarse_result,
    minimum_total = minimum_total,
    minimum_peak = minimum_peak
  )
}

scenario_results <- future_map(
  names(vax_scenarios),
  analyse_scenario,
  .options = furrr_options(seed = TRUE)
)

coarse_results <- map_dfr(
  scenario_results,
  "coarse_result"
) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)),
    threshold_met = total_hospitalizations <= hospitalization_threshold,
    peak_threshold_met = peak_hospitalizations <= peak_hospitalization_threshold
  )

minimum_duration <- map_dfr(
  scenario_results,
  "minimum_total"
) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)),
    hospitalization_threshold = hospitalization_threshold,
    margin_below_threshold = hospitalization_threshold - total_hospitalizations
  )

minimum_peak_duration <- map_dfr(
  scenario_results,
  "minimum_peak"
) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)),
    peak_hospitalization_threshold = peak_hospitalization_threshold,
    margin_below_threshold = peak_hospitalization_threshold - peak_hospitalizations
  )

plan(sequential)

threshold_data <- (tibble(
  scenario = "Actual Boosters + Actual NPIs"
  , study_start = study_start
  , study_end = study_end
  , cumulative_hospitalizations = hospitalization_threshold
  , peak_hospitalizations = peak_hospitalization_threshold
)
)

cumulative_label_line <- max(
  coarse_results$total_hospitalizations,
  na.rm = TRUE
) * 1.08

cumulative_label_y <- max(
  coarse_results$total_hospitalizations,
  na.rm = TRUE
) * 1.11

peak_label_line <- max(
  coarse_results$peak_hospitalizations,
  na.rm = TRUE
) * 1.08

peak_label_y <- max(
  coarse_results$peak_hospitalizations,
  na.rm = TRUE
) * 1.11


crossing_points <- (minimum_duration 
  |> filter(threshold_met) 
  |> mutate(
    label = paste0(minimum_strict_npi_days),
    label_x = case_when(
      vax_scenario == "earlier" ~ minimum_strict_npi_days - 1,
      vax_scenario == "faster rate" ~ minimum_strict_npi_days + 1,
      vax_scenario == "actual" ~ minimum_strict_npi_days + 2,
      vax_scenario == "slower rate" ~ minimum_strict_npi_days - 4,
      vax_scenario == "national rate" ~ minimum_strict_npi_days + 1,
      TRUE ~ as.numeric(minimum_strict_npi_days)
    ),
    label_y = case_when(
      vax_scenario %in% c("earlier", "slower rate") ~ cumulative_label_y * 1.04,
      TRUE ~ cumulative_label_y
    ),
    segment_yend = cumulative_label_line,
    label_hjust = 0.5
  )
)
peak_crossing_points <- (minimum_peak_duration 
  |> filter(threshold_met) 
  |> mutate(
   label = paste0(minimum_strict_npi_days),
   label_x = case_when(
     vax_scenario == "faster rate" ~ minimum_strict_npi_days - 2,
     vax_scenario == "earlier" ~ minimum_strict_npi_days + 1,
     vax_scenario == "actual" ~ minimum_strict_npi_days + 3,
     TRUE ~ as.numeric(minimum_strict_npi_days)
   ),
   label_y = peak_label_y,
   segment_yend = peak_label_line,
   label_hjust = 0.5
  )
)

figure4A <- (
  ggplot(coarse_results, aes(x = strict_npi_days,y = total_hospitalizations,colour = vax_scenario)) 
  + geom_hline(yintercept = hospitalization_threshold, linetype = "dashed",linewidth = 2.5,colour = "grey40") 
  + geom_line(linewidth = 1.5) 
  + geom_point(size = 4) 
  + geom_segment(
    data = crossing_points,
    aes(
      x = minimum_strict_npi_days,
      xend = minimum_strict_npi_days,
      y = hospitalization_threshold,
      yend = segment_yend,
      colour = vax_scenario
    ),
    linetype = "dotdash",
    linewidth = 2,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  + geom_point(data = crossing_points,aes(
    x = minimum_strict_npi_days,
    y = total_hospitalizations,
    colour = vax_scenario)
    , stroke = 1.8
    , shape = 8
    , size = 12
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + geom_text(
    data = crossing_points,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      colour = vax_scenario,
      hjust = label_hjust
    )
    , angle = 0
    , vjust = 0
    , size = 10
    , fontface = "bold"
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + scale_colour_manual(
    values = vax_colours,
    breaks = names(vax_scenarios),
    labels = vax_scen_legends,
    name = NULL,
    drop = FALSE
  )
  + scale_x_continuous(
    breaks = seq(0, max_strict_days, by = 14),
    limits = c(0, max_strict_days),
    expand = expansion(mult = c(0.005, 0.02))
  ) 
  + scale_y_continuous(
    labels = comma,
    limits = c(1500, NA),
    expand = expansion(mult = c(0.02, 0.25))
  ) 
  + labs(
    title = "Strict NPI duration required under different booster scenarios",
    x = "Strict NPI duration (days)",
    y = "Cumulative hospitalizations",
    tag = "A") 
  + guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 1, size = 4))) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 25, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.spacing.y = unit(0.8, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.tag = element_text(size = 30, face = "plain"),
    plot.margin = margin(12, 35, 12, 12)
  )
)


peak_y_min <- min(coarse_results$peak_hospitalizations, na.rm = TRUE)

figure4B <- (
  ggplot(coarse_results, aes(x = strict_npi_days,y = peak_hospitalizations,colour = vax_scenario)) 
  + geom_hline(yintercept = peak_hospitalization_threshold,linetype = "dashed",linewidth = 2.5,colour = "grey40") 
  + geom_line(linewidth = 1.5) 
  + geom_point(size = 4) 
  + geom_segment(
    data = peak_crossing_points,
    aes(
      x = minimum_strict_npi_days,
      xend = minimum_strict_npi_days,
      y = peak_hospitalization_threshold,
      yend = segment_yend,
      colour = vax_scenario
    ),
    linetype = "dotdash",
    linewidth = 2,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  + geom_point(data = peak_crossing_points,aes(
    x = minimum_strict_npi_days,
    y = peak_hospitalizations,
    colour = vax_scenario)
    , stroke = 1.8
    , shape = 8
    , size = 12
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + geom_text(
    data = peak_crossing_points,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      colour = vax_scenario,
      hjust = label_hjust
    )
    , angle = 0
    , vjust = 0
    , size = 10
    , fontface = "bold"
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + scale_colour_manual(
    values = vax_colours,
    breaks = names(vax_scenarios),
    labels = vax_scen_legends,
    name = NULL,
    drop = FALSE
  )
  + scale_x_continuous(
    breaks = seq(0, max_strict_days, by = 14),
    limits = c(0, max_strict_days),
    expand = expansion(mult = c(0.005, 0.02))
  ) 
  + scale_y_continuous(
    labels = comma,
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.25))
  ) 
  + labs(
    x = "Strict NPI duration (days)",
    y = "Peak hospital admissions",
    tag = "B") 
  + guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 1, size = 4))) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 25, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.spacing.y = unit(0.8, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.tag = element_text(size = 30, face = "plain"),
    plot.margin = margin(12, 35, 12, 12)
  )
)


figure4C <- (
  ggplot(coarse_results, aes(
    x = total_hospitalizations,
    y = peak_hospitalizations,
    colour = vax_scenario,
    group = vax_scenario
  )) 
  + geom_vline(
    xintercept = hospitalization_threshold,
    linetype = "dashed",
    linewidth = 2.5,
    colour = "grey40"
  )
  + geom_hline(
    yintercept = peak_hospitalization_threshold,
    linetype = "dashed",
    linewidth = 2.5,
    colour = "grey40"
  )
  + geom_line(linewidth = 1.5) 
  + geom_point(size = 4) 
  + scale_colour_manual(
    values = vax_colours,
    breaks = names(vax_scenarios),
    labels = vax_scen_legends,
    name = NULL,
    drop = FALSE
  )
  + scale_x_continuous(
    labels = comma,
    limits = c(NA, 4000), 
    expand = expansion(mult = c(0.02, 0.04))
  ) 
  + scale_y_continuous(
    labels = comma,
    limits = c(0, 85),
    expand = expansion(mult = c(0.02, 0.12))
  ) 
  + labs(
    x = "Cumulative hospitalizations",
    y = "Peak hospital admissions",
    tag = "C"
  ) 
  + guides(colour = guide_legend(
    override.aes = list(shape = 16, linetype = 1, size = 4)
  )) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 25, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.spacing.y = unit(0.8, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.tag = element_text(size = 30, face = "plain"),
    plot.margin = margin(12, 35, 12, 12)
  )
)

figure4 <- (
  figure4A / figure4B / figure4C
  + plot_layout(guides = "collect", heights = c(1, 1, 1))
  & theme(
    legend.position = "bottom",
    legend.text = element_text(size = 25)
  )
)

png("../figures/figure5.png",width = 4500,height = 5000,res = 300,bg = "white",type = "cairo")
print(figure4)
dev.off()