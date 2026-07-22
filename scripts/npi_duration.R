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
  "Actual Boosters" = "vax_actual_boosters",
  "Early Boosters"  = "vax_early_boosters",
  "Faster Boosters" = "vax_faster_rate",
  "Slower Boosters" = "vax_slower_rate",
  "National Rate"   = "vax_national_rate"
)

workers <- max(1L, min(length(vax_scenarios), future::availableCores()))
plan(multisession, workers = workers)

vax_colours <- c(
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "#1F78B4",
  "Faster Boosters" = "#FF7F00",
  "Slower Boosters" = "#CC79A7",
  "National Rate"   = "#33A02C"
)

if (!exists("burn_in_start")) {
  stop("The object 'burn_in_start' was not loaded by loadEnvironments().")
}

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
  
  if (length(missing_columns) > 0L) {
    stop(
      "The booster data are missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }
  
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
  
  if (length(endpoint) != 1L) {
    stop("The simulation did not produce a hospitalization value at the study endpoint.")
  }
  
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

actual_booster <- booster_list[["Actual Boosters"]]

hospitalization_thresholds <- simulate_hospitalizations(
  booster = actual_booster,
  beta = actual_beta
)

hospitalization_threshold <- hospitalization_thresholds$total_hospitalizations
peak_hospitalization_threshold <- hospitalization_thresholds$peak_hospitalizations

max_strict_days <- 150
coarse_durations <- unique(c(
  seq.int(0L, max_strict_days, by = coarse_step),
  max_strict_days
))

simulation_cache <- new.env(parent = emptyenv())

cache_key <- function(scenario, strict_days) {
  paste(scenario, strict_days, sep = "::")
}

simulate_duration <- function(scenario, strict_days) {
  key <- cache_key(scenario, strict_days)
  
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

run_coarse_sweep <- function() {
  future_map_dfr(names(vax_scenarios), function(scenario) {
    message("Simulating ", scenario, "...")
    
    map_dfr(coarse_durations, function(strict_days) {
      outcomes <- simulate_duration(scenario, strict_days)
      
      tibble(
        vax_scenario = scenario,
        strict_npi_days = strict_days,
        total_hospitalizations = outcomes$total_hospitalizations,
        peak_hospitalizations = outcomes$peak_hospitalizations
      )
    })
  }) |>
    mutate(
      vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)),
      threshold_met = total_hospitalizations <= hospitalization_threshold,
      peak_threshold_met = peak_hospitalizations <= peak_hospitalization_threshold
    )
}

find_minimum_duration <- function(scenario, coarse_results) {
  scenario_results <- coarse_results |>
    filter(vax_scenario == scenario) |>
    arrange(strict_npi_days)
  
  first_coarse_hit <- scenario_results |>
    filter(total_hospitalizations <= hospitalization_threshold) |>
    slice_head(n = 1L)
  
  if (nrow(first_coarse_hit) == 0L) {
    return(tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      total_hospitalizations = NA_real_,
      threshold_met = FALSE
    ))
  }
  
  upper_day <- first_coarse_hit$strict_npi_days[[1]]
  lower_day <- max(0L, upper_day - coarse_step + 1L)
  candidate_days <- seq.int(lower_day, upper_day)
  
  candidate_results <- map_dfr(candidate_days, function(strict_days) {
    outcomes <- simulate_duration(scenario, strict_days)
    
    tibble(
      strict_npi_days = strict_days,
      total_hospitalizations = outcomes$total_hospitalizations
    )
  }) |>
    filter(total_hospitalizations <= hospitalization_threshold) |>
    arrange(strict_npi_days)
  
  if (nrow(candidate_results) == 0L) {
    return(tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      total_hospitalizations = NA_real_,
      threshold_met = FALSE
    ))
  }
  
  candidate_results |>
    slice_head(n = 1L) |>
    transmute(
      vax_scenario = scenario,
      minimum_strict_npi_days = strict_npi_days,
      total_hospitalizations,
      threshold_met = TRUE
    )
}

find_minimum_peak_duration <- function(scenario, coarse_results) {
  scenario_results <- coarse_results |>
    filter(vax_scenario == scenario) |>
    arrange(strict_npi_days)
  
  first_coarse_hit <- scenario_results |>
    filter(peak_hospitalizations <= peak_hospitalization_threshold) |>
    slice_head(n = 1L)
  
  if (nrow(first_coarse_hit) == 0L) {
    return(tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      peak_hospitalizations = NA_real_,
      threshold_met = FALSE
    ))
  }
  
  upper_day <- first_coarse_hit$strict_npi_days[[1]]
  lower_day <- max(0L, upper_day - coarse_step + 1L)
  candidate_days <- seq.int(lower_day, upper_day)
  
  candidate_results <- map_dfr(candidate_days, function(strict_days) {
    outcomes <- simulate_duration(scenario, strict_days)
    
    tibble(
      strict_npi_days = strict_days,
      peak_hospitalizations = outcomes$peak_hospitalizations
    )
  }) |>
    filter(peak_hospitalizations <= peak_hospitalization_threshold) |>
    arrange(strict_npi_days)
  
  if (nrow(candidate_results) == 0L) {
    return(tibble(
      vax_scenario = scenario,
      minimum_strict_npi_days = NA_integer_,
      peak_hospitalizations = NA_real_,
      threshold_met = FALSE
    ))
  }
  
  candidate_results |>
    slice_head(n = 1L) |>
    transmute(
      vax_scenario = scenario,
      minimum_strict_npi_days = strict_npi_days,
      peak_hospitalizations,
      threshold_met = TRUE
    )
}

coarse_results <- run_coarse_sweep()

minimum_duration <- (future_map_dfr(names(vax_scenarios)
                                    , ~ find_minimum_duration(.x, coarse_results)) 
                     |> mutate(vax_scenario = factor(vax_scenario, levels = names(vax_scenarios))
                               , hospitalization_threshold = hospitalization_threshold
                               , margin_below_threshold = hospitalization_threshold - total_hospitalizations
                     )
)

minimum_peak_duration <- (future_map_dfr(names(vax_scenarios)
                                         , ~ find_minimum_peak_duration(.x, coarse_results)) 
                          |> mutate(vax_scenario = factor(vax_scenario, levels = names(vax_scenarios))
                                    , peak_hospitalization_threshold = peak_hospitalization_threshold
                                    , margin_below_threshold = peak_hospitalization_threshold - peak_hospitalizations
                          )
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

crossing_points <- (minimum_duration 
                    |> filter(threshold_met) 
                    |> mutate(label = paste0(minimum_strict_npi_days, " days"))
)

peak_crossing_points <- (minimum_peak_duration 
                         |> filter(threshold_met) 
                         |> mutate(label = paste0(minimum_strict_npi_days))
)

figure4A <- (
  ggplot(coarse_results, aes(x = strict_npi_days,y = total_hospitalizations,colour = vax_scenario)) 
  + geom_hline(yintercept = hospitalization_threshold,linetype = "dashed",linewidth = 2.5,colour = "red") 
  + geom_line(linewidth = 2) 
  + geom_point(size = 4) 
  + geom_segment(
    data = crossing_points,
    aes(
      x = minimum_strict_npi_days,
      xend = minimum_strict_npi_days,
      y = 1500,
      yend = hospitalization_threshold,
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
    , size = 10
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + geom_text(
    data = crossing_points,
    aes(x = minimum_strict_npi_days,y = total_hospitalizations - 1000,label = label,colour = vax_scenario)
    , angle = 90
    , vjust = -0.5
    , size = 10
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + scale_colour_manual(
    values = vax_colours,
    breaks = names(vax_scenarios),
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
    expand = expansion(mult = c(0.02, 0.12))
  ) 
  + labs(
    title = "Strict NPI duration required under different booster scenarios",
    x = "Strict NPI duration (days)",
    y = "Cumulative hospitalizations",
    tag = "A") 
  + guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 1, size = 4))) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 35, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    plot.title = element_text(size = 45, hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 35),
    legend.spacing.y = unit(0.8, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.tag = element_text(size = 45, face = "plain"),
    plot.margin = margin(12, 35, 12, 12)
  )
)


peak_y_min <- min(coarse_results$peak_hospitalizations, na.rm = TRUE)

figure4B <- (
  ggplot(coarse_results, aes(x = strict_npi_days,y = peak_hospitalizations,colour = vax_scenario)) 
  + geom_hline(yintercept = peak_hospitalization_threshold,linetype = "dashed",linewidth = 2.5,colour = "red") 
  + geom_line(linewidth = 2) 
  + geom_point(size = 4) 
  + geom_segment(
    data = peak_crossing_points,
    aes(
      x = minimum_strict_npi_days,
      xend = minimum_strict_npi_days,
      y = 0,
      yend = peak_hospitalization_threshold,
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
    , size = 10
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + geom_text(
    data = peak_crossing_points,
    aes(x = minimum_strict_npi_days,y = peak_hospitalizations - 0.18 * max(coarse_results$peak_hospitalizations, na.rm = TRUE),label = label,colour = vax_scenario)
    , angle = 90
    , vjust = -0.5
    , size = 10
    , inherit.aes = FALSE
    , show.legend = FALSE
  ) 
  + scale_colour_manual(
    values = vax_colours,
    breaks = names(vax_scenarios),
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
    expand = expansion(mult = c(0.02, 0.12))
  ) 
  + labs(
    x = "Strict NPI duration (days)",
    y = "Peak daily hospitalizations",
    tag = "B") 
  + guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 1, size = 4))) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 35, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 35),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    plot.title = element_text(size = 45, hjust = 0.5, face = "plain"),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 35),
    legend.spacing.y = unit(0.8, "cm"),
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.tag = element_text(size = 45, face = "plain"),
    plot.margin = margin(12, 35, 12, 12)
  )
)


figure4 <- (
  figure4A / figure4B
  + plot_layout(guides = "collect")
  & theme(legend.position = "bottom")
)

png("../figures/figure4.png",width = 10000,height = 7000,res = 300,bg = "white",type = "cairo")
print(figure4)
dev.off()