# Scenario question: How long must Strict NPI (strict NPIs) remain in place, per vaccination scenario,
# to keep peak concurrent hospitalizations below a threshold?
#
# Metric: peak concurrent hospitalizations (bed occupancy). flows.R has no hospital
# compartment, so occupancy is reconstructed from daily admissions (hosp_total) with a
# mean length of stay L:  H_t = (1 - 1/L) H_{t-1} + admissions_t; peak = max(H_t).

# Length of stay is uncertain, so the result is not tied to one value: admissions are
# simulated once, occupancy is reconstructed for a RANGE of L, and the plot shows the
# central curve with a length-of-stay band. 

library(macpan2)
library(tidyverse)
library(shellpipes)
library(ggthemes)
library(scales)

loadEnvironments()

calib     <- rdsRead("calibrate")
base_beta <- rdsRead("base")

strict_npi_constant_value <- 0.245

# observed N.L. peak concurrent hospitalizations used to anchor the scale (BA.2 wave; 25 for BA.1)
hosp_threshold   <- 47                
baseline_peak_hosp <- 47                
baseline_scenario  <- "Actual Boosters" # real-world baseline used to anchor

# length-of-stay SENSITIVITY
los_central <- 13
los_grid    <- c(7, 10, 13, 16, 20)

# Strict NPI durations in days: 0, 7, 14, ... (weekly), capped at the beta series length.
max_duration <- length(base_beta) - as.integer(step_start) + 1L
durations    <- seq(0L, 180L, by = 7L)
durations    <- durations[durations <= max_duration]

# vaccination scenarios 
vax_scenarios <- c(
  "Actual Boosters" = "vax_actual_boosters"
  , "Early Boosters"  = "vax_early_boosters"
  , "Faster Boosters" = "vax_faster_rate"
  , "Slower Boosters" = "vax_slower_rate"
  , "National Rate"   = "vax_national_rate"
)

load_booster <- function(stem) {
  e <- new.env(); load(paste0(stem, ".rda"), envir = e); get("new_booster", envir = e)
}

booster_list <- lapply(vax_scenarios, load_booster)

# hospital occupancy from admissions 
# concurrent occupancy H_t = (1 - 1/L) H_{t-1} + admissions_t; peak = max(H_t)
occ_peak <- function(admissions, los) {
  max(as.numeric(stats::filter(admissions, filter = 1 - 1 / los, method = "recursive", init = 0)))
}

# model specification 
base_spec <- (mp_optimized_spec(calib, "modified")
              |> mp_tmb_delete(phase = "before", at = 26L)
)

make_beta <- function(duration) {
  beta <- base_beta
  if (duration >= 1) {
    end <- min(as.integer(step_start) + duration - 1L, length(beta))
    beta[as.integer(step_start):end] <- strict_npi_constant_value
  }
  beta
}

# simulate one (scenario, duration) and return the daily admissions series
admissions_one <- function(new_booster, duration) {
  new_beta <- make_beta(duration)
  spec <- base_spec |> mp_tmb_insert(default = list(timeseries_beta = new_beta))
  cp <- as.integer(new_booster$day);
  vv <- as.numeric(new_booster$daily_rate)
  if (!any(cp == 0)) { cp <- c(0L, cp); vv <- c(0, vv) }
  spec$integers$booster_vac_changepoints <- as.integer(cp)
  spec$default$booster_vac_values        <- as.numeric(vv)
  
  w <- (spec
        |> mp_simulator(time_steps = sim_horizon, outputs = c("hosp_total"))
        |> mp_trajectory()
        |> arrange(time)
  )
  w$value
}

# simulate admissions once per scenario x duration 
sim_adm <- map_dfr(names(vax_scenarios), function(scenario) {
  booster <- booster_list[[scenario]]
  map_dfr(durations, function(days) {
    tibble(vax_scenario = scenario, duration = days
           , adm = list(admissions_one(booster, days)))
  })
})

# occupancy for each length of stay
sweep_los <- (map_dfr(los_grid, function(L) {
  sim_adm |> mutate(los = L, peak_hosp = map_dbl(adm, occ_peak, los = L))
})
|> select(-adm)
|> group_by(los)
|> mutate(peak_hosp = peak_hosp *
            (baseline_peak_hosp / peak_hosp[vax_scenario == baseline_scenario & duration == 0]))
|> ungroup()
|> mutate(vax_scenario = factor(vax_scenario, levels = names(vax_scenarios)))
)

# central curve + length-of-stay band (min-max across L)
sweep_central <- sweep_los |> filter(los == los_central)
sweep_band <- (sweep_los
               |> group_by(vax_scenario, duration)
               |> summarise(lo = min(peak_hosp), hi = max(peak_hosp), .groups = "drop")
)
plot_df <- sweep_central |> left_join(sweep_band, by = c("vax_scenario", "duration"))

# Strict NPI duration where each scenario first meets the threshold ----
duration_to_meet <- function(x, y, threshold) {
  o <- order(x); x <- x[o]; y <- y[o]
  if (all(y <= threshold)) return(min(x))
  if (all(y >  threshold)) return(NA_real_)
  j <- which(y <= threshold)[1]
  if (j == 1) return(x[1])
  x1 <- x[j - 1]; x2 <- x[j]; y1 <- y[j - 1]; y2 <- y[j]
  x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
}

# required duration for every L, then summarise to central + range
req_by_los <- (sweep_los
               |> group_by(vax_scenario, los)
               |> summarise(days = duration_to_meet(duration, peak_hosp, hosp_threshold), .groups = "drop")
)

required_tbl <- (req_by_los
                 |> group_by(vax_scenario)
                 |> summarise(
                   hosp_days_central = round(days[los == los_central])
                   , hosp_days_lo    = round(suppressWarnings(min(days, na.rm = TRUE)))
                   , hosp_days_hi    = round(suppressWarnings(max(days, na.rm = TRUE)))
                   , .groups = "drop"
                 )
                 |> mutate(across(c(hosp_days_lo, hosp_days_hi), ~ ifelse(is.infinite(.x), NA_real_, .x)))
)

# peak concurrent hospitalizations vs Strict NPI duration ----
cross_df <- (required_tbl
             |> filter(!is.na(hosp_days_central))
             |> mutate(strict_npi_days = hosp_days_central
                       , cross_lab = paste0("~", round(strict_npi_days))
                       , lab_y = hosp_threshold * 0.20
             )
)

vax_cols <- c(
  "Actual Boosters" = "#4d4d4d"
  , "Early Boosters" = "#e6550d"
  , "Faster Boosters" = "#3182bd"
  , "Slower Boosters" = "#d6699f"
  , "National Rate" = "springgreen3"
)

p <- (
  ggplot(plot_df, aes(x = duration, y = peak_hosp, colour = vax_scenario))
  + geom_ribbon(aes(ymin = lo, ymax = hi, fill = vax_scenario), alpha = 0.15, colour = NA)
  + geom_line(linewidth = 2.0)
  + geom_point(size = 4)
  + geom_hline(yintercept = hosp_threshold, colour = "red", linetype = "dashed", linewidth = 1.6)
  + geom_segment(data = cross_df, aes(x = hosp_days_lo, xend = hosp_days_hi, y = hosp_threshold, yend = hosp_threshold, colour = vax_scenario)
                 , linewidth = 1.4
                 , inherit.aes = FALSE
  )
  + geom_segment(data = cross_df, aes(x = strict_npi_days, xend = strict_npi_days, y = 0, yend = hosp_threshold, colour = vax_scenario)
                 , linetype = "dashed"
                 , linewidth = 0.7
                 , inherit.aes = FALSE
  )
  + geom_point(data = cross_df, aes(x = strict_npi_days, y = hosp_threshold, colour = vax_scenario), shape = 8, size = 10, stroke = 2, inherit.aes = FALSE)
  + geom_text(data = cross_df, aes(x = strict_npi_days, y = lab_y, label = cross_lab, colour = vax_scenario)
              , angle = 90
              , vjust = -0.3
              , size = 8
              , inherit.aes = FALSE
              , show.legend = FALSE)
  + scale_colour_manual(values = vax_cols, name = NULL, drop = FALSE)
  + scale_fill_manual(values = vax_cols, name = NULL, drop = FALSE, guide = "none")
  + scale_x_continuous(breaks = seq(0, 180, by = 14))
  + scale_y_continuous(labels = comma)
  + labs(
    title = "Peak hospitalizations across strict NPI duration and vaccination scenarios",
    x = "Strict NPI duration (days)", y = "Peak hospitalizations")
  + guides(colour = guide_legend(override.aes = list(shape = 8, linetype = 0, size = 4)))
  + theme_clean()
  + theme(
    axis.text.x = element_text(size = 30, hjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
    legend.position = c(0.90, 0.80),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 30),
    legend.spacing.y = unit(1.5, "cm"),        
    legend.key.spacing.y = unit(1.0, "cm"),    
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )
)

png("../figures/sensitivity_npi_duration_plot.png", width = 6000, height = 3600, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

rdsSave(sweep_los)
