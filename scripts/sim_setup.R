suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
  library(forcats)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

calibrator <- rdsRead("calibrate")

padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

# ---- helper functions ----
save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

panel_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

# ---- scenario colours ----
scenario_cols <- c(
  "National Rate"   = "springgreen3",
  "Slower Boosters" = "magenta",
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "blue",
  "Faster Boosters" = "orange"
)

scenario_dose_cols <- list(
  "Actual Boosters" = c("≤2 Doses" = "grey70",   "Booster" = "#4D4D4D"),
  "Slower Boosters" = c("≤2 Doses" = "#E8B6D2",  "Booster" = "magenta"),
  "Faster Boosters" = c("≤2 Doses" = "#FDBB84",  "Booster" = "orange"),
  "Early Boosters"  = c("≤2 Doses" = "#9ECAE1",  "Booster" = "blue")
)

booster_start_dates <- c(
  "Actual Boosters" = booster_start,
  "Early Boosters"  = early_start_date,
  "Faster Boosters" = booster_start,
  "Slower Boosters" = booster_start,
  "National Rate"   = booster_start
)

# ---- build vax specs ----
spec_actual <- mp_optimized_spec(calibrator, "modified")

slower_booster_changepoints <- as.integer(vax_slower_rate$day)
slower_booster_values       <- as.numeric(vax_slower_rate$daily_rate)
if (!any(slower_booster_changepoints == 0)) {
  slower_booster_changepoints <- c(0L, slower_booster_changepoints)
  slower_booster_values       <- c(0,  slower_booster_values)
}
spec_slower <- mp_optimized_spec(calibrator, "modified")
spec_slower$integers$booster_vac_changepoints <- as.integer(slower_booster_changepoints)
spec_slower$default$booster_vac_values        <- as.numeric(slower_booster_values)

faster_booster_changepoints <- as.integer(vax_faster_rate$day)
faster_booster_values       <- as.numeric(vax_faster_rate$daily_rate)
if (!any(faster_booster_changepoints == 0)) {
  faster_booster_changepoints <- c(0L, faster_booster_changepoints)
  faster_booster_values       <- c(0,  faster_booster_values)
}
spec_faster <- mp_optimized_spec(calibrator, "modified")
spec_faster$integers$booster_vac_changepoints <- as.integer(faster_booster_changepoints)
spec_faster$default$booster_vac_values        <- as.numeric(faster_booster_values)

early_booster_changepoints <- as.integer(vax_early_rate$day)
early_booster_values       <- as.numeric(vax_early_rate$daily_rate)
if (!any(early_booster_changepoints == 0)) {
  early_booster_changepoints <- c(0L, early_booster_changepoints)
  early_booster_values       <- c(0,  early_booster_values)
}
spec_early <- mp_optimized_spec(calibrator, "modified")
spec_early$integers$booster_vac_changepoints <- as.integer(early_booster_changepoints)
spec_early$default$booster_vac_values        <- as.numeric(early_booster_values)

national_booster_changepoints <- as.integer(can_boosterframe$day)
national_booster_values       <- as.numeric(can_boosterframe$daily_rate)
if (!any(national_booster_changepoints == 0)) {
  national_booster_changepoints <- c(0L, national_booster_changepoints)
  national_booster_values       <- c(0,  national_booster_values)
}
spec_natl <- mp_optimized_spec(calibrator, "modified")
spec_natl$integers$booster_vac_changepoints <- as.integer(national_booster_changepoints)
spec_natl$default$booster_vac_values        <- as.numeric(national_booster_values)

# ---- beta updater ----
update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

# ---- scenario lists and ordering ----
vax_scenario_specs <- list(
  "National Rate"   = spec_natl,
  "Slower Boosters" = spec_slower,
  "Actual Boosters" = spec_actual,
  "Early Boosters"  = spec_early,
  "Faster Boosters" = spec_faster
)

policy_list <- list(
  "Actual ALS"  = padded_beta_base,
  "ALS-4 Only"  = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "National Rate",
  "Slower Boosters",
  "Actual Boosters",
  "Early Boosters",
  "Faster Boosters"
)

model_outputs <- c(
  "incS", "incv2", "incv3", "inc",
  "inc_leq2dose", "inc_booster",
  "cum_inc_leq2dose", "cum_inc_booster", "cum_inc_total",
  "hosp_leq2dose", "hosp_booster", "hosp_total",
  "death_leq2dose", "death_booster", "death_total",
  "cum_hosp_leq2dose", "cum_hosp_booster", "cum_hosp_total",
  "cum_death_leq2dose", "cum_death_booster", "cum_death_total",
  "cum_booster"
)

diagnostic_outputs <- c(
  "S1", "V2", "V3", "E1", "E2", "E3", "A1", "A2",
  "A3", "I1", "I2", "I3", "R1", "R2", "R3", "inc",
  "inc_leq2dose", "inc_booster", "p1", "p2", "p3",
  "cum_booster", "cum_double_vac"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name    = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

saveVars(
  fig_dir,
  save_png, max_or_na, panel_tag,
  scenario_cols, scenario_dose_cols, booster_start_dates,
  update_vax_spec_beta,
  vax_scenario_specs, policy_list,
  policy_order, vax_order,
  model_outputs, diagnostic_outputs,
  sim_grid
)