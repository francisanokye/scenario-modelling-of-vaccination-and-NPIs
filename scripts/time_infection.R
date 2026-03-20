# ======================================================================
# Duration-to-target study (FAST + CORRECT)
#  - Finds the minimum number of days an ALS policy must be maintained
#    from Omicron start (date0) so that peak daily infections within
#    [date0, date_end] is <= peak_target, given observed vaccination pace.
#  - Links vaccination pace (0/50/100%) and booster timing (baseline vs -35d)
#    to policy duration and intensity (No-ALS, ALS-2, ALS-4).
#
# Key speed + correctness choices:
#   (1) Inject timeseries_beta into the TMB model ONCE per grid cell using
#       the proven update_beta() (mp_tmb_delete/insert).
#   (2) During the duration search, only overwrite spec$default$timeseries_beta
#       (fast) and simulate ONLY "inc" (fast).
#   (3) Use bracket + bisection (no coarse scan) and memoize peaks.
#   (4) After min durations are found, run a small number of full sims
#       to generate informative plots (release-date vlines, coverage, bars).
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
  library(grid)  # unit()
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

# ---- core dates and mapping ----
date0       <- as.Date("2021-12-17")  # Omicron start
date_end    <- as.Date("2022-05-22")  # evaluation end
offset0     <- 150L                   # date = date0 + (time - offset0)
sim_horizon <- 500L

# ---- target outcome ----
peak_target <- 1261

# ---- search controls ----
duration_max_days  <- 160L    # max duration from Omicron start to keep ALS policy
bisection_max_iter <- 14L     # bisection iterations

# ---- scenarios ----
als_levels <- c("No-ALS Only", "ALS-2 Only", "ALS-4 Only")

# vaccination pace scenarios (keep 0/50/100 only)
vax_scales <- c("0%" = 0.0, "50%" = 0.5, "100%" = 1.0)
vax_order  <- names(vax_scales)

# booster timing scenarios: baseline vs early (Nov 8 = 35 days earlier than Dec 13)
booster_shift_days <- c("Baseline boosters" = 0L, "Early boosters (−35d)" = 35L)

# ---- paths ----
fig_dir <- "../figures"
out_dir <- "../output"
dur_dir <- file.path(fig_dir, "als_duration_policy_linked_to_vax")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dur_dir, recursive = TRUE, showWarnings = FALSE)

csv_summary <- file.path(out_dir, "als_duration_summary.csv")
csv_detail  <- file.path(out_dir, "als_duration_detail.csv")
csv_report  <- file.path(out_dir, "als_duration_report_table.csv")

save_png <- function(path, plot_obj, width, height, res = 260, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

# ---- inputs ----
calibrator  <- rdsRead("calibrate")
fitserodata <- rdsRead("fitsero.rds")
time_steps  <- max(fitserodata$time)

beta_noals <- rdsRead("policy_neutral")
beta_als2  <- rdsRead("beta_als2")
beta_als4  <- rdsRead("beta_als4")

beta_policy_list <- list(
  "No-ALS Only" = beta_noals,
  "ALS-2 Only"  = beta_als2,
  "ALS-4 Only"  = beta_als4
)

# ---- base spec ----
spec_fixed <- mp_optimized_spec(calibrator, "modified")
spec_fixed$default$time_steps <- time_steps

# ---- cache pristine vaccination series once ----
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

stopifnot(
  !is.null(spec_fixed$default$double_vac_values_base),
  !is.null(spec_fixed$default$booster_vac_values_base),
  !is.null(spec_fixed$integers$booster_vac_changepoints_base)
)

# ---- index mapping ----
idx_start <- as.integer(offset0)

# ---- proven beta updater (KEEP) ----
update_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

# ---- scale both dose-2 and boosters (pace/coverage scaling) ----
set_vax_scale <- function(spec, scale = 1.0) {
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

# ---- shift booster schedule earlier by N days (keeps values, shifts changepoints) ----
make_early_booster <- function(spec, shift_days = 0L) {
  cf <- spec
  if (shift_days == 0L) return(cf)
  
  cp0  <- cf$integers$booster_vac_changepoints_base
  val0 <- cf$default$booster_vac_values_base
  stopifnot(length(cp0) == length(val0), length(cp0) > 0)
  
  cp <- pmax(1L, as.integer(cp0) - as.integer(shift_days))
  o  <- order(cp)
  
  cp_s  <- cp[o]
  val_s <- as.numeric(val0)[o]
  
  # keep last if duplicates (clipping to 1 can duplicate)
  keep_last <- !duplicated(cp_s, fromLast = TRUE)
  cp_s  <- cp_s[keep_last]
  val_s <- val_s[keep_last]
  
  cf$integers$booster_vac_changepoints <- cp_s
  cf$default$booster_vac_values        <- val_s
  
  # set these as the new "base" for downstream scaling
  cf$integers$booster_vac_changepoints_base <- cp_s
  cf$default$booster_vac_values_base        <- val_s
  
  cf
}

# Construct beta(t) such that:
#  - policy beta is applied from date0 for duration_days
#  - after that, beta switches to No-ALS beta
compose_beta_duration <- function(beta_policy, beta_release, duration_days) {
  beta_cf <- beta_release
  
  start_i <- idx_start
  end_i   <- min(length(beta_cf), start_i + as.integer(duration_days) - 1L)
  
  if (end_i >= start_i) {
    beta_cf[start_i:end_i] <- beta_policy[start_i:end_i]
  }
  beta_cf
}

# =============================================================================
# FAST + CORRECT SEARCH ENGINE
# =============================================================================

# Peak simulator for "inc" only
peak_from_spec <- function(spec) {
  tr <- spec |>
    mp_simulator(sim_horizon, outputs = "inc") |>
    mp_trajectory()
  
  tr |>
    mutate(date = date0 + (time - offset0)) |>
    filter(date >= date0, date <= date_end) |>
    summarise(pk = max(value, na.rm = TRUE), .groups = "drop") |>
    pull(pk)
}

# Memoised peak evaluator for a given scenario template (with timeseries_beta injected!)
make_peak_evaluator <- function(spec_template_injected, beta_policy) {
  cache <- new.env(parent = emptyenv())
  
  function(duration_days) {
    key <- as.character(duration_days)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    
    beta_cf <- compose_beta_duration(beta_policy, beta_noals, duration_days)
    
    # overwrite ONLY the injected numeric vector (fast, correct)
    spec_run <- spec_template_injected
    spec_run$default$timeseries_beta <- as.numeric(beta_cf)
    
    pk <- peak_from_spec(spec_run)
    
    assign(key, pk, envir = cache)
    pk
  }
}

# Bracket + bisection: find minimum duration achieving peak <= target
find_min_duration_fast <- function(spec_template_injected, beta_policy) {
  
  peak_at <- make_peak_evaluator(spec_template_injected, beta_policy)
  
  # check 1 day
  pk1 <- peak_at(1L)
  if (pk1 <= peak_target) {
    det <- tibble(duration_days = 1L,
                  peak_daily_infections = pk1,
                  meets_target = TRUE)
    return(list(min_days = 1L, detail = det))
  }
  
  # expanding bracket
  hi <- 2L
  pk_hi <- peak_at(hi)
  while (hi < duration_max_days && pk_hi > peak_target) {
    hi <- min(duration_max_days, hi * 2L)
    pk_hi <- peak_at(hi)
    if (hi == duration_max_days) break
  }
  
  # not achievable
  if (pk_hi > peak_target) {
    sample_days <- unique(pmin(duration_max_days, c(1L,2L,4L,8L,16L,32L,64L,128L,duration_max_days)))
    det <- tibble(duration_days = sample_days) |>
      mutate(peak_daily_infections = map_dbl(duration_days, peak_at),
             meets_target = peak_daily_infections <= peak_target)
    return(list(min_days = NA_integer_, detail = det))
  }
  
  lo <- as.integer(floor(hi / 2L))
  iter <- 0L
  while (lo + 1L < hi && iter < bisection_max_iter) {
    mid <- as.integer(floor((lo + hi) / 2L))
    pk_mid <- peak_at(mid)
    if (pk_mid <= peak_target) hi <- mid else lo <- mid
    iter <- iter + 1L
  }
  
  # detail curve: sparse + dense near solution
  dense  <- unique(pmax(1L, pmin(duration_max_days, hi + (-10:10))))
  sparse <- unique(pmin(duration_max_days, c(1L,2L,4L,8L,16L,32L,64L,128L,duration_max_days)))
  days <- sort(unique(c(sparse, dense)))
  
  det <- tibble(duration_days = days) |>
    mutate(peak_daily_infections = map_dbl(duration_days, peak_at),
           meets_target = peak_daily_infections <= peak_target)
  
  list(min_days = hi, detail = det)
}

# ---- run grid across policy intensity (ALS), vaccination pace, booster timing ----
run_duration_study_fast <- function() {
  
  grid <- tidyr::crossing(
    als_scen  = als_levels,
    vax_lab   = names(vax_scales),
    boost_lab = names(booster_shift_days)
  ) |>
    mutate(
      vax_scale  = unname(vax_scales[vax_lab]),
      shift_days = unname(booster_shift_days[boost_lab])
    )
  
  res <- pmap(grid, function(als_scen, vax_lab, boost_lab, vax_scale, shift_days) {
    
    # Build scenario template spec ONCE per grid cell (vax + boosters)
    spec_tmp <- spec_fixed |>
      make_early_booster(shift_days = shift_days) |>
      set_vax_scale(scale = vax_scale)
    
    # CRITICAL: inject timeseries_beta into TMB ONCE (correctness)
    spec_tmp_injected <- update_beta(spec_tmp, beta_noals)
    
    beta_policy <- beta_policy_list[[als_scen]]
    out <- find_min_duration_fast(spec_tmp_injected, beta_policy)
    
    detail <- out$detail |>
      mutate(
        als_scen = als_scen,
        vax = factor(vax_lab, levels = vax_order),
        boosters = boost_lab,
        vax_scale = vax_scale,
        booster_shift_days = shift_days
      )
    
    tibble(
      als_scen = als_scen,
      vax = factor(vax_lab, levels = vax_order),
      boosters = boost_lab,
      vax_scale = vax_scale,
      booster_shift_days = shift_days,
      min_duration_days = out$min_days
    ) |>
      mutate(detail = list(detail))
  })
  
  summary <- bind_rows(res) |>
    mutate(
      als_scen = factor(als_scen, levels = als_levels),
      boosters = factor(boosters, levels = names(booster_shift_days)),
      achieved = !is.na(min_duration_days)
    ) |>
    arrange(boosters, vax, als_scen)
  
  detail <- bind_rows(map(res, "detail")) |>
    mutate(
      als_scen = factor(als_scen, levels = als_levels),
      boosters = factor(boosters, levels = names(booster_shift_days))
    ) |>
    arrange(boosters, vax, als_scen, duration_days)
  
  list(summary = summary, detail = detail)
}

# =============================================================================
# RUN STUDY
# =============================================================================
dur_out <- run_duration_study_fast()
summary_df <- dur_out$summary
detail_df  <- dur_out$detail

write.csv(summary_df, csv_summary, row.names = FALSE)
write.csv(detail_df,  csv_detail,  row.names = FALSE)

print(summary_df)

# =============================================================================
# PLOTS (5) + REPORT TABLE
# =============================================================================
max_d <- duration_max_days

# Helper for plotting NA (not achieved)
summary_plot_df <- summary_df |>
  mutate(
    bar_height = if_else(is.na(min_duration_days), as.numeric(max_d), as.numeric(min_duration_days)),
    bar_label  = if_else(is.na(min_duration_days), paste0("> ", max_d, " d"), paste0(min_duration_days, " d"))
  )

# ---------------- Plot 1: Peak vs duration curves (per ALS), with min-day vline ----------
main_vax   <- "100%"
main_boost <- "Baseline boosters"

main_detail <- detail_df |>
  filter(vax == main_vax, boosters == main_boost)

main_vlines <- summary_df |>
  filter(vax == main_vax, boosters == main_boost, !is.na(min_duration_days)) |>
  transmute(als_scen, xint = min_duration_days)

p1 <- ggplot(main_detail, aes(duration_days, peak_daily_infections)) +
  geom_hline(yintercept = peak_target, linetype = 4, linewidth = 1.0) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(shape = meets_target), size = 2.2) +
  geom_vline(data = main_vlines, aes(xintercept = xint), linetype = 3, linewidth = 1.0) +
  facet_wrap(~als_scen, nrow = 1, scales = "free_y") +
  scale_shape_manual(values = c(`FALSE`=1, `TRUE`=16), name = NULL,
                     labels = c("Does not meet target", "Meets target")) +
  scale_x_continuous(breaks = pretty_breaks(7), limits = c(1, max_d)) +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste0("Minimum duration from Omicron start to keep peak daily infections \u2264 ", comma(peak_target)),
    subtitle = paste0("Observed vaccination pace (100%) | ", main_boost),
    x = "Days ALS is maintained from Dec 17, 2021 (then switches to No-ALS)",
    y = "Peak daily infections (Dec 17, 2021–May 22, 2022)"
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size=22, hjust=0.5, face="plain"),
    plot.subtitle = element_text(size=16, hjust=0.5),
    strip.text = element_text(size=18),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16),
    legend.position = "top"
  )

save_png(file.path(dur_dir, "P1_peak_vs_duration_observedpace.png"), p1, 5600, 2100, res=300)

# ---------------- Plot 2: Heatmap of minimum duration by coverage (vax) x ALS ----------
p2 <- ggplot(summary_df, aes(vax, als_scen, fill = min_duration_days)) +
  geom_tile(color = "white", linewidth = 0.8) +
  facet_wrap(~boosters, nrow = 1) +
  scale_fill_gradient(low = "grey90", high = "grey20", na.value = "grey95") +
  geom_text(aes(label = if_else(is.na(min_duration_days), paste0("> ", max_d), as.character(min_duration_days))),
            size = 5) +
  labs(
    title = "Minimum ALS duration required to meet the target peak (days)",
    subtitle = paste0("Target: peak daily infections \u2264 ", comma(peak_target)),
    x = "Vaccination pace (scaled)", y = NULL, fill = "Days"
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size=22, hjust=0.5, face="plain"),
    plot.subtitle = element_text(size=16, hjust=0.5),
    strip.text = element_text(size=16),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16),
    legend.position = "right"
  )

save_png(file.path(dur_dir, "P2_heatmap_min_duration.png"), p2, 5200, 1800, res=300)

# ---------------- Plot 3: Bar chart (min duration) grouped by vax within boosters ----------
p3 <- ggplot(summary_plot_df, aes(als_scen, bar_height)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = bar_label), vjust = -0.4, size = 5) +
  facet_grid(boosters ~ vax) +
  scale_y_continuous(breaks = pretty_breaks(7), limits = c(0, max_d*1.15)) +
  labs(
    title = "Minimum duration needed (or not achieved) across vaccination pace and booster timing",
    x = NULL, y = "Minimum duration (days)"
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size=22, hjust=0.5, face="plain"),
    strip.text = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=16)
  )

save_png(file.path(dur_dir, "P3_bar_min_duration_grid.png"), p3, 5200, 2600, res=300)

# ---------------- Plot 4: Time series at the minimum duration (observed pace) ----------
traj_inc_only <- function(spec) {
  spec |>
    mp_simulator(sim_horizon, outputs = "inc") |>
    mp_trajectory() |>
    mutate(date = date0 + (time - offset0)) |>
    filter(date >= date0, date <= date_end) |>
    transmute(date = as.Date(date), inc = as.numeric(value))
}

simulate_at_min <- function(als_scen, min_days) {
  d_use <- ifelse(is.na(min_days), max_d, min_days)
  
  spec_tmp <- spec_fixed |>
    make_early_booster(shift_days = 0L) |>
    set_vax_scale(scale = 1.0)
  
  # inject baseline first (correctness)
  spec_tmp <- update_beta(spec_tmp, beta_noals)
  
  beta_policy <- beta_policy_list[[als_scen]]
  beta_cf <- compose_beta_duration(beta_policy, beta_noals, d_use)
  
  spec_tmp$default$timeseries_beta <- as.numeric(beta_cf)
  
  inc <- traj_inc_only(spec_tmp) |>
    mutate(
      als_scen = als_scen,
      duration_days = d_use,
      switch_date = date0 + (d_use - 1L)
    )
  inc
}

min_obs <- summary_df |>
  filter(vax == "100%", boosters == "Baseline boosters") |>
  select(als_scen, min_duration_days)

ts_df <- map2_dfr(min_obs$als_scen, min_obs$min_duration_days, simulate_at_min) |>
  mutate(als_scen = factor(als_scen, levels = als_levels))

switch_df <- ts_df |> distinct(als_scen, switch_date, duration_days)

p4 <- ggplot(ts_df, aes(date, inc)) +
  geom_hline(yintercept = peak_target, linetype = 4, linewidth = 1.0) +
  geom_line(linewidth = 1.2) +
  geom_vline(data = switch_df, aes(xintercept = switch_date), linetype = 3, linewidth = 1.0) +
  facet_wrap(~als_scen, nrow = 1, scales = "free_y") +
  scale_x_date(expand = c(0,0), date_breaks="2 weeks", date_labels="%b %d") +
  scale_y_continuous(labels=comma) +
  labs(
    title = "Daily infections under the minimum required duration (observed vaccination pace)",
    subtitle = "Vertical dashed line indicates when the policy is released to No-ALS",
    x = NULL, y = "Daily infections"
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size=22, hjust=0.5, face="plain"),
    plot.subtitle = element_text(size=16, hjust=0.5),
    strip.text = element_text(size=18),
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=16)
  )

save_png(file.path(dur_dir, "P4_timeseries_at_min_duration_observed.png"), p4, 5600, 2100, res=300)

# ---------------- Plot 5: Vaccination coverage proxy at switch date (dose-2 + boosters) ----------
traj_cov <- function(spec) {
  spec |>
    mp_simulator(sim_horizon, outputs = c("double_vac","booster_shot")) |>
    mp_trajectory() |>
    mutate(date = date0 + (time - offset0)) |>
    filter(date >= date0, date <= date_end) |>
    transmute(date = as.Date(date), matrix, value = as.numeric(value))
}

extract_cov_at_switch <- function(als_scen, min_days) {
  d_use <- ifelse(is.na(min_days), max_d, min_days)
  
  spec_tmp <- spec_fixed |>
    make_early_booster(shift_days = 0L) |>
    set_vax_scale(scale = 1.0)
  
  # inject baseline then overwrite beta vector
  spec_tmp <- update_beta(spec_tmp, beta_noals)
  
  beta_policy <- beta_policy_list[[als_scen]]
  beta_cf <- compose_beta_duration(beta_policy, beta_noals, d_use)
  spec_tmp$default$timeseries_beta <- as.numeric(beta_cf)
  
  tr <- traj_cov(spec_tmp)
  
  switch_date <- date0 + (d_use - 1L)
  
  out <- tr |>
    filter(date == switch_date) |>
    select(matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value)
  
  tibble(
    als_scen = als_scen,
    duration_days = d_use,
    switch_date = switch_date,
    double_vac = out$double_vac,
    booster_shot = out$booster_shot
  )
}

cov_df <- map2_dfr(min_obs$als_scen, min_obs$min_duration_days, extract_cov_at_switch) |>
  mutate(als_scen = factor(als_scen, levels = als_levels)) |>
  pivot_longer(cols = c(double_vac, booster_shot),
               names_to = "coverage_type", values_to = "coverage_value") |>
  mutate(
    coverage_type = recode(coverage_type,
                           double_vac = "Dose-2 vaccination (model input proxy)",
                           booster_shot = "Booster uptake (model input proxy)")
  )

p5 <- ggplot(cov_df, aes(als_scen, coverage_value, fill = coverage_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, na.rm = TRUE) +
  geom_text(aes(label = comma(round(coverage_value))),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 4) +
  labs(
    title = "Vaccination coverage proxies at the policy-release day (observed pace)",
    subtitle = "These are the vaccination inputs used by the calibrated model at the switch day",
    x = NULL, y = "Value (model scale)", fill = NULL
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size=22, hjust=0.5, face="plain"),
    plot.subtitle = element_text(size=16, hjust=0.5),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=16),
    legend.position = "top",
    legend.text = element_text(size=12)
  )

save_png(file.path(dur_dir, "P5_coverage_at_switch_observed.png"), p5, 5200, 2200, res=300)

# ---------------- Report table ----------------
report_tbl <- summary_df |>
  mutate(
    target = peak_target,
    achieved = !is.na(min_duration_days),
    min_duration_days = if_else(is.na(min_duration_days), duration_max_days, min_duration_days),
    release_date = date0 + (min_duration_days - 1L)
  ) |>
  select(boosters, vax, als_scen, achieved, min_duration_days, release_date, target) |>
  arrange(boosters, vax, als_scen)

print(report_tbl)
write.csv(report_tbl, csv_report, row.names = FALSE)
