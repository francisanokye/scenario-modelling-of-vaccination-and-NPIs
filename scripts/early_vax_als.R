
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# config 
# dates, output dir, and time offset used to align model time with calendar time
date0      <- as.Date("2021-12-15")
date_end   <- as.Date("2022-05-22")
start_date <- date0
dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
offset0 <- 150

# inputs
# calibrated model, multipliers, seroprevalence, and time horizon
calibrator    <- rdsRead("calibrate")
multipliers   <- rdsRead("build_multipliers")
seroprevdata  <- rdsRead("seroprevdata")
time_steps    <- max(seroprevdata$time)

# time-varying beta series for scenarios
padded_base_beta <- rdsRead("base_beta")
padded_beta0_hat <- rdsRead("policy_neutral")

# baseline spec 
# optimized spec with fitted parameters
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# cache pristine vaccination series once (defensive)
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

# scenarios 
# vaccination scale factors and labels
vax_scales <- c("25%" = 0.25, "50%" = 0.50, "75%" = 0.75)
vax_fac    <- factor(names(vax_scales), levels = names(vax_scales))

# als levels for facets
als_levels <- c("Observed ALS", "No-ALS Only")
als_fac    <- factor(als_levels, levels = als_levels)

# scale daily vaccination rates
set_vax <- function(spec, scale) {
  if (is.null(spec$default$double_vac_values_base))
    spec$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
  if (is.null(spec$default$booster_vac_values_base))
    spec$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

# replace beta(t) in a spec
update_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

# simulate a single series at a given vaccination scale
simulate_series <- function(spec, scale) {
  cf <- set_vax(spec, scale)
  cf |>
    mp_simulator(time_steps, outputs = "inc") |>
    mp_trajectory() |>
    filter(matrix == "inc") |>
    mutate(date = date0 + (time - offset0)) |>
    filter(date >= date0, date <= date_end) |>
    transmute(date, inc = value)
}

# simulate across vaccination scales for one beta(t) and tag with its als label
simulate_vax_set <- function(spec, new_beta, als_label) {
  stopifnot(als_label %in% levels(als_fac))
  spec_cf <- update_beta(spec, new_beta)
  purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
    simulate_series(spec_cf, scale_val) |>
      mutate(
        date     = as.Date(date),
        inc      = as.integer(round(inc)),
        vax      = factor(vlab, levels = levels(vax_fac)),
        als_scen = factor(als_label, levels = levels(als_fac)),
        series   = "counterfactuals"
      ) |>
      select(date, inc, vax, als_scen, series)
  })
}

# observed baseline replicated across panels (vaccination kept at 100% for all vax labels)
baseline_vax_scales <- c("25%" = 1.0, "50%" = 1.0, "75%" = 1.0)
baseline_vax_levels <- names(baseline_vax_scales)
baseline_vax_fac    <- factor(baseline_vax_levels, levels = baseline_vax_levels)

baseline <- purrr::imap_dfr(baseline_vax_scales, function(scl, vlbl) {
  simulate_series(spec_fixed, scl) |>
    mutate(
      vax    = factor(vlbl, levels = levels(baseline_vax_fac)),
      series = "observed",
      inc    = as.integer(inc)
    ) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>
    select(date, inc, vax, als_scen, series)
}) |> arrange(als_scen)

# list of beta(t) by als scenario
beta_list <- list(
  "Observed ALS" = padded_base_beta,
  "No-ALS Only"  = padded_beta0_hat
)

# run all beta scenarios for a given spec
run_all_betas <- function(base_spec, beta_list) {
  stopifnot(length(beta_list) > 0, all(names(beta_list) %in% levels(als_fac)))
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_vax_set(base_spec, beta_vec, als_label)
  }) |> arrange(als_scen, date)
}

# early booster shots (-28 days) 
# clone calibrated spec and shift booster changepoints earlier by 28 days
cf_early28 <- spec_fixed
cp  <- as.integer(cf_early28$integers$booster_vac_changepoints_base) - 28L
val <- as.numeric(cf_early28$default$booster_vac_values_base)

# clamp negatives to 0 (cannot start before t = 0)
cp <- pmax(0L, cp)

# sort changepoints and corresponding values
order_index <- order(cp)
cp  <- cp[order_index]
val <- val[order_index]

# write back into the scenario spec
cf_early28$integers$booster_vac_changepoints <- cp
cf_early28$default$booster_vac_values        <- val

# simulate all scenarios 
df_all <- run_all_betas(cf_early28, beta_list)

# combined plotting frame
plot_df <- bind_rows(df_all, baseline) |>
  mutate(
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = factor(paste0(als_scen), levels = als_levels)
  )

# summaries (peaks and excess cases) 
cf_peaks <- plot_df |>
  filter(series == "counterfactuals", date >= start_date, date <= date_end) |>
  group_by(als_scen, panel_title) |>
  summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx            = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date       = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups         = "drop"
  )

cf_only <- plot_df |>
  filter(series == "counterfactuals", date >= start_date, date <= date_end) |>
  select(date, als_scen, panel_title, inc_cf = inc)

obs_only <- plot_df |>
  filter(series == "observed", date >= start_date, date <= date_end) |>
  select(date, als_scen, panel_title, inc_obs = inc)

diff_totals <- cf_only |>
  left_join(obs_only, by = c("date", "als_scen", "panel_title")) |>
  mutate(diff_inc = inc_cf - inc_obs) |>
  group_by(als_scen, panel_title) |>
  summarise(excess_cases = sum(diff_inc, na.rm = TRUE), .groups = "drop")

cf_ann <- cf_peaks |>
  left_join(diff_totals, by = c("als_scen", "panel_title")) |>
  mutate(
    ann_lab = paste0(
      "peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "peak cases: ", comma(round(peak_infections)), "\n",
      "excess cases: ", comma(round(excess_cases))
    )
  ) |>
  select(als_scen, panel_title, peak_infections, peak_date, excess_cases, ann_lab)

ann_pos <- plot_df |>
  group_by(als_scen, panel_title) |>
  summarise(
    xpos = min(date, na.rm = TRUE) + 1 + (max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) * 0.1,
    ypos = 4000,
    .groups = "drop"
  )

ann_df <- cf_ann |>
  inner_join(ann_pos, by = c("als_scen", "panel_title")) |>
  droplevels()

panel_letters <- LETTERS[seq_along(als_levels)]
panel_labels  <- tibble(
  panel_title = factor(als_levels, levels = als_levels),
  panel_label = panel_letters
)

plot_df <- left_join(plot_df, panel_labels, by = "panel_title")
ann_df  <- left_join(ann_df,  panel_labels, by = "panel_title")

# three panels for manuscript 
cf_df <- plot_df |>
  filter(
    series == "counterfactuals",
    date >= date0, date <= date_end,
    vax %in% c("25%", "50%", "75%"),
    als_scen %in% levels(als_scen)
  )

# single observed baseline line (one per date)
obs_df <- plot_df |>
  filter(series == "observed", date >= date0, date <= date_end) |>
  arrange(date) |>
  group_by(date) |>
  summarise(inc = first(inc), .groups = "drop") |>
  mutate(series = "Baseline")

global_top <- max(c(cf_df$inc, obs_df$inc), na.rm = TRUE)

# panel labels at a common y across facets
label_pos <- cf_df |>
  distinct(als_scen) |>
  arrange(als_scen) |>
  mutate(
    x_lab        = min(cf_df$date, na.rm = TRUE) + 3,
    y_lab        = global_top * 0.99,
    panel_label  = LETTERS[seq_len(n())]
  )

pp <- ggplot(data = plot_df) +
  geom_line(data = cf_df, aes(date, inc, color = vax), linewidth = 0.5, alpha = 1) +
  geom_point(data = obs_df, aes(date, inc, color = series), size = 0.25, alpha = 1) +
  facet_wrap(~ als_scen, nrow = 1, ncol = 2, axes = "all", scales = "fixed") +
  geom_text(data = label_pos, aes(x = x_lab, y = y_lab, label = panel_label),hjust = 0, vjust = 1, size = 3, fontface = "bold", inherit.aes = FALSE) +
  scale_x_date(limits = c(date0, date_end), expand = c(0, 0),date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c("Baseline" = "red", "25%" = "#9467bd", "50%" = "#2ca02c", "75%" = "#1f77b4"),
    breaks = c("Baseline", "25%", "50%", "75%"),
    name   = NULL
  ) +
  labs(
    title = "Early boosters with/without ALS and varying vaccination rates",
    x = NULL, y = "estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x    = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.y    = element_text(size = 5),
    axis.title.y   = element_text(size = 5),
    plot.title     = element_text(size = 5, hjust = 0.5, face = "plain"),
    strip.text     = element_text(size = 5, face = "plain", hjust = 0.5),
    legend.position = "top",
    legend.background = element_blank(),
    legend.text      = element_text(size = 5),
    panel.spacing    = unit(0.7, "lines"),
    plot.background  = element_blank(),
    strip.background = element_blank()
  )

png("../figures/early_vax_als.png",width = 1300, height = 800, res = 300, bg = "white", type = "cairo")
print(pp)
dev.off()

