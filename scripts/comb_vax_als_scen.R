
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2); library(grid)
  library(ggh4x); library(scales); library(tibble)
  library(broom.mixed); library(zoo); library(patchwork)
})

set.seed(2025)
loadEnvironments()
options(dplyr.summarise.inform = FALSE)

## ------------------------- config ------------------------------------------

date0            <- as.Date("2021-12-17")
date_end         <- as.Date("2022-07-09")
offset0          <- 150
upper_plot_time  <- 300
start_date       <- date0
basis_cols       <- 5
tau_days         <- 2.0
keep_vax         <- c("0%","25%","50%","75%","100%")
sim_horizon      <- 500

fig_dir          <- "../figures"
diag_dir         <- file.path(fig_dir, "diagnostics_comb_vax_als")

csv_out          <- "../output/comb_vax_als_scen_summary.csv"
csv_out_bars     <- "../output/comb_vax_als_bars.csv"

## ------------------------- inputs ------------------------------------------

calibrator    <- rdsRead("calibrate")
multipliers   <- rdsRead("build_multipliers")
seroprevdata  <- rdsRead("seroprevdata")
fitserodata   <- rdsRead("fitsero.rds")
time_steps    <- max(fitserodata$time)

padded_base_beta     <- rdsRead("base_beta")
padded_beta0_hat     <- rdsRead("policy_neutral")
padded_beta_cf_ALS2  <- rdsRead("beta_als2")
padded_beta_cf_ALS4  <- rdsRead("beta_als4")

spec_fixed <- mp_optimized_spec(calibrator, "modified")

## ------------------------- plotting vars -----------------------------------

lim_dat   <- \(x) dplyr::filter(x, time <= upper_plot_time)
state_var <- c("S1","E1","V2","I1","A1","R1","E2","V3","I2","A2","R2","E3","I3","A3","R3")

## ------------------------- io helpers --------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

slugify <- function(x) {
  x |>
    tolower() |>
    gsub("%", "pct", x = _) |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

ensure_dir(fig_dir)
ensure_dir(diag_dir)
ensure_dir(dirname(csv_out))
ensure_dir(dirname(csv_out_bars))

## ------------------------- helpers -----------------------------------------

set_vax <- function(spec, scale) {
  if (is.null(spec$default$double_vac_values_base))
    spec$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
  if (is.null(spec$default$booster_vac_values_base))
    spec$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
  
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

simulate_series <- function(spec, scale) {
  cf <- set_vax(spec, scale)
  cf |>
    mp_simulator(sim_horizon, outputs = "inc") |>
    mp_trajectory() |>
    dplyr::filter(matrix == "inc", time >= offset0) |>
    dplyr::mutate(date = date0 + (time - offset0)) |>
    dplyr::filter(date >= date0) |>
    dplyr::transmute(date = as.Date(date), inc = value)
}

save_diagnostics <- function(spec, vax_label, als_label, beta_label) {
  p_diag <- spec |>
    mp_simulator(time_steps, outputs = state_var) |>
    mp_trajectory() |>
    lim_dat() |>
    ggplot() +
    aes(time, value) +
    geom_line() +
    facet_wrap(~matrix, scales = "free") +
    labs(
      title = paste0("Diagnostics: ", als_label, " | vax = ", vax_label, " | beta = ", beta_label),
      x     = NULL,
      y     = "Values"
    ) +
    theme_clean() +
    theme(
      axis.text.x       = element_text(size = 15, angle = 0, hjust = 1),
      axis.text.y       = element_text(size = 15),
      axis.title.y      = element_text(size = 15),
      plot.title        = element_text(size = 20, hjust = 0.5),
      strip.text        = element_text(size = 15, face = "plain"),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.background = element_blank(),
      legend.text       = element_text(size = 15),
      panel.spacing     = unit(1.5, "cm"),
      plot.background   = element_blank()
    )
  
  fn <- file.path(
    diag_dir,
    paste0(
      "diagnostic_",
      slugify(als_label), "_vax_", slugify(vax_label),
      "_beta_", slugify(beta_label),
      ".png"
    )
  )
  
  save_png(fn, p_diag, width = 5000, height = 4000)
}

## ------------------------- factors & labels --------------------------------

als_levels <- c("No-ALS Only", "ALS-4 Only")  # extend as needed
als_fac    <- factor(als_levels, levels = als_levels)

vax_scales <- c("0%" = 0, "25%" = 0.25, "50%" = 0.50, "75%" = 0.75, "100%" = 1.0)
vax_levels <- names(vax_scales)
vax_fac    <- factor(vax_levels, levels = vax_levels)

baseline_vax_scales <- setNames(rep(1.0, length(vax_levels)), vax_levels)
baseline_vax_fac    <- factor(names(baseline_vax_scales), levels = names(baseline_vax_scales))

## ------------------------- observed baseline -------------------------------

baseline <- purrr::imap_dfr(baseline_vax_scales, function(scl, vlbl) {
  simulate_series(spec_fixed, scl) |>
    dplyr::mutate(
      vax    = factor(vlbl, levels = levels(baseline_vax_fac)),
      series = "observed",
      inc    = as.integer(round(inc))
    ) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>
    dplyr::select(date, inc, vax, als_scen, series)
}) |>
  dplyr::arrange(als_scen, vax, date)

## ------------------------- replace beta(t) ---------------------------------

update_beta <- function(new_beta) {
  cp_spec <- spec_fixed
  stopifnot(length(new_beta) == cp_spec$default$time_steps)
  
  cp_spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

simulate_vax_set <- function(new_beta, als_label, beta_label) {
  stopifnot(als_label %in% levels(als_fac))
  
  spec_cf0 <- update_beta(new_beta)
  
  purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
    spec_cf <- set_vax(spec_cf0, scale_val)
    
    # diagnostics per (beta scenario, als label, vax scale)
    save_diagnostics(spec_cf, vax_label = vlab, als_label = als_label, beta_label = beta_label)
    
    spec_cf |>
      mp_simulator(sim_horizon, outputs = "inc") |>
      mp_trajectory() |>
      dplyr::filter(matrix == "inc", time >= offset0) |>
      dplyr::mutate(date = date0 + (time - offset0)) |>
      dplyr::filter(date >= date0) |>
      dplyr::transmute(date = as.Date(date), inc = as.integer(round(value))) |>
      dplyr::mutate(
        vax      = factor(vlab, levels = levels(vax_fac)),
        als_scen = factor(als_label, levels = levels(als_fac)),
        series   = "counterfactuals"
      ) |>
      dplyr::select(date, inc, vax, als_scen, series)
  })
}

beta_list <- list(
  #"Implemented ALS" = padded_base_beta,
  #"ALS-2 Only"      = padded_beta_cf_ALS2,
  "ALS-4 Only"      = padded_beta_cf_ALS4,
  "No-ALS Only"     = padded_beta0_hat
)

# if als_levels is a subset, only run those and name them consistently
beta_list <- beta_list[names(beta_list) %in% levels(als_fac)]

run_all_betas <- function(beta_list) {
  stopifnot(length(beta_list) > 0)
  
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_vax_set(beta_vec, als_label = als_label, beta_label = als_label)
  }) |>
    dplyr::arrange(als_scen, vax, date)
}

## ------------------------- simulate all scenarios --------------------------

df_all <- run_all_betas(beta_list)

vax_order <- c("0%","25%","50%","75%","100%")

panel_levels <- as.vector(outer(vax_order, als_levels, FUN = \(v, a) paste0(v, " Vax || ", a)))

plot_df <- dplyr::bind_rows(df_all, baseline) |>
  dplyr::mutate(
    date        = as.Date(date),
    vax         = factor(vax, levels = vax_order),
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = factor(paste0(vax, " Vax || ", als_scen), levels = panel_levels)
  )

## ------------------------- peak & excess summaries -------------------------

cf_peaks <- plot_df |>
  dplyr::filter(
    series == "counterfactuals",
    date   >= start_date,
    vax    %in% keep_vax,
    als_scen %in% als_levels
  ) |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx            = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date       = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups         = "drop"
  )

diff_totals <- plot_df |>
  dplyr::filter(
    date   >= start_date,
    vax    %in% keep_vax,
    als_scen %in% als_levels,
    series %in% c("counterfactuals","observed")
  ) |>
  dplyr::select(date, vax, als_scen, panel_title, series, inc) |>
  tidyr::pivot_wider(names_from = series, values_from = inc) |>
  dplyr::mutate(diff_inc = counterfactuals - observed) |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(excess_cases = sum(diff_inc, na.rm = TRUE), .groups = "drop")

cf_ann <- cf_peaks |>
  dplyr::left_join(diff_totals, by = c("vax","als_scen","panel_title")) |>
  dplyr::mutate(
    ann_lab = paste0(
      "peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "peak cases: ", scales::comma(round(peak_infections)), "\n",
      "excess cases: ", scales::comma(round(excess_cases))
    )
  ) |>
  dplyr::select(vax, als_scen, panel_title, peak_infections, peak_date, excess_cases, ann_lab)

ann_pos <- plot_df |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(
    xpos = min(date, na.rm = TRUE) + 1 + (max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) * 0.40,
    ypos = 14500,
    .groups = "drop"
  )

ann_df <- cf_ann |>
  dplyr::inner_join(ann_pos, by = c("vax","als_scen","panel_title")) |>
  dplyr::mutate(panel_title = factor(panel_title, levels = panel_levels)) |>
  droplevels()

panel_letters <- LETTERS[seq_along(panel_levels)]
panel_labels  <- tibble(
  panel_title = factor(panel_levels, levels = panel_levels),
  panel_label = panel_letters
)

plot_df <- dplyr::left_join(plot_df, panel_labels, by = "panel_title")
ann_df  <- dplyr::left_join(ann_df,  panel_labels, by = "panel_title")

## ------------------------- summary table -----------------------------------

summ_tab <- dplyr::bind_rows(df_all, baseline) |>
  dplyr::mutate(title = paste(vax, als_scen, sep = " - ")) |>
  dplyr::filter(date <= date_end) |>
  dplyr::group_by(title) |>
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups               = "drop"
  )

print(summ_tab)
write.csv(summ_tab, file = csv_out, row.names = FALSE)

## ------------------------- 20-panel time-series plot -----------------------

p <- ggplot(data = plot_df) +
  geom_line(
    data = dplyr::filter(plot_df, series == "counterfactuals"),
    aes(date, inc, color = "Counterfactual"),
    linewidth = 1.5
  ) +
  geom_point(
    data = dplyr::filter(plot_df, series == "observed"),
    aes(date, inc, color = "Actual"),
    size = 1.5
  ) +
  facet_wrap(~ panel_title, nrow = 4, ncol = 5, scales = "fixed") +
  geom_text(
    data = ann_df,
    aes(x = xpos, y = ypos, label = ann_lab),
    inherit.aes = FALSE,
    size = 6.5,
    lineheight = 0.95,
    hjust = 0,
    vjust = 1
  ) +
  geom_text(
    data = dplyr::distinct(plot_df, panel_title, panel_label),
    aes(
      x = min(plot_df$date, na.rm = TRUE),
      y = max(plot_df$inc,  na.rm = TRUE),
      label = panel_label
    ),
    hjust = -0.2,
    vjust = 1.2,
    size  = 8,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_vline(
    aes(xintercept = as.Date("2022-03-14")),
    color = "black",
    linetype  = 4,
    linewidth = 1.0
  ) +
  geom_vline(
    xintercept = as.Date("2022-05-22"),
    colour    = "purple",
    linetype  = 4,
    linewidth = 1.5
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c("Actual" = "red", "Counterfactual" = "blue"),
    name   = NULL
  ) +
  labs(
    title = "Counterfactual scenarios (blue) vs estimated infections (red)",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 22, face = "plain"),
    panel.spacing     = unit(2, "cm"),
    plot.background   = element_blank()
  )

save_png(file.path(fig_dir, "comb_vax_als_scenario.png"), p, width = 8500, height = 5800)

## ------------------------- 4-panel time-series by ALS ----------------------

cf_df <- plot_df |>
  dplyr::filter(
    series  == "counterfactuals",
    date    >= date0,
    vax     %in% vax_order,
    als_scen %in% levels(als_scen)
  )

obs_df <- plot_df |>
  dplyr::filter(series == "observed", date >= date0) |>
  dplyr::arrange(date) |>
  dplyr::group_by(date) |>
  dplyr::summarise(inc = first(inc), .groups = "drop") |>
  dplyr::mutate(series = "Actual")

global_top <- max(c(cf_df$inc, obs_df$inc), na.rm = TRUE)

label_pos <- cf_df |>
  dplyr::distinct(als_scen) |>
  dplyr::arrange(als_scen) |>
  dplyr::mutate(
    x_lab       = min(cf_df$date, na.rm = TRUE) + 3,
    y_lab       = global_top * 0.99,
    panel_label = LETTERS[seq_len(dplyr::n())]
  )

pp <- ggplot(data = plot_df) +
  geom_line(
    data = cf_df,
    aes(date, inc, color = vax),
    linewidth = 1.2,
    alpha = 1
  ) +
  geom_point(
    data = obs_df,
    aes(date, inc, color = series),
    size = 1.2,
    alpha = 1
  ) +
  facet_wrap(~ als_scen, nrow = 2, ncol = 2, scales = "fixed") +
  geom_text(
    data = label_pos,
    aes(x = x_lab, y = y_lab, label = panel_label),
    hjust = 0,
    vjust = 1,
    size  = 8,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d %Y") +
  scale_y_continuous(labels = scales::comma) +
  geom_vline(
    aes(xintercept = as.Date("2022-03-14")),
    color = "black",
    linetype  = 4,
    linewidth = 1.0
  ) +
  geom_vline(
    xintercept = as.Date("2022-05-22"),
    colour    = "purple",
    linetype  = 4,
    linewidth = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Actual" = "red",
      "0%"     = "grey50",
      "25%"    = "#9467bf",
      "50%"    = "#2ca02c",
      "75%"    = "#1f77b4",
      "100%"   = "orange"
    ),
    breaks = c("Actual","0%","25%","50%","75%","100%"),
    name   = NULL
  ) +
  labs(
    title = "Combined ALS and vaccination scenarios",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 20, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    panel.spacing     = unit(2.5, "cm"),
    plot.background   = element_blank(),
    strip.background  = element_rect(fill = "grey90", color = "grey90")
  )

save_png(file.path(fig_dir, "comb_vax_als_scen1.png"), pp, width = 5500, height = 3500)

## ------------------ stacked barplot: total infections ----------------------

cutoff_date <- as.Date("2022-03-14")

bar_df <- df_all |>
  dplyr::mutate(series = "Counterfactual") |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    vaccination = as.character(vax),
    period      = dplyr::if_else(date <= cutoff_date, "Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    vaccination = factor(vaccination, levels = c("0%","25%","50%","75%","100%")),
    als_scen    = factor(als_scen, levels = als_levels),
    period      = factor(period, levels = c("March 15 - Nov 30, 2022","Up to March 14, 2022"))
  ) |>
  dplyr::group_by(als_scen, vaccination, period) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(type = "normal")

write.csv(bar_df, csv_out_bars, row.names = FALSE)

n_panels <- length(unique(bar_df$als_scen))

panel_lab_df <- tibble::tibble(
  als_scen = levels(bar_df$als_scen)[seq_len(n_panels)]
) |>
  dplyr::mutate(
    panel_label = dplyr::case_when(
      as.character(als_scen) %in% c("No-ALS Only") ~ "E",
      as.character(als_scen) %in% c("ALS-4 Only")  ~ "H",
      TRUE ~ NA_character_
    )
  )

ordered_als_levels <- panel_lab_df |>
  dplyr::filter(!is.na(panel_label)) |>
  dplyr::arrange(panel_label) |>
  dplyr::pull(als_scen)

bar_df <- bar_df |>
  dplyr::mutate(als_scen = factor(als_scen, levels = ordered_als_levels))

panel_lab_df <- panel_lab_df |>
  dplyr::mutate(als_scen = factor(als_scen, levels = ordered_als_levels))

p_bar <- ggplot(bar_df, aes(x = vaccination, y = total_infections, fill = period)) +
  geom_col(width = 0.75, position = "stack") +
  facet_wrap(~ als_scen, nrow = 1, scales = "fixed") +
  scale_fill_manual(
    values = c(
      "Up to March 14, 2022"    = "navy",
      "March 15 - Nov 30, 2022" = "goldenrod"
    ),
    breaks = c("Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    name = NULL
  ) +
  geom_text(
    data = panel_lab_df,
    aes(x = -Inf, y = Inf, label = panel_label),
    inherit.aes = FALSE,
    hjust = 1.5,
    vjust = -0.9,
    size  = 7,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 20, 20, 80)) +
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total infections under both ALS and vaccination scenarios",
    x     = "Vaccination rate",
    y     = "Total infections (people)"
  ) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 18),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 20, face = "plain"),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    plot.background   = element_blank(),
    panel.spacing.x   = unit(2.5, "cm")
  )

save_png(file.path(fig_dir, "comb_vax_als_bar.png"), p_bar, width = 4800, height = 2000)
