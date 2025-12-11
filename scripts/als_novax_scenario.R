
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

## ------------------------- config -----------------------------------------
date0       <- as.Date("2021-12-17")
date_end    <- as.Date("2022-07-09")
start_date  <- date0
offset0     <- 150
sim_horizon <- 500

dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
csv_out <- "../output/als_novax_scenario_summary.csv"

## ------------------------- inputs -----------------------------------------
calibrator   <- rdsRead("calibrate")
multipliers  <- rdsRead("build_multipliers")
seroprevdata <- rdsRead("seroprevdata")
fitserodata <- rdsRead("fitsero.rds")
time_steps   <- max(fitserodata$time)

## time-varying betas for scenarios
padded_base_beta    <- rdsRead("base_beta")
padded_beta0_hat    <- rdsRead("policy_neutral")
padded_beta_cf_ALS2 <- rdsRead("beta_als2")
padded_beta_cf_ALS4 <- rdsRead("beta_als4")

## ------------------------- baseline spec -----------------------------------
spec_fixed <- mp_optimized_spec(calibrator, "modified")

## cache pristine vaccination series once (defensive)
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

## ------------------------- scenarios ---------------------------------------
als_levels <- c("No-ALS Only", "ALS-2 Only","Implemented ALS","ALS-4 Only")
als_fac    <- factor(als_levels, levels = als_levels)

## ------------------------- helpers (spec-aware) ----------------------------
simulate_series <- function(spec) {
  spec %>%
    mp_simulator(sim_horizon, outputs = "inc") %>%
    mp_trajectory() %>%
    dplyr::filter(matrix == "inc") %>%
    dplyr::filter(time >= offset0) %>%
    dplyr::mutate(date = date0 + (time - offset0)) %>%
    dplyr::filter(date >= date0) %>%
    dplyr::transmute(date, inc = value)
}

update_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec %>%
    mp_tmb_delete(phase = "before", at = 24L) %>%
    mp_tmb_delete(phase = "before", at = 24L) %>%
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

## ------------------------- Actual (observed) -------------------------------
actual <- simulate_series(spec_fixed) %>%
  dplyr::mutate(
    series = "Actual",
    inc    = as.integer(round(inc))
  ) %>%
  tidyr::crossing(tibble(als_scen = "Actual")) %>%
  dplyr::select(date, inc, als_scen, series) %>%
  dplyr::arrange(als_scen, date)

## ------------------------- remove vax (no vaccination) ---------------------
no_vax_spec <- spec_fixed
no_vax_spec$default$double_vac_values  <- 0 * no_vax_spec$default$double_vac_values_base
no_vax_spec$default$booster_vac_values <- 0 * no_vax_spec$default$booster_vac_values_base

## scenario-specific counterfactuals (no-vax spec; only beta differs)
simulate_set <- function(base_spec, new_beta, als_label) {
  stopifnot(als_label %in% levels(als_fac))
  spec_cf <- update_beta(base_spec, new_beta)
  simulate_series(spec_cf) %>%
    dplyr::mutate(
      date     = as.Date(date),
      inc      = as.integer(round(inc)),
      als_scen = factor(als_label, levels = levels(als_fac)),
      series   = "counterfactuals"
    ) %>%
    dplyr::select(date, inc, als_scen, series)
}

beta_list <- list(
  "Implemented ALS" = padded_base_beta,
  "ALS-2 Only"   = padded_beta_cf_ALS2,
  "ALS-4 Only"   = padded_beta_cf_ALS4,
  "No-ALS Only"  = padded_beta0_hat
)

run_all_betas <- function(base_spec, beta_list) {
  stopifnot(
    length(beta_list) > 0,
    all(names(beta_list) %in% levels(als_fac))
  )
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_set(base_spec, beta_vec, als_label)
  }) %>%
    dplyr::arrange(als_scen, date)
}

df_all <- run_all_betas(no_vax_spec, beta_list)

## ------------------------- combined plotting frame -------------------------
plot_df <- dplyr::bind_rows(df_all, actual) %>%
  dplyr::mutate(
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = factor(als_scen, levels = als_levels)
  )

## ------------------------- summaries ---------------------------------------
cf_peaks <- plot_df %>%
  dplyr::filter(series == "counterfactuals", date >= start_date) %>%
  dplyr::group_by(als_scen, panel_title) %>%
  dplyr::summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups = "drop"
  )

## excess cases (counterfactuals vs observed)
cf_only <- plot_df %>%
  dplyr::filter(series == "counterfactuals", date >= start_date) %>%
  dplyr::select(date, als_scen, panel_title, inc_cf = inc)

obs_only <- plot_df %>%
  dplyr::filter(series == "Actual", date >= start_date) %>%
  dplyr::select(date, als_scen, panel_title, inc_obs = inc)

diff_totals <- cf_only %>%
  dplyr::left_join(obs_only, by = c("date", "als_scen", "panel_title")) %>%
  dplyr::mutate(diff_inc = inc_cf - inc_obs) %>%
  dplyr::group_by(als_scen, panel_title) %>%
  dplyr::summarise(
    excess_cases = sum(diff_inc, na.rm = TRUE),
    .groups      = "drop"
  )

## annotation labels & positions
cf_ann <- cf_peaks %>%
  dplyr::left_join(diff_totals, by = c("als_scen", "panel_title")) %>%
  dplyr::mutate(
    ann_lab = paste0(
      "Peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "Peak cases: ", scales::comma(round(peak_infections))#, "\n",
      #"Excess cases: ", scales::comma(round(excess_cases))
    )
  ) %>%
  dplyr::select(als_scen, panel_title, peak_infections, peak_date, excess_cases, ann_lab)

ann_pos <- plot_df %>%
  dplyr::group_by(als_scen, panel_title) %>%
  dplyr::summarise(
    xpos = min(date, na.rm = TRUE) + 1 +
      (max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) * 0.35,
    ypos = 10000,
    .groups = "drop"
  )

ann_df <- cf_ann %>%
  dplyr::inner_join(ann_pos, by = c("als_scen", "panel_title")) %>%
  droplevels()

panel_letters <- LETTERS[seq_along(als_levels)]
panel_labels  <- tibble(
  panel_title = factor(als_levels, levels = als_levels),
  panel_label = panel_letters
)

plot_df <- dplyr::left_join(plot_df, panel_labels, by = "panel_title")
ann_df  <- dplyr::left_join(ann_df,  panel_labels, by = "panel_title")

## ------------------------- overlay figure -----------------------------------
cf_df <- plot_df %>%
  dplyr::filter(series == "counterfactuals", date >= date0) %>%
  dplyr::arrange(als_scen, date)

obs_df <- plot_df %>%
  dplyr::filter(series == "Actual", date >= date0) %>%
  dplyr::arrange(date) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(inc = dplyr::first(inc), .groups = "drop") |>
  mutate(als_scen = "Actual")

pal <- c(
  "Actual"       = "red",
  "Implemented ALS" = "purple",
  "No-ALS Only"  = "#0072B2",
  "ALS-2 Only"   = "#009E73",
  "ALS-4 Only"   = "#E69F00"
)

legend_breaks <- c("Actual", als_levels)

# ------------------------- summaries -------------------------------------
summ_tab <- dplyr::bind_rows(df_all, actual) |>
  dplyr::filter(date <= as.Date("2022-05-22")) |>
  dplyr::group_by(als_scen) |>
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups               = "drop"
  )

print(summ_tab)
write.csv(summ_tab, file = csv_out, row.names = FALSE)

## ------------------------- multi-panel time-series plot (p) ----------------

p <- ggplot() +
  geom_line(
    data = cf_df,
    aes(date, inc, color = "Counterfactual"),
    linewidth = 1.2
  ) +
  geom_point(
    data = obs_df,
    aes(date, inc, color = "Actual"),
    size = 1.2,
    linetype  = 2
  ) +
  facet_wrap(~ panel_title, nrow = 2, ncol = 2, scales = "fixed") +
  geom_text(
    data = ann_df,
    aes(x = xpos, y = ypos, label = ann_lab),
    inherit.aes = FALSE,
    size = 5.5,
    lineheight = 0.95,
    hjust = 0,
    vjust = 1
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date("2022-03-14")),
    color = "black",
    linetype  = 4,
    linewidth = 1.0
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c("Actual" = "red", "Counterfactual" = "#1b9e77"),
    name   = NULL
  ) +
  labs(
    title = "Estimated daily infections under ALS scenarios without vaccination",
    x     = NULL,
    y     = "estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    plot.title       = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 18, face = "plain"),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.background= element_blank(),
    legend.text      = element_text(size = 18),
    panel.spacing    = unit(1.5, "cm"),
    plot.background  = element_blank()
  )

png("../figures/als_novax_scenario.png",width = 4000, height = 3500, res = 300,bg = "white", type = "cairo")
print(p)
dev.off()

# ------------------------- plotting --------------------------------------

pp <- ggplot() +
  geom_line(
    data  = cf_df,
    aes(x = date, y = inc, color = als_scen, group = als_scen),
    linewidth = 1.2
  ) +
  geom_point(data = obs_df, aes(date, inc, color = "Actual"), size = 1.2, alpha = 1) +
  scale_x_date(
    expand      = c(0, 0),
    date_breaks = "2 weeks",
    date_labels = "%b %d %Y"
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = pal, breaks = legend_breaks, name = NULL) +
  labs(
    title = "Estimated daily infections under ALS scenarios without vaccination",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  geom_vline(
    aes(xintercept = as.Date("2022-03-14")),
    colour    = "black",
    linetype  = 4,
    linewidth = 1.0
  ) +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 16, hjust = 1, angle = 45),
    axis.title.x     = element_blank(),
    axis.text.y      = element_text(size = 16),
    axis.title.y     = element_text(size = 16),
    plot.title       = element_text(size = 20, hjust = 0.5, face = "plain"),
    legend.position  = c(0.75, 0.65),
    legend.direction = "vertical",
    legend.box.just  = "left",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 18, face = "plain"),
    legend.background= element_blank(),
    plot.background  = element_blank()
  )

print(pp)

png("../figures/als_novax_scenario1.png",width = 4000, height = 1800, res = 300,bg = "white", type = "cairo")
print(pp)
dev.off()

## ------------------ stacked barplot: total infections ----------------------

cutoff_date <- as.Date("2022-03-14")

cf_bar <- cf_df |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    period = dplyr::if_else(date <= cutoff_date,"Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    series  = "Counterfactual"
  ) |>
  dplyr::group_by(als_scen, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

actual_bar <- obs_df |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    period = dplyr::if_else(date <= cutoff_date,"Up to March 14, 2022","March 15 - Nov 30, 2022"),
    series  = "Actual"
  ) |>
  dplyr::group_by(als_scen, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

bar_df <- dplyr::bind_rows(cf_bar, actual_bar) |>
  dplyr::mutate(
    als_scen    = factor(als_scen, levels = c("No-ALS Only", "ALS-2 Only","Implemented ALS","ALS-4 Only","Actual")),
    period = factor(period, levels = c("March 15 - Nov 30, 2022", "Up to March 14, 2022")),fill_group = interaction(series, period, drop = TRUE)
  )

write.csv(bar_df, "../output/als_novax_scenario_bars.csv", row.names = FALSE) 

fill_vals <- c(
  "Counterfactual.Up to March 14, 2022"    = "violet",
  "Counterfactual.March 15 - Nov 30, 2022" =  "gray70",
  "Actual.Up to March 14, 2022"            = "navy",
  "Actual.March 15 - Nov 30, 2022"         = "#f55f74"
)

fill_labs <- c(
  "Actual.Up to March 14, 2022"            = "Up to March 14, 2022 (Actual)",
  "Actual.March 15 - Nov 30, 2022"         = "March 15 – Nov 30, 2022 (Actual)",
  "Counterfactual.Up to March 14, 2022"    = "Up to March 14, 2022 (Counterfactual)",
  "Counterfactual.March 15 - Nov 30, 2022" = "March 15 – Nov 30, 2022 (Counterfactual)"
)

p_bar <- ggplot() +
  geom_col(data = bar_df |> dplyr::filter(series == "Counterfactual"), 
           aes(x = als_scen, y = total_infections, fill = period), width = 0.75, position = "stack") +
  scale_fill_manual(name = "Counterfactuals",values = c("Up to March 14, 2022" = "violet","March 15 - Nov 30, 2022" = "gray70"),
                    breaks = c("Up to March 14, 2022", "March 15 - Nov 30, 2022")) +
  # reset fill scale: allows a second independent legend
  ggnewscale::new_scale_fill() +
  geom_col(
    data = bar_df |> dplyr::filter(series == "Actual"), aes(x = als_scen, y = total_infections, fill = period), width = 0.75, position = "stack") +
  scale_fill_manual(name = "Actual", values = c("Up to March 14, 2022" = "navy","March 15 - Nov 30, 2022" = "#f55f74"),
                    breaks = c("Up to March 14, 2022", "March 15 - Nov 30, 2022"))+
  
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Cumulative infections under ALS scenarios with no vaccination",
    x     = "",
    y     = "Total infections (people)"
  ) +
  
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 18),
    axis.text.y      = element_text(size = 18),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    plot.title       = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 0, face = "plain"),
    strip.background = element_blank(),
    legend.position  = "top",
    legend.direction = "vertical",
    legend.box       = "horizontal",
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 18),
    legend.background= element_blank(),
    panel.spacing.x  = unit(2, "cm")
  )

png("../figures/als_novax_bar.png",width = 4500, height = 2000, res = 300, bg = "white")
print(p_bar)
dev.off()
