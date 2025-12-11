
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2); library(grid)
  library(scales); library(tibble)
  library(broom.mixed); library(zoo); library(patchwork)
})

set.seed(2025)
loadEnvironments()
options(dplyr.summarise.inform = FALSE)

## ------------------------- config ------------------------------------------

date0       <- as.Date("2021-12-17")
date_end    <- as.Date("2022-07-09")
offset0         <- 150
upper_plot_time <- 300
start_date      <- date0
sim_horizon     <- 500

csv_out_summary <- "../output/vax_only_scenario_summary.csv"
csv_out_bars    <- "../output/vax_only_scenario_bars.csv"

## ------------------------- inputs ------------------------------------------

calibrator   <- rdsRead("calibrate")
multipliers  <- rdsRead("build_multipliers")
seroprevdata <- rdsRead("seroprevdata")
fitserodata <- rdsRead("fitsero.rds")
time_steps   <- max(fitserodata$time)

padded_beta0_hat <- rdsRead("policy_neutral")
spec_fixed       <- mp_optimized_spec(calibrator, "modified")

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

simulate_series <- function(spec, scale = 1) {
  cf <- set_vax(spec, scale)
  cf |>
    mp_simulator(sim_horizon, outputs = "inc") |>
    mp_trajectory() |>
    dplyr::filter(matrix == "inc", time >= offset0) |>
    dplyr::mutate(date = date0 + (time - offset0)) |>
    dplyr::filter(date >= date0) |>
    dplyr::transmute(date = as.Date(date), inc = value)
}

update_beta <- function(new_beta) {
  cp_spec <- spec_fixed
  stopifnot(length(new_beta) == cp_spec$default$time_steps)
  cp_spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

## ------------------------- factors & labels --------------------------------

als_levels <- c("No-ALS Only")
vax_scales <- c("0%" = 0, "25%" = 0.25, "50%" = 0.50, "75%" = 0.75, "100%" = 1.0)

vax_cf_levels <- names(vax_scales)
keep_vax      <- vax_cf_levels

panel_levels <- paste0(vax_cf_levels, " Vax & No-ALS Only")

## ------------------------- actual (observed) -------------------------------

actual <- simulate_series(spec_fixed) |>
  dplyr::mutate(
    series = "Actual",
    inc    = as.integer(round(inc))
  )

## ------------------------- simulate counterfactuals ------------------------

simulate_vax_set <- function(new_beta, als_label) {
  spec_cf <- update_beta(new_beta)
  purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
    simulate_series(spec_cf, scale_val) |>
      dplyr::mutate(
        date       = as.Date(date),
        inc        = as.integer(round(inc)),
        vax        = factor(vlab, levels = vax_cf_levels),
        als_scen   = factor(als_label, levels = als_levels),
        series     = "Counterfactual",
        panel_title = factor(
          paste0(vax, " Vax & ", als_scen),
          levels = panel_levels
        )
      ) |>
      dplyr::select(date, inc, vax, als_scen, series, panel_title)
  })
}

beta_list <- list("No-ALS Only" = padded_beta0_hat)

df_all <- purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
  simulate_vax_set(beta_vec, als_label)
}) |>
  dplyr::arrange(als_scen, vax, date)

cf_df_full <- df_all

## ------------------------- actual replicated to all panels -----------------

actual_panel_df <- actual |>
  dplyr::select(date, inc, series) |>
  tidyr::crossing(tibble(
    vax      = factor(vax_cf_levels, levels = vax_cf_levels),
    als_scen = factor(als_levels, levels = als_levels)
  )) |>
  dplyr::mutate(
    panel_title = factor(
      paste0(vax, " Vax & ", als_scen),
      levels = panel_levels
    )
  )

## ------------------------- peak & excess summaries -------------------------

actual_series <- actual |>
  dplyr::select(date, inc) |>
  dplyr::rename(actual_inc = inc)

cf_peaks <- cf_df_full |>
  dplyr::filter(
    date     >= start_date,
    vax      %in% keep_vax,
    als_scen %in% als_levels
  ) |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx            = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date       = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups         = "drop"
  )

diff_totals <- cf_df_full |>
  dplyr::filter(
    date     >= start_date,
    vax      %in% keep_vax,
    als_scen %in% als_levels
  ) |>
  dplyr::left_join(actual_series, by = "date") |>
  dplyr::mutate(diff_inc = inc - actual_inc) |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(
    excess_cases = sum(diff_inc, na.rm = TRUE),
    .groups      = "drop"
  )

cf_ann <- cf_peaks |>
  dplyr::left_join(diff_totals, by = c("vax", "als_scen", "panel_title")) |>
  dplyr::mutate(
    ann_lab = paste0(
      "peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "peak cases: ", scales::comma(round(peak_infections)), "\n",
      "excess cases: ", scales::comma(round(excess_cases))
    )
  )

ann_pos <- cf_df_full |>
  dplyr::group_by(vax, als_scen, panel_title) |>
  dplyr::summarise(
    xpos = min(date, na.rm = TRUE) +
      1 + (max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) * 0.45,
    ypos = max(inc, na.rm = TRUE) * 0.9,
    .groups = "drop"
  )

ann_df <- cf_ann |>
  dplyr::inner_join(ann_pos, by = c("vax", "als_scen", "panel_title"))

## ------------------------- summary table -----------------------------------

summ_df <- dplyr::bind_rows(
  cf_df_full |> dplyr::mutate(series = "Counterfactual"),
  actual |>
    tidyr::crossing(tibble(
      vax      = "Actual",
      als_scen = "No-ALS Only"
    )) |>
    dplyr::mutate(series = "Actual")
)

summ_tab <- summ_df |>
  dplyr::mutate(
    title = paste(vax, als_scen, sep = " - ")
  ) |>
  dplyr::filter(date <= date_end) |>
  dplyr::group_by(title, series) |>
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups               = "drop"
  )

print(summ_tab)
write.csv(summ_tab, file = csv_out_summary, row.names = FALSE)

## ------------------------- multi-panel time-series plot (p) ----------------

p <- ggplot() +
  geom_line(
    data = cf_df_full,
    aes(date, inc, color = "Counterfactual"),
    linewidth = 1.2
  ) +
  geom_point(
    data = actual_panel_df,
    aes(date, inc, color = "Actual"),
    size = 1.2,
    linetype  = 2
  ) +
  facet_wrap(~ panel_title, nrow = 3, ncol = 2, scales = "fixed") +
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
    title = "Estimated daily infections under vaccination-only scenario without NPIs (No-ALS)",
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

png("../figures/vax_scenarios.png",width = 4000, height = 3500, res = 300,bg = "white", type = "cairo")
print(p)
dev.off()

#####---------------------------------------------------------------------------------------

pal <- c(
  "Actual"       = "red",
  "0%" = "gray49",
  "25%" = "purple",
  "50%"  = "#0072B2",
  "75%"   = "#009E73",
  "100%"   = "#E69F00"
)

legend_breaks <- c("Actual",vax_cf_levels)

pp <- ggplot() +
  geom_line(
    data  = cf_df_full,
    aes(x = date, y = inc, color = vax, group = vax),
    linewidth = 1.2
  ) +
  geom_point(
    data = actual_panel_df,
    aes(date, inc, color = "Actual"),
    size = 1.2,
    linetype  = 2
  ) +
  scale_x_date(
    expand      = c(0, 0),
    date_breaks = "2 weeks",
    date_labels = "%b %d %Y"
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = pal, breaks = legend_breaks, name = NULL) +
  labs(
    title = "Estimated daily infections under vaccination-only scenario without NPIs (No-ALS)",
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
    legend.position  = c(0.8, 0.75),
    legend.direction = "horizontal",
    legend.box.just  = "left",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 18, face = "plain"),
    legend.background= element_blank(),
    plot.background  = element_blank()
  )

print(pp)

png("../figures/vax_scenario1.png",width = 4000, height = 1800, res = 300,bg = "white", type = "cairo")
print(pp)
dev.off()

## ------------------ stacked barplot: total infections ----------------------

cutoff_date <- as.Date("2022-03-14")

cf_bar <- df_all |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    vaccination = as.character(vax),
    period = dplyr::if_else(date <= cutoff_date,"Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    series  = "Counterfactual",
    als_scen = "No-ALS Only"
  ) |>
  dplyr::group_by(als_scen, vaccination, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

actual_bar <- actual |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    vaccination = "Actual",
    period = dplyr::if_else(date <= cutoff_date,"Up to March 14, 2022","March 15 - Nov 30, 2022"),
    series  = "Actual",
    als_scen = "No-ALS Only"
  ) |>
  dplyr::group_by(als_scen, vaccination, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

bar_df <- dplyr::bind_rows(cf_bar, actual_bar) |>
  dplyr::mutate(
    als_scen    = factor(als_scen, levels = "No-ALS Only"),
    vaccination = factor(vaccination, levels = c("0%", "25%", "50%", "75%", "100%", "Actual")),
    period = factor(period, levels = c("March 15 - Nov 30, 2022", "Up to March 14, 2022")),fill_group = interaction(series, period, drop = TRUE)
    )

write.csv(bar_df, csv_out_bars, row.names = FALSE)

fill_vals <- c(
  "Counterfactual.Up to March 14, 2022"    = "#1b9e77",
  "Counterfactual.March 15 - Nov 30, 2022" =  "gold",
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
           aes(x = vaccination, y = total_infections, fill = period), width = 0.75, position = "stack") +
  scale_fill_manual(name = "Counterfactuals",values = c("Up to March 14, 2022" = "#1b9e77","March 15 - Nov 30, 2022" = "gold"),
    breaks = c("Up to March 14, 2022", "March 15 - Nov 30, 2022")) +
  # reset fill scale: allows a second independent legend
  ggnewscale::new_scale_fill() +
  geom_col(
    data = bar_df |> dplyr::filter(series == "Actual"), aes(x = vaccination, y = total_infections, fill = period), width = 0.75, position = "stack") +
  scale_fill_manual(name = "Actual", values = c("Up to March 14, 2022" = "navy","March 15 - Nov 30, 2022" = "#f55f74"),
                    breaks = c("Up to March 14, 2022", "March 15 - Nov 30, 2022"))+
  facet_wrap(~ als_scen, nrow = 1, scales = "fixed") +
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(labels = scales::comma) +
  
  labs(
    title = "cumulative infections under vaccination-only scenario without ALS",
    x     = "vaccination rate",
    y     = "total infections (people)"
  ) +
  
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 18),
    axis.text.y      = element_text(size = 18),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    plot.title       = element_text(size = 24, hjust = 0.5, face = "plain"),
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

png("../figures/vax_scenarios_bar.png",width = 4500, height = 2000, res = 300, bg = "white", type = "cairo")
print(p_bar)
dev.off()


