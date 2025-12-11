suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

## --------------------------- config ---------------------------------------

date0       <- as.Date("2021-12-17")
date_end    <- as.Date("2022-07-09")
start_date <- date0
offset0    <- 150                      # align model time with calendar time
sim_horizon <- 500
csv_out    <- "../output/early_vax_als_summary.csv"

dir.create("../figures", recursive = TRUE, showWarnings = FALSE)

## --------------------------- inputs ---------------------------------------

calibrator   <- rdsRead("calibrate")
multipliers  <- rdsRead("build_multipliers")
seroprevdata <- rdsRead("seroprevdata")
fitserodata <- rdsRead("fitsero.rds")
time_steps   <- max(fitserodata$time)

# time-varying beta series
padded_base_beta    <- rdsRead("base_beta")
padded_beta0_hat    <- rdsRead("policy_neutral")
padded_beta_cf_ALS2 <- rdsRead("beta_als2")
padded_beta_cf_ALS4 <- rdsRead("beta_als4")

# optimized spec with fitted parameters
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# cache pristine vaccination series once (defensive)
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

## ------------------------- scenario factors -------------------------------

# vaccination scale factors and labels
vax_scales <- c("0%" = 0, "25%" = 0.25, "50%" = 0.50, "75%" = 0.75, "100%" = 1.0)
vax_fac    <- factor(names(vax_scales), levels = names(vax_scales))
vax_order  <- levels(vax_fac)

# als levels (facet order)
als_levels <- c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")
als_fac    <- factor(als_levels, levels = als_levels)

## ---------------------- vaccination & beta -----------------------

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
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

# simulate a single series at a given vaccination scale
simulate_series <- function(spec, scale) {
  cf <- set_vax(spec, scale)
  cf %>%
    mp_simulator(sim_horizon, outputs = "inc") %>%
    mp_trajectory() %>%
    dplyr::filter(matrix == "inc") %>%
    dplyr::filter(time >= offset0) %>%
    dplyr::mutate(date = date0 + (time - offset0)) %>%
    dplyr::filter(date >= date0) %>%
    dplyr::transmute(date, inc = value)
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

## ---------------------- observed (actual) series ------------------------

# replicate observed trajectory across all vax levels and als panels
actual_vax_scales <- c("0%" = 1.0, "25%" = 1.0, "50%" = 1.0, "75%" = 1.0, "100%" = 1.0)
actual_vax_fac    <- factor(names(actual_vax_scales), levels = names(actual_vax_scales))

actual <- purrr::imap_dfr(actual_vax_scales, function(scl, vlbl) {
  simulate_series(spec_fixed, scl) |>
    mutate(
      vax    = factor(vlbl, levels = levels(actual_vax_fac)),
      series = "observed",
      inc    = as.integer(inc)
    ) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>
    select(date, inc, vax, als_scen, series)
}) |>
  arrange(als_scen, date)

## ---------------------- beta list by ALS scenario -------------------------

beta_list <- list(
  "Implemented ALS" = padded_base_beta,
  "ALS-2 Only"      = padded_beta_cf_ALS2,
  "ALS-4 Only"      = padded_beta_cf_ALS4,
  "No-ALS Only"     = padded_beta0_hat
)

run_all_betas <- function(base_spec, beta_list) {
  stopifnot(length(beta_list) > 0, all(names(beta_list) %in% levels(als_fac)))
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_vax_set(base_spec, beta_vec, als_label)
  }) |>
    arrange(als_scen, date)
}

##--------------------------- early boosters --------------------------------

# clone calibrated spec and shift booster changepoints earlier by 28 days
cf_early28 <- spec_fixed
cp  <- as.integer(cf_early28$integers$booster_vac_changepoints_base) - 28L
val <- as.numeric(cf_early28$default$booster_vac_values_base)

# clamp negatives to 0 (cannot start before t = 0) and sort
cp           <- pmax(0L, cp)
order_index  <- order(cp)
cp           <- cp[order_index]
val          <- val[order_index]

cf_early28$integers$booster_vac_changepoints <- cp
cf_early28$default$booster_vac_values        <- val

##---------------------- simulate all early-booster scenarios & merge with actual ---------------
df_all <- run_all_betas(cf_early28, beta_list)

panel_levels <- c(
  "0% Vax & No-ALS Only",
  "25% Vax & No-ALS Only",
  "50% Vax & No-ALS Only",
  "75% Vax & No-ALS Only",
  "100% Vax & No-ALS Only",
  "0% Vax & ALS-2 Only",
  "25% Vax & ALS-2 Only",
  "50% Vax & ALS-2 Only",
  "75% Vax & ALS-2 Only",
  "100% Vax & ALS-2 Only",
  "0% Vax & Implemented ALS",
  "25% Vax & Implemented ALS",
  "50% Vax & Implemented ALS",
  "75% Vax & Implemented ALS",
  "100% Vax & Implemented ALS",
  "0% Vax & ALS-4 Only",
  "25% Vax & ALS-4 Only",
  "50% Vax & ALS-4 Only",
  "75% Vax & ALS-4 Only",
  "100% Vax & ALS-4 Only"
)

plot_df <- bind_rows(df_all, actual) |>
  mutate(
    date        = as.Date(date),
    vax         = factor(vax, levels = vax_order),
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = paste0(vax, " Vax & ", als_scen),
    panel_title = factor(panel_title, levels = panel_levels)
  )

## ------------------------- peak & excess summaries ------------------------

cf_peaks <- plot_df |>
  filter(series == "counterfactuals", date >= start_date) |>
  group_by(als_scen, panel_title) |>
  summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx            = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date       = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups         = "drop"
  )

cf_only <- plot_df |>
  filter(series == "counterfactuals", date >= start_date) |>
  select(date, als_scen, panel_title, inc_cf = inc)

obs_only <- plot_df |>
  filter(series == "observed", date >= start_date) |>
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
    xpos = min(date, na.rm = TRUE) + 1 +
      (max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) * 0.35,
    ypos = 14500,
    .groups = "drop"
  )

ann_df <- cf_ann |>
  inner_join(ann_pos, by = c("als_scen", "panel_title")) |>
  droplevels()

## --------------------- panel letters (20-panel figure) --------------------

panel_letters <- LETTERS[seq_along(panel_levels)]

panel_labels <- tibble(
  panel_title = factor(panel_levels, levels = panel_levels),
  panel_label = panel_letters
)

plot_df <- left_join(plot_df, panel_labels, by = "panel_title")
ann_df  <- left_join(ann_df,  panel_labels, by = "panel_title")

## --------------------------- summary table --------------------------------

summ_tab <- bind_rows(df_all, actual) |>
  mutate(title = paste(vax, als_scen, sep = " - ")) |>
  filter(date <= date_end) |>
  group_by(title) |>
  summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups               = "drop"
  )

print(summ_tab)
write.csv(summ_tab, file = csv_out, row.names = FALSE)

## ---------------------- Figure Sx: 20-panel early-booster grid ------------

p <- ggplot(data = plot_df) +
  geom_line(data = filter(plot_df, series == "counterfactuals"),aes(date, inc, color = "Counterfactual"),linewidth = 1.5) +
  geom_point(data = filter(plot_df, series == "observed"),aes(date, inc, color = "Actual"),size = 1.5) +
  facet_wrap(~ panel_title, nrow = 4, ncol = 5, scales = "fixed") +
  geom_text(data = ann_df,aes(x = xpos, y = ypos, label = ann_lab),inherit.aes = FALSE,size = 7.0, lineheight = 0.95, hjust = 0, vjust = 1) +
  geom_text(data = distinct(plot_df, panel_title, panel_label),aes(x = min(plot_df$date, na.rm = TRUE),
    y = max(plot_df$inc,  na.rm = TRUE), label = panel_label), hjust = -0.2, vjust = 1.2, size = 8, fontface = "bold",inherit.aes = FALSE) +
  geom_vline(
    aes(xintercept = as.Date("2022-03-14")),
    color = "black", linetype  = 4, linewidth = 1.0
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c("Actual" = "red", "Counterfactual" = "#2ca02c"),
    name   = NULL
  ) +
  labs(
    title = "Early booster schedule: counterfactual infections vs estimated baseline",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x   = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 18),
    axis.title.y  = element_text(size = 18),
    plot.title    = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text    = element_text(size = 22, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box = "horizontal",
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 22, face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    plot.background   = element_blank()
  )

png("../figures/early_vax_als_scenario.png",width = 8000, height = 5800, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

## ---------------------- Figure X: 4-panel ALS × vax curves ----------------

# counterfactual series only
cf_df <- plot_df |>
  filter(
    series == "counterfactuals",
    date   >= date0,
    vax    %in% c("0%", "25%", "50%", "75%", "100%"),
    als_scen %in% levels(als_scen)
  )

# single observed baseline line
obs_df <- plot_df |>
  filter(series == "observed", date >= date0) |>
  arrange(date) |>
  group_by(date) |>
  summarise(inc = first(inc), .groups = "drop") |>
  mutate(series = "Actual")

global_top <- max(c(cf_df$inc, obs_df$inc), na.rm = TRUE)

# panel labels (A–D) by ALS scenario
label_pos <- cf_df |>
  distinct(als_scen) |>
  arrange(als_scen) |>
  mutate(
    x_lab       = min(cf_df$date, na.rm = TRUE) + 3,
    y_lab       = global_top * 0.99,
    panel_label = LETTERS[seq_len(n())]
  )

pp <- ggplot(data = plot_df) +
  geom_line(data = cf_df,aes(date, inc, color = vax),linewidth = 1.2, alpha = 1) +
  geom_point(data = obs_df,aes(date, inc, color = series),size = 1.2, alpha = 1) +
  facet_wrap(~ als_scen, nrow = 2, ncol = 2, scales = "fixed") +
  geom_text(data = label_pos,aes(x = x_lab, y = y_lab, label = panel_label),
    hjust = 0, vjust = 1, size = 8, fontface = "bold",inherit.aes = FALSE) +
  geom_vline(aes(xintercept = as.Date("2022-03-14")),color = "black", linetype  = 4, linewidth = 1.0) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d %Y") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c(
      "Actual" = "red",
      "0%"     = "gray",
      "25%"    = "#9467bf",
      "50%"    = "#2ca02c",
      "75%"    = "#1f77b4",
      "100%"   = "orange"
    ),
    breaks = c("Actual", "0%", "25%", "50%", "75%", "100%"),
    name   = NULL
  ) +
  labs(
    title = "Early booster schedule: infections by ALS scenario and vaccination rate",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x   = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 1),
    axis.title.y  = element_text(size = 20),
    plot.title    = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text    = element_text(size = 22, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box = "horizontal",
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    panel.spacing     = unit(2.5, "cm"),
    plot.background   = element_blank(),
    strip.background  = element_blank()
  )

png("../figures/early_vax_als1.png",width = 5500, height = 3500, res = 300, bg = "white", type = "cairo")
print(pp)
dev.off()

## ------------------ stacked bar: total infections by period ---------------

cutoff_date <- as.Date("2022-03-14")

bar_df <- df_all |>
  mutate(series = "Counterfactual") |>
  filter(date >= date0, date <= date_end) |>
  mutate(
    vaccination = if_else(series == "Actual", "Actual", as.character(vax)),
    period      = if_else(date <= cutoff_date,"Up to March 14, 2022","March 15 - Nov 30, 2022"),
    vaccination = factor(vaccination, levels = c("Actual", "0%", "25%", "50%", "75%", "100%")),
    als_scen    = factor(als_scen, levels = als_levels),
    period      = factor(period, levels = c("March 15 - Nov 30, 2022","Up to March 14, 2022"))
        ) |>
  group_by(als_scen, vaccination, period) |>
  summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop") |>
  mutate(type = "early")

write.csv(bar_df, "../output/early_boosters_bars.csv", row.names = FALSE)

p_bar <- ggplot(bar_df, aes(x = vaccination, y = total_infections, fill = period)) +
  geom_col(width = 0.75, position = "stack") +
  facet_wrap(~ als_scen, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c("Up to March 14, 2022" = "brown","March 15 - Nov 30, 2022" = "gold"),name = NULL) +
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(labels = comma) +
  labs(title = "Total infections under both ALS and early booster scenarios", x = "Vaccination rate", y = "Total infections (people)") +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 17, angle = 0, hjust = 0.5),
    axis.text.y      = element_text(size = 17),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    plot.title       = element_text(size = 20, hjust = 0.5, face = "plain"),
    strip.text       = element_text(size = 18, face = "plain"),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.box.just  = "left",
    legend.background = element_blank(),
    legend.text      = element_text(size = 16),
    plot.background  = element_blank(),
    panel.spacing.x  = unit(2, "cm")
  )

# number of panels (als_scen facets)
n_panels <- length(unique(bar_df$als_scen))

# panel letters starting at E (5th letter)
panel_letters <- LETTERS[5:(5 + n_panels - 1)]

# label dataframe, keeping facet order
panel_lab_df <- tibble::tibble(
  als_scen     = levels(bar_df$als_scen)[seq_len(n_panels)],
  panel_label  = panel_letters
)

p_bar <- p_bar +
  geom_text(
    data = panel_lab_df,
    aes(
      x = -Inf,   
      y = Inf,    
      label = panel_label
    ),
    inherit.aes = FALSE,
    hjust = -0.6,   # push left outside panel
    vjust = 1.4,    # slightly below top border
    size  = 7,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +       # allow drawing outside panel region
  theme(
    plot.margin = margin(20, 20, 20, 80)  # extra left space for labels
  )

png("../figures/early_vax_als_bar.png",width = 4500, height = 2000, res = 300, bg = "white", type = "cairo")
print(p_bar)
dev.off()
