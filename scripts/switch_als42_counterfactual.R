# counterfactuals via policy-neutral transmission scaling with robust injection/fallback
# all comments are lowercase and concise

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2)
  library(ggh4x); library(scales); library(tibble);
  library(tidyverse); library(broom.mixed); library(zoo)
})

set.seed(2025)
loadEnvironments()
options(dplyr.summarise.inform = FALSE)

date0    <- as.Date("2021-12-15")
date_end <- as.Date("2022-05-22")
offset0  <- 150
upper_plot_time <- 300
start_date <- date0
basis_cols <- 5
tau_days <- 2.0
keep_vax   <- c("25%","50%","75%")

# load saved inputs
calibrator <- rdsRead("calibrate")
multipliers   <- rdsRead("build_multipliers")
seroprevdata <- rdsRead("seroprevdata")
time_steps   <- max(seroprevdata$time)

padded_beta_ALS_ALS2 <- rdsRead("switch_als4_als2")

# initialize calibrated model spec with optimized parameters
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# beta_hat(t) from calibrator and extract beta
beta_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  dplyr::filter(matrix == "beta_thing", time >= offset0, time <= upper_plot_time) |>
  mutate(date = seq(start_date, by = "day", length.out = length(value)), beta_hat = value) |>
  filter(date >= date0, date <= date_end) |>
  select(date, beta_hat)

# scale daily c=vaccination rates
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
  cf <- spec %>% set_vax(scale)
  cf %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    filter(matrix == "inc") %>%
    mutate(date = date0 + (time - offset0)) %>%
    filter(date >= date0, date <= date_end) %>%
    transmute(date, inc = value)
}

# ALS scenarios required
als_cols <- c("ALS-2 ==> ALS-2")
als_fac  <- factor(als_cols, levels = als_cols)

# use baseline rate (observed @ 100%) and replicate infection values across vaccination scales
baseline_vax_scales <- c("25%" = 1.0, "50%" = 1.0, "75%" = 1.0)
baseline_vax_levels <- names(baseline_vax_scales)
baseline_vax_fac    <- factor(baseline_vax_levels, levels = baseline_vax_levels)

baseline <- purrr::imap_dfr(baseline_vax_scales, function(scl, vlbl) {
  simulate_series(spec_fixed, scl) |>
    dplyr::mutate(
      vax = factor(vlbl, levels = levels(baseline_vax_fac)),
      series = "observed"
    ) |>
    mutate(inc = as.integer(inc)) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>
    dplyr::select(date, inc, vax, als_scen, series)
}) |>
  arrange(als_scen)

####################################### Counterfactual Scenarios #######################################

vax_scales <- c("25%" = 0.25, "50%" = 0.50, "75%" = 0.75)
vax_levels <- names(vax_scales)
vax_fac    <- factor(vax_levels, levels = vax_levels)

# factors (if not already defined)
vax_scales <- c("25%"=0.25,"50%"=0.50,"75%"=0.75)
vax_fac    <- factor(names(vax_scales), levels = names(vax_scales))
als_levels <- c("ALS-2 ==> ALS-2")
als_fac    <- factor(als_levels, levels = als_levels)

update_beta <- function(new_beta){
  cp_spec <- spec_fixed 
  stopifnot(length(new_beta) == cp_spec$default$time_steps)
  cp_spec = cp_spec |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
  return(cp_spec)
}

## simulate for one updated beta and tag its ALS label
simulate_vax_set <- function(new_beta, als_label) {
  stopifnot(als_label %in% levels(als_fac))
  spec_cf <- update_beta(new_beta)
  purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
    simulate_series(spec_cf, scale_val) |>
      dplyr::mutate(
        date     = as.Date(date),
        inc      = as.integer(round(inc)),
        vax      = factor(vlab, levels = levels(vax_fac)),
        als_scen = factor(als_label, levels = levels(als_fac)),
        series   = "counterfactuals"
      ) |>
      dplyr::select(date, inc, vax, als_scen, series)
  })
}

beta_list <- list("ALS-2 ==> ALS-2" = padded_beta_ALS_ALS2)

run_all_betas <- function(beta_list) {
  stopifnot(length(beta_list) > 0, all(names(beta_list) %in% levels(als_fac)))
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_vax_set(beta_vec, als_label)
  }) |>
    dplyr::arrange(als_scen, vax, date)
}

# run counterfactuals
df_all <- run_all_betas(beta_list)

# combine baseline and counterfactuals
plot_df <- bind_rows(df_all, baseline) |>
  mutate(
    als_scen = factor(als_scen, levels = c("ALS-2 ==> ALS-2")),
    vax = factor(vax, levels = c("25%", "50%", "75%")),
    panel_title = paste0("vaccination ", vax)
  )

# peaks from counterfactuals
cf_peaks <- plot_df %>%
  filter(
    series == "counterfactuals",
    date >= start_date, date <= date_end,
    vax %in% keep_vax, als_scen %in% als_cols
  ) %>%
  group_by(vax, als_scen, panel_title) %>%
  summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups = "drop"
  )

# calculate excess cases between counterfactuals & observed
diff_totals <- plot_df %>%
  filter(
    date >= start_date, date <= date_end,
    vax %in% keep_vax, als_scen %in% als_cols,
    series %in% c("counterfactuals","observed")
  ) %>%
  select(date, vax, als_scen, panel_title, series, inc) %>%
  pivot_wider(names_from = series, values_from = inc) %>%
  mutate(diff_inc = counterfactuals - observed) %>%
  group_by(vax, als_scen, panel_title) %>%
  summarise(excess_cases = sum(diff_inc, na.rm = TRUE), .groups = "drop")

# combine and build label (only excess_cases replaced)
cf_ann <- cf_peaks %>%
  left_join(diff_totals, by = c("vax","als_scen","panel_title")) %>%
  mutate(
    ann_lab = paste0(
      "Peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "Peak cases: ", comma(round(peak_infections)), "\n",
      "Excess cases: ", comma(round(excess_cases))
    )
  ) %>%
  select(vax, als_scen, panel_title, peak_infections, peak_date, excess_cases, ann_lab)


# ---- desired facet order (No-ALS → ALS-2 → ALS-4 within each vax)
vax_order <- c("25%","50%","75%")
als_order <- c("ALS-2 ==> ALS-2")

panel_levels <- c(
  "vaccination 25%",
  "vaccination 50%",
  "vaccination 75%"
)

plot_df <- plot_df %>%
  mutate(
    date       = as.Date(date),
    vax        = factor(vax, levels = vax_order),
    als_scen   = factor(als_scen, levels = als_order),
    panel_title= paste0("vaccination ", vax)
  ) %>%
  mutate(panel_title = factor(panel_title, levels = panel_levels)) 

# build annotations after the factor reset so they align to panels
ann_pos <- plot_df %>%
  group_by(vax, als_scen, panel_title) %>%
  summarise(
    xpos = min(date, na.rm=TRUE) + 1 + (max(date, na.rm=TRUE) - min(date, na.rm=TRUE)) * 0.15,
    ypos = 8500, .groups = "drop"
  )
ann_df <- cf_ann %>%
  inner_join(ann_pos, by=c("vax","als_scen","panel_title")) %>%
  mutate(panel_title = factor(panel_title, levels = panel_levels)) %>%
  droplevels()

# panel letters in this enforced order
panel_letters <- LETTERS[seq_along(panel_levels)]
panel_labels  <- tibble(panel_title = factor(panel_levels, levels = panel_levels), panel_label = panel_letters)
plot_df <- left_join(plot_df, panel_labels, by="panel_title")
ann_df  <- left_join(ann_df,  panel_labels, by="panel_title")

p <- ggplot(data = plot_df) +
  geom_line(data = filter(plot_df, series == "counterfactuals"),aes(date, inc, color = "Counterfactual"), linewidth = 1.5) +
  geom_point(data = filter(plot_df, series == "observed"),aes(date, inc, color = "Baseline"), size = 1.5) +
  facet_wrap(~ panel_title, nrow = 3, ncol = 3, axes = "all", scales = "fixed")+
  geom_text(data = ann_df, aes(x = xpos, y = ypos, label = ann_lab),inherit.aes = FALSE, size = 3.5, lineheight = 0.95, hjust = 0, vjust = 1) +
  geom_text(data = distinct(plot_df, panel_title, panel_label),aes(x = min(plot_df$date, na.rm=TRUE), y = max(plot_df$inc,  na.rm=TRUE),
                                                                   label = panel_label),hjust = -0.2, vjust = 1.2, size = 8, fontface = "bold",inherit.aes = FALSE) +
  scale_x_date(limits = c(date0, date_end), expand = c(0, 0),date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("Baseline" = "red", "Counterfactual" = "blue")) +
  labs(title = "Early Strong Restrictions Before Relaxation",x = NULL, y = "Estimated infections", color = NULL) +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 14, hjust = 0.5, face = "plain"),
    strip.text   = element_text(size = 14, face = "plain", hjust = 0.5),
    legend.position = "top",
    legend.background= element_blank(),
    legend.text      = element_text(size = 14, face = "plain"),
    panel.spacing   = unit(0.7, "lines"),
    plot.background = element_blank()
    #strip.background = element_blank()
  )

print(p)
###################### three (3) panels for manuscript ####################

cf_df <- plot_df |>
  filter(series == "counterfactuals",
         date >= date0, date <= date_end,
         vax %in% c("25%","50%","75%"),
         als_scen %in% levels(als_scen)
  )

# single observed baseline line (one per date)
obs_df <- plot_df |>
  filter(series == "observed",date >= date0, date <= date_end) |>
  arrange(date) |>
  group_by(date) |>
  summarise(inc = first(inc), .groups = "drop") |>
  mutate(series = "Baseline")

global_top <- max(c(cf_df$inc, obs_df$inc), na.rm = TRUE)

# panel labels placed at a common y across facets
label_pos <- cf_df |>
  distinct(als_scen) |>
  arrange(als_scen) |>
  mutate(
    x_lab = min(cf_df$date, na.rm = TRUE) + 3,
    y_lab = global_top * 0.99,                
    panel_label = LETTERS[seq_len(n())]
  )

# ALS scenario)
pp <- ggplot(data = plot_df) +
  geom_line(data = cf_df,aes(date, inc, color = vax),linewidth = 1.2, alpha = 1) +
  geom_point(data = obs_df,aes(date, inc, color = series),size = 1.2, alpha = 1) +
  facet_wrap(~ als_scen, nrow = 1, ncol = 3, axes = "all", scales = "fixed") +
  geom_text(data = label_pos,aes(x = x_lab, y = y_lab, label = panel_label),hjust = 0, vjust = 1, size = 8, fontface = "bold",inherit.aes = FALSE) +  
  scale_x_date(limits = c(date0, date_end),expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c("Baseline" = "red", "25%" = "#9467bd", "50%" = "#2ca02c", "75%" = "#1f77b4"),
    breaks = c("Baseline","25%","50%","75%"),
    name = NULL
  ) +
  labs(title = "Early Strong Restrictions Before Relaxation", x = NULL, y = "Estimated infections") +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(hjust = 0.5, face = "plain"),
    #strip.text   = element_text(size = 14, face = "plain", hjust = 0.5),
    legend.position = "top",
    legend.background = element_blank(),
    legend.text      = element_text(size = 14),
    panel.spacing    = unit(0.7, "lines"),
    plot.background  = element_blank(),
    strip.background = element_blank()
  )

print(pp)

print(p)

png("../figures/switch_als4_als2.png",width = 3200, height = 1500, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

png("../figures/switch_als4_als2_1.png",width = 2500, height = 1500, res = 300, bg = "white", type = "cairo")
print(pp)
dev.off()

