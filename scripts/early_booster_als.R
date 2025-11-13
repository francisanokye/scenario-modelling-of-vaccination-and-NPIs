suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
})
set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ------------------------- config -------------------------
date0      <- as.Date("2021-12-15")
date_end   <- as.Date("2022-05-22")
start_date <- date0
dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
offset0 <- 150

# ------------------------- inputs -------------------------
calibrator    <- rdsRead("calibrate")
multipliers   <- rdsRead("build_multipliers")
seroprevdata  <- rdsRead("seroprevdata")
time_steps    <- max(seroprevdata$time)

# time-varying betas for scenarios
padded_base_beta    <- rdsRead("base_beta") # Observed ALS
padded_beta0_hat    <- rdsRead("policy_neutral") # No-ALS Only
padded_beta_cf_ALS2 <- rdsRead("beta_als2")      # ALS-2 Only
padded_beta_cf_ALS4 <- rdsRead("beta_als4")      # ALS-4 Only

# ------------------------- baseline spec -------------------------
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# cache pristine vaccination series once (defensive)
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

# ------------------------- scenarios -------------------------
als_levels <- c("Observed ALS", "No-ALS Only","ALS-2 Only","ALS-4 Only")
als_fac    <- factor(als_levels, levels = als_levels)

# ---- helpers (spec-aware) ----
simulate_series <- function(spec) {
  spec %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    filter(matrix == "inc") %>%
    mutate(date = date0 + (time - offset0)) %>%
    filter(date >= date0, date <= date_end) %>%
    transmute(date, inc = value)
}

update_beta <- function(spec, new_beta){
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_delete(phase = "before", at = 28L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

# ------------------------- observed & counterfactuals -------------------------
# observed (replicated across panels)
baseline <- simulate_series(spec_fixed) |>
  mutate(series = "observed",
         inc = as.integer(round(inc))) |>
  tidyr::crossing(tibble(als_scen = als_fac)) |>
  select(date, inc, als_scen, series) |>
  arrange(als_scen, date)

# scenario-specific counterfactuals (use observed vaccination; only beta differs)
simulate_set <- function(base_spec, new_beta, als_label) {
  stopifnot(als_label %in% levels(als_fac))
  spec_cf <- update_beta(base_spec, new_beta)
  simulate_series(spec_cf) |>
    mutate(
      date     = as.Date(date),
      inc      = as.integer(round(inc)),
      als_scen = factor(als_label, levels = levels(als_fac)),
      series   = "counterfactuals"
    ) |>
    select(date, inc, als_scen, series)
}

beta_list <- list(
  "Observed ALS" = padded_base_beta,
  "ALS-2 Only"  = padded_beta_cf_ALS2,
  "ALS-4 Only"  = padded_beta_cf_ALS4,
  "No-ALS Only" = padded_beta0_hat
)

run_all_betas <- function(base_spec, beta_list) {
  stopifnot(length(beta_list) > 0, all(names(beta_list) %in% levels(als_fac)))
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    simulate_set(base_spec, beta_vec, als_label)
  }) |> arrange(als_scen, date)
}

######################################### Early Booster Shots -28 days  ######################################### 

# clone calibrated spec for the counterfactual
cf_early28 <- spec_fixed

# shift booster timing
cp  <- as.integer(cf_early28$integers$booster_vac_changepoints_base) - 28L   
val <- as.numeric(cf_early28$default$booster_vac_values_base)

# clamp negatives to 0 (cannot start before t=0)
cp <- pmax(0L, cp)

# clean up changepoints and values and sort in ascending order
order_index <- order(cp)
cp  <- cp[order_index]
val <- val[order_index]

# write back into the scenario spec
cf_early28$integers$booster_vac_changepoints <- cp
cf_early28$default$booster_vac_values        <- val

df_all <- run_all_betas(cf_early28, beta_list)

# combined plotting frame
plot_df <- bind_rows(df_all, baseline) |>
  mutate(
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = factor(paste0(als_scen), levels = als_levels)
  )

# ------------------------- summaries -------------------------
# peaks from counterfactuals
cf_peaks <- plot_df %>%
  filter(series == "counterfactuals",
         date >= start_date, date <= date_end) %>%
  group_by(als_scen, panel_title) %>%
  summarise(
    peak_infections = max(inc, na.rm = TRUE),
    .idx = if (all(is.na(inc))) NA_integer_ else which.max(replace_na(inc, -Inf)),
    peak_date = if (is.na(.idx)) as.Date(NA) else date[.idx],
    .groups = "drop"
  )

# excess cases (counterfactuals vs observed), join to avoid list-cols
cf_only <- plot_df %>%
  filter(series == "counterfactuals",
         date >= start_date, date <= date_end) %>%
  select(date, als_scen, panel_title, inc_cf = inc)

obs_only <- plot_df %>%
  filter(series == "observed",
         date >= start_date, date <= date_end) %>%
  select(date, als_scen, panel_title, inc_obs = inc)

diff_totals <- cf_only %>%
  left_join(obs_only, by = c("date","als_scen","panel_title")) %>%
  mutate(diff_inc = inc_cf - inc_obs) %>%
  group_by(als_scen, panel_title) %>%
  summarise(excess_cases = sum(diff_inc, na.rm = TRUE), .groups = "drop")

# annotation labels & positions
cf_ann <- cf_peaks %>%
  left_join(diff_totals, by = c("als_scen","panel_title")) %>%
  mutate(
    ann_lab = paste0(
      "Peak date: ", format(peak_date, "%b %d, %Y"), "\n",
      "Peak cases: ", comma(round(peak_infections)), "\n",
      "Excess cases: ", comma(round(excess_cases))
    )
  ) %>%
  select(als_scen, panel_title, peak_infections, peak_date, excess_cases, ann_lab)

ann_pos <- plot_df %>%
  group_by(als_scen, panel_title) %>%
  summarise(
    xpos = min(date, na.rm=TRUE) + 1 + (max(date, na.rm=TRUE) - min(date, na.rm=TRUE)) * 0.1,
    ypos = 4000,
    .groups = "drop"
  )
ann_df <- cf_ann %>%
  inner_join(ann_pos, by=c("als_scen","panel_title")) %>%
  droplevels()

panel_letters <- LETTERS[seq_along(als_levels)]
panel_labels  <- tibble(panel_title = factor(als_levels, levels = als_levels),
                        panel_label = panel_letters)
plot_df <- left_join(plot_df, panel_labels, by="panel_title")
ann_df  <- left_join(ann_df,  panel_labels, by="panel_title")

# ------------------------- Figure 1: three facets -------------------------
p <- ggplot(data = plot_df) +
  geom_line(data = filter(plot_df, series == "counterfactuals"),aes(date, inc, color = "Counterfactual"), linewidth = 1.2) +
  geom_line(data = filter(plot_df, series == "observed"),aes(date, inc, color = "Baseline"), linewidth = 1.2) +
  facet_wrap(~ panel_title, nrow = 2, ncol = 2, scales = "free_y") +
  geom_text(data = ann_df,aes(x = xpos, y = ypos, label = ann_lab),inherit.aes = FALSE, size = 4.5, lineheight = 0.95, hjust = 0, vjust = 1) +
  geom_text(data = distinct(plot_df, panel_title, panel_label),aes(x = min(plot_df$date, na.rm=TRUE),y = max(plot_df$inc,  na.rm=TRUE),label = panel_label),
            hjust = -0.2, vjust = -1.2, size = 5, fontface = "bold", inherit.aes = FALSE) +
  scale_x_date(limits = c(date0, date_end), expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("Baseline" = "red", "Counterfactual" = "blue"), name = NULL) +
  labs(title = "Early Boosters and ALS counterfactuals (blue) vs Observed Baseline (red)",x = NULL, y = "Estimated infections") +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 14, hjust = 0.5),
    strip.text   = element_text(size = 14),
    legend.position = "top",
    legend.background= element_blank(),
    legend.text      = element_text(size = 14),
    panel.spacing   = unit(0.7, "lines"),
    plot.background = element_blank()
  )
print(p)

# ------------------------- Figure 2: single-panel overlay -------------------------
cf_df <- plot_df %>%
  filter(series == "counterfactuals",
         date >= date0, date <= date_end) %>%
  arrange(als_scen, date)

obs_df <- plot_df %>%
  filter(series == "observed",
         date >= date0, date <= date_end) %>%
  arrange(date) %>%
  group_by(date) %>%
  summarise(inc = first(inc), .groups = "drop") %>%
  mutate(series = "Baseline")

pal <- c(
  "Baseline"     = "red",          
  "No-ALS Only"  = "#0072B2",      
  "ALS-2 Only"   = "#009E73",      
  "ALS-4 Only"   = "#E69F00",      
  "Observed ALS" = "#CC79A7"       
)
legend_breaks <- c("Baseline", als_levels)

pp <- ggplot() +
  geom_line(data = cf_df,aes(x = date, y = inc, color = als_scen, group = als_scen),linewidth = 0.5) +
  geom_line(data = obs_df,aes(x = date, y = inc, color = "Baseline"),linewidth = 0.5) +
  scale_x_date(limits = c(date0, date_end), expand = c(0, 0),date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = pal, breaks = legend_breaks, name = NULL) +
  labs(title = "Early Boosters and ALS counterfactuals (blue) vs Observed Baseline (red)",x = NULL, y = "Estimated infections") +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 1)),
    axis.text.y  = element_text(size = 6),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    plot.title   = element_text(size = 6, hjust = 0.5, face = "plain", margin = margin(b = 2)),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.text      = element_text(size = 6),
    legend.key.height = unit(0.2, "lines"),   # reduce legend vertical padding
    legend.key.width  = unit(0.6, "lines"),   # shrink horizontal spacing between keys
    legend.margin     = margin(t = -5),       # pull legend closer to plot
    panel.spacing     = unit(0.3, "lines"),   # tighten space between panels
    plot.margin       = margin(t = 2, r = 2, b = 2, l = 2),  # reduce outer margins
    plot.background   = element_blank()
  )

print(pp)

# ------------------------- save -------------------------
# png("../figures/early_booster_als.png",  width = 4000, height = 2000, res = 300, bg = "white", type = "cairo"); print(p);  dev.off()
# png("../figures/early_booster_als_1.png",width = 1500, height = 800, res = 300, bg = "white", type = "cairo"); print(pp); dev.off()

png("../figures/early_booster_als.png",width = 4000, height = 2000, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

png("../figures/early_booster_als_1.png",width = 1300, height = 800, res = 300, bg = "white", type = "cairo")
print(pp)
dev.off()
