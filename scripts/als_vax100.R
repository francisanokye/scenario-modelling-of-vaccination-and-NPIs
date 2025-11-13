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
date0      <- as.Date("2021-12-15")
date_end   <- as.Date("2022-05-22")
start_date <- date0
dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
offset0 <- 150

# inputs 
calibrator    <- rdsRead("calibrate")
multipliers   <- rdsRead("build_multipliers")
seroprevdata  <- rdsRead("seroprevdata")
time_steps    <- max(seroprevdata$time)

# time-varying betas for scenarios
padded_base_beta    <- rdsRead("base_beta") 
padded_beta0_hat    <- rdsRead("policy_neutral") 
padded_beta_cf_ALS2 <- rdsRead("beta_als2")      
padded_beta_cf_ALS4 <- rdsRead("beta_als4")      

# baseline spec 
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# cache pristine vaccination series once (defensive)
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

# scenarios 
als_levels <- c("No-ALS Only","ALS-2 Only","ALS-4 Only")
als_fac    <- factor(als_levels, levels = als_levels)

# helpers (spec-aware)
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

# observed (replicated across panels)
baseline <- simulate_series(spec_fixed) |>
  mutate(series = "observed", inc = as.integer(round(inc))) |>
  tidyr::crossing(tibble(als_scen = als_fac)) |>
  select(date, inc, als_scen, series) |>
  arrange(als_scen, date)

############################# 100% vaccintaion ######################
# vax100_spec <- spec_fixed
# 
# # compute total current vaccinated (2-dose and booster)
# N <- vax100_spec$default$N
# sum_double <- sum(vax100_spec$default$double_vac_values_base, na.rm = TRUE)
# sum_booster <- sum(vax100_spec$default$booster_vac_values_base, na.rm = TRUE)
# 
# # scale factors so that total vaccination equals total population
# s2 <- N / sum_double
# s3 <- N / sum_booster
# 
# # apply scaled daily rates
# vax100_spec$default$double_vac_values  <- s2 * vax100_spec$default$double_vac_values_base
# vax100_spec$default$booster_vac_values <- s3 * vax100_spec$default$booster_vac_values_base

###################################################

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

df_all <- run_all_betas(spec_fixed, beta_list)

# combined plotting frame
plot_df <- bind_rows(df_all, baseline) |>
  mutate(
    als_scen    = factor(als_scen, levels = als_levels),
    panel_title = factor(paste0(als_scen), levels = als_levels)
  )

# summaries
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

# Figure 2: single-panel overlay 
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
  "ALS-4 Only"   = "#E69F00"      
)
legend_breaks <- c("Baseline", als_levels)

pp <- ggplot() +
  geom_line(data = cf_df,aes(x = date, y = inc, color = als_scen, group = als_scen),linewidth = 0.5) +
  geom_line(data = obs_df,aes(x = date, y = inc, color = "Baseline"),linewidth = 0.5) +
  scale_x_date(limits = c(date0, date_end), expand = c(0, 0),date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = pal, breaks = legend_breaks, name = NULL) +
  labs(title = "100% Vaccination and ALS counterfactuals vs baseline",x = NULL, y = "Estimated infections") +
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 1)),
    axis.text.y  = element_text(size = 6),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    plot.title   = element_text(size = 6, hjust = 0.5, face = "plain", margin = margin(b = 2)),
    legend.position = c(0.75, 0.75),
    legend.background = element_blank(),
    legend.text      = element_text(size = 6),
    legend.spacing.y = unit(0.8, 'cm'),
    legend.key.height = unit(0.5, "lines"),  
    legend.key.width  = unit(0.6, "lines"),   
    legend.margin     = margin(t = -5),       
    panel.spacing     = unit(0.3, "lines"),   
    plot.margin       = margin(t = 2, r = 2, b = 2, l = 2),  
    plot.background   = element_blank()
  )

png("../figures/als_vax_100_1.png",width = 1300, height = 800, res = 300, bg = "white", type = "cairo")
print(pp)
dev.off()
