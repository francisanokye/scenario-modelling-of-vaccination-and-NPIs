# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(ggthemes)
#   library(shellpipes)
#   library(macpan2)
#   library(scales)
#   library(grid)
#   library(lubridate)
# })
# 
# set.seed(2025)
# 
# options(macpan2_log_dir = ".")
# options(dplyr.summarise.inform = FALSE)
# options(scipen = 9999)
# 
# loadEnvironments()
# 
# init_date      <- as.Date("2021-07-24")
# analysis_start <- as.Date("2021-12-15")
# analysis_end   <- as.Date("2022-05-22")
# converge_date  <- as.Date("2022-03-14")
# 
# fig_dir <- "../figures"
# dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# 
# sim_horizon <- 530
# new_V2_init <- 277707
# 
# save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
#   png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
#   print(plot_obj)
#   dev.off()
#   invisible(path)
# }
# 
# calibrator <- rdsRead("calibrate")
# 
# padded_beta_noals <- rdsRead("policy_neutral")
# padded_beta_base  <- rdsRead("base_beta")
# padded_beta_als4  <- rdsRead("beta_als4")
# 
# stopifnot(
#   length(padded_beta_noals) == length(padded_beta_base),
#   length(padded_beta_als4)  == length(padded_beta_base),
#   exists("booster_inputs_df")
# )
# 
# beta_check_df <- tibble(
#   time = seq_along(padded_beta_base),
#   date = init_date + time,
#   beta_actual = as.numeric(padded_beta_base),
#   beta_noals  = as.numeric(padded_beta_noals),
#   beta_als4   = as.numeric(padded_beta_als4)
# ) |>
#   filter(date >= analysis_start, date <= analysis_end) |>
#   mutate(
#     period = if_else(date < converge_date, "Before Mar 14", "From Mar 14 onward"),
#     diff_noals_actual = beta_noals - beta_actual,
#     diff_als4_actual  = beta_als4 - beta_actual
#   )
# 
# beta_diagnostic <- beta_check_df |>
#   group_by(period) |>
#   summarise(
#     n_days = n(),
#     mean_actual = mean(beta_actual, na.rm = TRUE),
#     mean_noals  = mean(beta_noals, na.rm = TRUE),
#     mean_als4   = mean(beta_als4, na.rm = TRUE),
#     max_abs_diff_noals_actual = max(abs(diff_noals_actual), na.rm = TRUE),
#     max_abs_diff_als4_actual  = max(abs(diff_als4_actual), na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(beta_diagnostic)
# 
# if (
#   max(abs(beta_check_df$beta_als4 - beta_check_df$beta_actual), na.rm = TRUE) < 1e-10
# ) {
#   warning("Actual NPIs and Strict NPIs use identical beta values across the analysis window.")
# }
# 
# if (
#   max(abs(
#     beta_check_df$beta_als4[beta_check_df$date >= converge_date] -
#     beta_check_df$beta_actual[beta_check_df$date >= converge_date]
#   ), na.rm = TRUE) < 1e-10
# ) {
#   message("Strict NPIs match Actual NPIs from March 14 onward.")
# }
# 
# if (
#   max(abs(
#     beta_check_df$beta_noals[beta_check_df$date >= converge_date] -
#     beta_check_df$beta_actual[beta_check_df$date >= converge_date]
#   ), na.rm = TRUE) < 1e-10
# ) {
#   message("No NPIs match Actual NPIs from March 14 onward.")
# }
# 
# booster_start_dates <- c(
#   "Actual Boosters"           = as.Date("2021-11-08"),
#   "Early Boosters (-21 Days)" = as.Date("2021-10-04"),
#   "Faster Rate"               = as.Date("2021-11-08"),
#   "Slower Rate"               = as.Date("2021-11-08"),
#   "National Rate"             = as.Date("2021-11-08")
# )
# 
# generate_scenario_boosters <- function(scenario_name) {
#   
#   if (!scenario_name %in% names(booster_start_dates)) {
#     stop("No booster start date defined for scenario: ", scenario_name)
#   }
#   
#   start_day <- as.integer(
#     booster_start_dates[[scenario_name]] - init_date
#   )
#   
#   scenario_daily <- booster_inputs_df |>
#     filter(scenario == scenario_name) |>
#     transmute(
#       changepoints = as.integer(day),
#       vac_values   = as.numeric(booster_vac_value)
#     ) |>
#     filter(changepoints >= start_day) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
#   
#   if (nrow(scenario_daily) == 0) {
#     stop("No booster_inputs_df rows found for scenario: ", scenario_name)
#   }
#   
#   scenario_daily <- bind_rows(
#     tibble(changepoints = 0L, vac_values = 0),
#     scenario_daily
#   ) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
#   
#   active_rows <- scenario_daily |>
#     filter(vac_values > 0)
#   
#   if (nrow(active_rows) == 0) {
#     stop("All booster values are zero for scenario: ", scenario_name)
#   }
#   
#   print(
#     active_rows |>
#       summarise(
#         scenario = scenario_name,
#         first_booster_day = min(changepoints),
#         first_booster_date = init_date + first_booster_day,
#         last_booster_day = max(changepoints),
#         last_booster_date = init_date + last_booster_day,
#         total_boosters = sum(vac_values),
#         max_daily_boosters = max(vac_values)
#       )
#   )
#   
#   list(
#     changepoints = scenario_daily$changepoints,
#     vac_values   = scenario_daily$vac_values
#   )
# }
# 
# set_booster_vax <- function(booster_days, booster_vals) {
#   spec <- mp_optimized_spec(calibrator, "modified")
#   spec$integers$booster_vac_changepoints <- as.integer(booster_days)
#   spec$default$booster_vac_values <- as.numeric(booster_vals)
#   spec
# }
# 
# update_vax_spec_beta <- function(spec, new_beta) {
#   
#   stopifnot(length(new_beta) == spec$default$time_steps)
#   
#   spec |>
#     mp_tmb_delete(phase = "before", at = 26L) |>
#     mp_tmb_delete(phase = "before", at = 26L) |>
#     mp_tmb_insert(default = list(timeseries_beta = new_beta))
# }
# 
# spec_actual <- generate_scenario_boosters("Actual Boosters") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_early <- generate_scenario_boosters("Early Boosters (-21 Days)") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_faster <- generate_scenario_boosters("Faster Rate") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_slower <- generate_scenario_boosters("Slower Rate") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_natl <- generate_scenario_boosters("National Rate") |>
#   with(set_booster_vax(changepoints, vac_values)) |>
#   mp_tmb_insert(default = list(V2_init = new_V2_init))
# 
# vax_scenario_specs <- list(
#   "Slower Rate"               = spec_slower,
#   "Actual Boosters"           = spec_actual,
#   "Faster Rate"               = spec_faster,
#   "Early Boosters (-21 Days)" = spec_early,
#   "National Rate"             = spec_natl
# )
# 
# policy_list <- list(
#   "No-ALS Only" = padded_beta_noals,
#   "Actual ALS"  = padded_beta_base,
#   "ALS-4 Only"  = padded_beta_als4
# )
# 
# policy_order <- c(
#   "No-ALS Only",
#   "Actual ALS",
#   "ALS-4 Only"
# )
# 
# severity_fill_cols <- c(
#   "Slower rate"               = "#CC79A7",
#   "Actual boosters"           = "#4D4D4D",
#   "Faster rate"               = "#0072B2",
#   "Early Boosters (-21 Days)" = "#D55E00",
#   "National rate"             = "springgreen"
# )
# 
# vax_short_labs <- c(
#   "Slower rate"               = "Slower\nboosters",
#   "Actual boosters"           = "Actual\nboosters",
#   "Faster rate"               = "Faster\nboosters",
#   "Early Boosters (-21 Days)" = "Earlier\nboosters",
#   "National rate"             = "National\nrate"
# )
# 
# sim_grid <- expand.grid(
#   policy_name = names(policy_list),
#   vax_name    = names(vax_scenario_specs),
#   stringsAsFactors = FALSE
# )
# 
# all_sim_results <- vector("list", nrow(sim_grid))
# 
# for (i in seq_len(nrow(sim_grid))) {
#   
#   spec_i <- update_vax_spec_beta(
#     vax_scenario_specs[[sim_grid$vax_name[i]]],
#     policy_list[[sim_grid$policy_name[i]]]
#   )
#   
#   all_sim_results[[i]] <- spec_i |>
#     mp_simulator(
#       sim_horizon,
#       outputs = c("hosp_total", "death_total")
#     ) |>
#     mp_trajectory() |>
#     mutate(
#       date = init_date + time,
#       policy = sim_grid$policy_name[i],
#       vax_scenario = sim_grid$vax_name[i]
#     )
# }
# 
# severity_df <- bind_rows(all_sim_results) |>
#   mutate(
#     value = as.numeric(value),
#     policy = factor(policy, levels = policy_order),
#     vax_scenario = recode(
#       vax_scenario,
#       "Slower Rate"               = "Slower rate",
#       "Actual Boosters"           = "Actual boosters",
#       "Faster Rate"               = "Faster rate",
#       "Early Boosters (-21 Days)" = "Early Boosters (-21 Days)",
#       "National Rate"             = "National rate"
#     )
#   ) |>
#   pivot_wider(
#     names_from = matrix,
#     values_from = value,
#     values_fill = 0
#   ) |>
#   filter(
#     date >= analysis_start,
#     date <= analysis_end
#   ) |>
#   arrange(policy, vax_scenario, date) |>
#   group_by(policy, vax_scenario) |>
#   mutate(
#     cum_hosp_total  = cumsum(hosp_total),
#     cum_death_total = cumsum(death_total)
#   ) |>
#   ungroup()
# 
# diagnostic_summary <- severity_df |>
#   group_by(policy, vax_scenario) |>
#   summarise(
#     hospitalizations = max(cum_hosp_total, na.rm = TRUE),
#     deaths = max(cum_death_total, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(diagnostic_summary, n = Inf)
# 
# figure3_df <- bind_rows(
#   severity_df |>
#     group_by(policy, vax_scenario) |>
#     summarise(
#       severe_outcome = max(cum_hosp_total, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     mutate(outcome = "Hospitalizations"),
#   
#   severity_df |>
#     group_by(policy, vax_scenario) |>
#     summarise(
#       severe_outcome = max(cum_death_total, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     mutate(outcome = "Deaths")
# ) |>
#   mutate(
#     policy_label = recode(
#       as.character(policy),
#       "No-ALS Only" = "No NPIs",
#       "Actual ALS"  = "Actual NPIs",
#       "ALS-4 Only"  = "Strict NPIs"
#     ),
#     policy_label = factor(
#       policy_label,
#       levels = c("No NPIs", "Actual NPIs", "Strict NPIs")
#     ),
#     outcome = factor(
#       outcome,
#       levels = c("Hospitalizations", "Deaths")
#     )
#   )
# 
# x_order_df <- figure3_df |>
#   filter(outcome == "Hospitalizations") |>
#   group_by(policy_label) |>
#   arrange(desc(severe_outcome), .by_group = TRUE) |>
#   mutate(x_pos = row_number()) |>
#   select(policy_label, vax_scenario, x_pos)
# 
# figure3_df <- figure3_df |>
#   left_join(x_order_df, by = c("policy_label", "vax_scenario")) |>
#   mutate(
#     x_lab = vax_short_labs[vax_scenario],
#     value_label = comma(round(severe_outcome, 0))
#   )
# 
# figure3 <- ggplot(
#   figure3_df,
#   aes(
#     x = x_pos,
#     y = severe_outcome,
#     fill = vax_scenario
#   )
# ) +
#   geom_col(
#     colour = "black",
#     width = 0.72,
#     linewidth = 0.25
#   ) +
#   geom_text(
#     aes(label = value_label),
#     vjust = -0.35,
#     size = 10,
#     colour = "black"
#   ) +
#   facet_grid(
#     outcome ~ policy_label,
#     scales = "free_y"
#   ) +
#   scale_fill_manual(
#     values = severity_fill_cols,
#     breaks = names(severity_fill_cols),
#     labels = c(
#       "Slower rate" = "Slower boosters",
#       "Actual boosters" = "Actual boosters",
#       "Faster rate" = "Faster boosters",
#       "Early Boosters (-21 Days)" = "Earlier boosters",
#       "National rate" = "National rate"
#     )
#   ) +
#   scale_x_continuous(
#     breaks = 1:5,
#     labels = NULL,
#     expand = expansion(mult = c(0.04, 0.04))
#   ) +
#   scale_y_continuous(
#     labels = comma,
#     expand = expansion(mult = c(0.08, 0.18))
#   ) +
#   labs(
#     title = "Severe outcomes across NPI and booster scenarios",
#     x = NULL,
#     y = "Severe outcomes",
#     fill = NULL
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 35),
#     axis.ticks.x = element_blank(),
#     axis.text.y = element_text(size = 35),
#     axis.title.y = element_text(size = 35, colour = "black"),
#     strip.text = element_text(size = 35, face = "plain"),
#     strip.background = element_blank(),
#     legend.position = "top",
#     legend.text = element_text(size = 35),
#     legend.background = element_blank(),
#     legend.spacing.x = unit(0.45, "cm"),
#     plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
#     panel.grid.major.x = element_blank(),
#     panel.spacing = unit(2, "cm"),
#     plot.margin = margin(t = 10, r = 20, b = 55, l = 20),
#     plot.background = element_blank()
#   ) +
#   geom_text(
#     data = figure3_df |>
#       filter(outcome == "Deaths"),
#     aes(
#       y = -Inf,
#       label = x_lab
#     ),
#     vjust = 1.15,
#     size = 8,
#     lineheight = 0.85,
#     inherit.aes = TRUE
#   ) +
#   coord_cartesian(clip = "off")
# 
# save_png(
#   file.path(fig_dir, "figure2.png"),
#   figure3,
#   width = 10000,
#   height = 5000
# )





# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(ggthemes)
#   library(shellpipes)
#   library(macpan2)
#   library(scales)
#   library(grid)
#   library(lubridate)
# })
# 
# set.seed(2025)
# 
# options(macpan2_log_dir = ".")
# options(dplyr.summarise.inform = FALSE)
# options(scipen = 9999)
# 
# loadEnvironments()
# 
# analysis_start <- as.Date("2021-12-15")
# offset <- 150L
# init_date <- analysis_start - offset
# 
# analysis_end   <- as.Date("2022-12-15") #as.Date("2022-05-22")
# converge_date  <- as.Date("2022-03-14")
# 
# actual_booster_start <- as.Date("2021-11-08")
# early_shift_days <- -14
# early_start_date <- actual_booster_start + days(early_shift_days)
# 
# early_scenario_name <- paste0("Early Boosters (", early_shift_days, " Days)")
# early_plot_label <- paste0("Early boosters (", early_shift_days, " Days)")
# 
# stopifnot(
#   offset == as.integer(analysis_start - init_date),
#   analysis_start + (0 - offset) == init_date,
#   actual_booster_start >= init_date,
#   early_start_date >= init_date
# )
# 
# fig_dir <- "../figures"
# dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# 
# sim_horizon <- 500
# new_V2_init <- 277707
# 
# save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
#   png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
#   print(plot_obj)
#   dev.off()
#   invisible(path)
# }
# 
# calibrator <- rdsRead("calibrate")
# 
# padded_beta_noals <- rdsRead("policy_neutral")
# padded_beta_base  <- rdsRead("base_beta")
# padded_beta_als4  <- rdsRead("beta_als4")
# 
# stopifnot(
#   length(padded_beta_noals) == length(padded_beta_base),
#   length(padded_beta_als4) == length(padded_beta_base),
#   exists("booster_inputs_df")
# )
# 
# beta_check_df <- tibble(
#   time = seq_along(padded_beta_base) - 1L,
#   date = analysis_start + (time - offset),
#   beta_actual = as.numeric(padded_beta_base),
#   beta_noals  = as.numeric(padded_beta_noals),
#   beta_als4   = as.numeric(padded_beta_als4)
# ) |>
#   filter(date >= analysis_start, date <= analysis_end) |>
#   mutate(
#     period = if_else(date < converge_date, "Before Mar 14", "From Mar 14 onward"),
#     diff_noals_actual = beta_noals - beta_actual,
#     diff_als4_actual  = beta_als4 - beta_actual
#   )
# 
# print(
#   beta_check_df |>
#     group_by(period) |>
#     summarise(
#       n_days = n(),
#       mean_actual = mean(beta_actual, na.rm = TRUE),
#       mean_noals = mean(beta_noals, na.rm = TRUE),
#       mean_als4 = mean(beta_als4, na.rm = TRUE),
#       max_abs_diff_noals_actual = max(abs(diff_noals_actual), na.rm = TRUE),
#       max_abs_diff_als4_actual = max(abs(diff_als4_actual), na.rm = TRUE),
#       .groups = "drop"
#     )
# )
# 
# booster_start_dates <- c(
#   "Actual Boosters" = actual_booster_start,
#   setNames(early_start_date, early_scenario_name),
#   "Faster Rate" = actual_booster_start,
#   "Slower Rate" = actual_booster_start,
#   "National Rate" = actual_booster_start
# )
# 
# generate_scenario_boosters <- function(scenario_name) {
#   
#   scenario_match <- booster_inputs_df |>
#     distinct(scenario) |>
#     filter(tolower(scenario) == tolower(scenario_name)) |>
#     pull(scenario)
#   
#   if (length(scenario_match) != 1) {
#     stop(
#       "Scenario name not found: ", scenario_name,
#       "\nAvailable scenarios are:\n",
#       paste(unique(booster_inputs_df$scenario), collapse = "\n")
#     )
#   }
#   
#   scenario_name_actual <- scenario_match
#   
#   start_lookup <- tibble(
#     scenario = names(booster_start_dates),
#     start_date = as.Date(booster_start_dates),
#     scenario_clean = tolower(names(booster_start_dates))
#   ) |>
#     filter(scenario_clean == tolower(scenario_name_actual))
#   
#   if (nrow(start_lookup) != 1) {
#     stop(
#       "No booster start date defined for scenario: ", scenario_name_actual,
#       "\nAvailable start-date scenarios are:\n",
#       paste(names(booster_start_dates), collapse = "\n")
#     )
#   }
#   
#   start_day <- as.integer(start_lookup$start_date - init_date)
#   
#   scenario_daily <- booster_inputs_df |>
#     filter(scenario == scenario_name_actual) |>
#     transmute(
#       changepoints = as.integer(day),
#       vac_values = as.numeric(booster_vac_value)
#     ) |>
#     filter(changepoints >= start_day) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
#   
#   if (nrow(scenario_daily) == 0) {
#     stop("No booster input rows found for scenario: ", scenario_name_actual)
#   }
#   
#   scenario_daily <- bind_rows(
#     tibble(changepoints = 0L, vac_values = 0),
#     scenario_daily
#   ) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
#   
#   active_rows <- scenario_daily |>
#     filter(vac_values > 0)
#   
#   if (nrow(active_rows) == 0) {
#     stop("All booster values are zero for scenario: ", scenario_name_actual)
#   }
#   
#   print(
#     active_rows |>
#       summarise(
#         scenario = scenario_name_actual,
#         intended_start_date = start_lookup$start_date,
#         first_booster_day = min(changepoints),
#         first_booster_date = init_date + first_booster_day,
#         last_booster_day = max(changepoints),
#         last_booster_date = init_date + last_booster_day,
#         total_boosters = sum(vac_values),
#         max_daily_boosters = max(vac_values)
#       )
#   )
#   
#   list(
#     changepoints = scenario_daily$changepoints,
#     vac_values = scenario_daily$vac_values
#   )
# }
# 
# set_booster_vax <- function(booster_days, booster_vals) {
#   stopifnot(length(booster_days) == length(booster_vals))
#   stopifnot(all(diff(booster_days) >= 0))
#   
#   spec <- mp_optimized_spec(calibrator, "modified")
#   spec$integers$booster_vac_changepoints <- as.integer(booster_days)
#   spec$default$booster_vac_values <- as.numeric(booster_vals)
#   spec
# }
# 
# update_vax_spec_beta <- function(spec, new_beta) {
#   stopifnot(length(new_beta) == spec$default$time_steps)
#   
#   spec |>
#     mp_tmb_delete(phase = "before", at = 26L) |>
#     mp_tmb_delete(phase = "before", at = 26L) |>
#     mp_tmb_insert(default = list(timeseries_beta = new_beta))
# }
# 
# spec_actual <- generate_scenario_boosters("Actual Boosters") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_early <- generate_scenario_boosters(early_scenario_name) |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_faster <- generate_scenario_boosters("Faster Rate") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_slower <- generate_scenario_boosters("Slower Rate") |>
#   with(set_booster_vax(changepoints, vac_values))
# 
# spec_natl <- generate_scenario_boosters("National Rate") |>
#   with(set_booster_vax(changepoints, vac_values)) |>
#   mp_tmb_insert(default = list(V2_init = new_V2_init))
# 
# vax_scenario_specs <- list(
#   "Actual Boosters" = spec_actual,
#   "Slower Rate" = spec_slower,
#   "Faster Rate" = spec_faster,
#   "National Rate" = spec_natl
# )
# 
# vax_scenario_specs[[early_scenario_name]] <- spec_early
# 
# policy_list <- list(
#   "No-ALS Only" = padded_beta_noals,
#   "Actual ALS" = padded_beta_base,
#   "ALS-4 Only" = padded_beta_als4
# )
# 
# policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")
# 
# severity_fill_cols <- c(
#   "National rate" = "springgreen",
#   "Slower rate" = "#CC79A7",
#   "Actual boosters" = "#4D4D4D",
#   setNames("#D55E00", early_plot_label),
#   "Faster rate" = "#0072B2"
# )
# 
# vax_short_labs <- c(
#   "National rate" = "National\nrate",
#   "Slower rate" = "Slower\nboosters",
#   "Actual boosters" = "Actual\nboosters",
#   setNames("Earlier\nboosters", early_plot_label),
#   "Faster rate" = "Faster\nboosters"
# )
# 
# vax_order <- c(
#   "National rate",
#   "Slower rate",
#   "Actual boosters",
#   early_plot_label,
#   "Faster rate"
# )
# 
# sim_grid <- expand.grid(
#   policy_name = names(policy_list),
#   vax_name = names(vax_scenario_specs),
#   stringsAsFactors = FALSE
# )
# 
# model_outputs <- c(
#   "incS",
#   "incv2",
#   "incv3",
#   "inc_leq2dose",
#   "inc_booster",
#   "cum_inc_leq2dose",
#   "cum_inc_booster",
#   "cum_inc_total",
#   "hosp_leq2dose",
#   "hosp_booster",
#   "hosp_total",
#   "death_leq2dose",
#   "death_booster",
#   "death_total",
#   "cum_hosp_leq2dose",
#   "cum_hosp_booster",
#   "cum_hosp_total",
#   "cum_death_leq2dose",
#   "cum_death_booster",
#   "cum_death_total"
# )
# 
# all_sim_results <- vector("list", nrow(sim_grid))
# 
# for (i in seq_len(nrow(sim_grid))) {
#   
#   spec_i <- update_vax_spec_beta(
#     vax_scenario_specs[[sim_grid$vax_name[i]]],
#     policy_list[[sim_grid$policy_name[i]]]
#   )
#   
#   all_sim_results[[i]] <- spec_i |>
#     mp_simulator(
#       time_steps = sim_horizon,
#       outputs = model_outputs
#     ) |>
#     mp_trajectory() |>
#     mutate(
#       date = analysis_start + (time - offset),
#       policy = sim_grid$policy_name[i],
#       vax_scenario = sim_grid$vax_name[i]
#     )
# }
# 
# 
# 
# 
# severity_df <- bind_rows(all_sim_results) |>
#   mutate(
#     value = as.numeric(value),
#     policy = factor(policy, levels = policy_order),
#     vax_scenario = recode(
#       vax_scenario,
#       "National Rate" = "National rate",
#       "Slower Rate" = "Slower rate",
#       "Actual Boosters" = "Actual boosters",
#       "Faster Rate" = "Faster rate",
#       .default = vax_scenario
#     ),
#     vax_scenario = if_else(
#       tolower(vax_scenario) == tolower(early_scenario_name),
#       early_plot_label,
#       vax_scenario
#     ),
#     vax_scenario = factor(vax_scenario, levels = vax_order)
#   ) |>
#   filter(date >= analysis_start, date <= analysis_end)
# 
# diagnostic_summary <- severity_df |>
#   group_by(policy, vax_scenario) |>
#   summarise(
#     inc_leq2dose = max(value[matrix == "cum_inc_leq2dose"], na.rm = TRUE),
#     inc_booster = max(value[matrix == "cum_inc_booster"], na.rm = TRUE),
#     hospitalizations = max(value[matrix == "cum_hosp_total"], na.rm = TRUE),
#     deaths = max(value[matrix == "cum_death_total"], na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(diagnostic_summary, n = Inf)
# 
# figure3_df <- bind_rows(
#   severity_df |>
#     filter(matrix == "cum_hosp_total") |>
#     group_by(policy, vax_scenario) |>
#     summarise(
#       severe_outcome = max(value, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     mutate(outcome = "Hospitalizations"),
#   
#   severity_df |>
#     filter(matrix == "cum_death_total") |>
#     group_by(policy, vax_scenario) |>
#     summarise(
#       severe_outcome = max(value, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     mutate(outcome = "Deaths")
# ) |>
#   mutate(
#     policy_label = recode(
#       as.character(policy),
#       "No-ALS Only" = "No NPIs",
#       "Actual ALS" = "Actual NPIs",
#       "ALS-4 Only" = "Strict NPIs"
#     ),
#     policy_label = factor(policy_label, levels = c("No NPIs", "Actual NPIs", "Strict NPIs")),
#     outcome = factor(outcome, levels = c("Hospitalizations", "Deaths")),
#     vax_scenario = factor(vax_scenario, levels = vax_order)
#   )
# 
# x_order_df <- figure3_df |>
#   filter(outcome == "Hospitalizations") |>
#   group_by(policy_label) |>
#   arrange(desc(severe_outcome), .by_group = TRUE) |>
#   mutate(x_pos = row_number()) |>
#   ungroup() |>
#   select(policy_label, vax_scenario, x_pos)
# 
# figure3_df <- figure3_df |>
#   left_join(x_order_df, by = c("policy_label", "vax_scenario")) |>
#   mutate(
#     x_lab = vax_short_labs[as.character(vax_scenario)],
#     value_label = comma(round(severe_outcome, 0))
#   )
# 
# figure3 <- ggplot(
#   figure3_df,
#   aes(
#     x = x_pos,
#     y = severe_outcome,
#     fill = vax_scenario
#   )
# ) +
#   geom_col(
#     colour = "black",
#     width = 0.72,
#     linewidth = 0.25
#   ) +
#   geom_text(
#     aes(label = value_label),
#     vjust = -0.35,
#     size = 10,
#     colour = "black"
#   ) +
#   facet_grid(
#     outcome ~ policy_label,
#     scales = "free_y"
#   ) +
#   scale_fill_manual(
#     values = severity_fill_cols,
#     breaks = vax_order,
#     labels = c(
#       "National rate" = "National rate",
#       "Slower rate" = "Slower boosters",
#       "Actual boosters" = "Actual boosters",
#       setNames("Earlier boosters", early_plot_label),
#       "Faster rate" = "Faster boosters"
#     ),
#     drop = FALSE
#   ) +
#   scale_x_continuous(
#     breaks = 1:5,
#     labels = NULL,
#     expand = expansion(mult = c(0.04, 0.04))
#   ) +
#   scale_y_continuous(
#     labels = comma,
#     expand = expansion(mult = c(0.08, 0.18))
#   ) +
#   labs(
#     title = "Severe outcomes across NPI and booster scenarios",
#     x = NULL,
#     y = "Severe outcomes",
#     fill = NULL
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 35),
#     axis.ticks.x = element_blank(),
#     axis.text.y = element_text(size = 35),
#     axis.title.y = element_text(size = 35, colour = "black"),
#     strip.text = element_text(size = 35, face = "plain"),
#     strip.background = element_blank(),
#     legend.position = "top",
#     legend.text = element_text(size = 35),
#     legend.background = element_blank(),
#     legend.spacing.x = unit(0.45, "cm"),
#     plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
#     panel.grid.major.x = element_blank(),
#     panel.spacing = unit(2, "cm"),
#     plot.margin = margin(t = 10, r = 20, b = 55, l = 20),
#     plot.background = element_blank()
#   ) +
#   geom_text(
#     data = figure3_df |> filter(outcome == "Deaths"),
#     aes(
#       y = -Inf,
#       label = x_lab
#     ),
#     vjust = 1.15,
#     size = 8,
#     lineheight = 0.85,
#     inherit.aes = TRUE
#   ) +
#   coord_cartesian(clip = "off")
# 
# save_png(
#   file.path(fig_dir, "figure2.png"),
#   figure3,
#   width = 10000,
#   height = 5000
# )




suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(lubridate)
})

set.seed(2025)

options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

analysis_start <- as.Date("2021-12-15")
offset <- 150L
init_date <- analysis_start - offset

sim_horizon <- 530L
analysis_end <- init_date + sim_horizon - 1L

actual_booster_start <- as.Date("2021-11-08")
new_V2_init <- 277707

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  init_date + offset == analysis_start,
  analysis_start + (0 - offset) == init_date
)

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

calibrator <- rdsRead("calibrate")

padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

stopifnot(
  length(padded_beta_noals) == length(padded_beta_base),
  length(padded_beta_als4) == length(padded_beta_base),
  exists("booster_inputs_df")
)

available_early <- booster_inputs_df |>
  distinct(scenario) |>
  filter(grepl("^Early Boosters", scenario, ignore.case = TRUE)) |>
  pull(scenario)

if (length(available_early) != 1) {
  stop(
    "Expected exactly one Early Boosters scenario, but found:\n",
    paste(available_early, collapse = "\n")
  )
}

early_scenario_name <- available_early

early_shift_days <- as.integer(
  gsub(".*\\((-?[0-9]+) Days\\).*", "\\1", early_scenario_name)
)

early_start_date <- actual_booster_start + days(early_shift_days)

early_plot_label <- paste0(
  "Early boosters (",
  early_shift_days,
  " Days)"
)

stopifnot(early_start_date >= init_date)

booster_start_dates <- c(
  "Actual Boosters" = actual_booster_start,
  setNames(early_start_date, early_scenario_name),
  "Faster Rate" = actual_booster_start,
  "Slower Rate" = actual_booster_start,
  "National Rate" = actual_booster_start
)

generate_scenario_boosters <- function(scenario_name) {
  
  scenario_match <- booster_inputs_df |>
    distinct(scenario) |>
    filter(tolower(scenario) == tolower(scenario_name)) |>
    pull(scenario)
  
  if (length(scenario_match) != 1) {
    stop(
      "Scenario name not found: ", scenario_name,
      "\nAvailable scenarios are:\n",
      paste(unique(booster_inputs_df$scenario), collapse = "\n")
    )
  }
  
  scenario_name_actual <- scenario_match
  
  start_lookup <- tibble(
    scenario = names(booster_start_dates),
    start_date = as.Date(booster_start_dates),
    scenario_clean = tolower(names(booster_start_dates))
  ) |>
    filter(scenario_clean == tolower(scenario_name_actual))
  
  if (nrow(start_lookup) != 1) {
    stop(
      "No booster start date defined for scenario: ",
      scenario_name_actual
    )
  }
  
  start_day <- as.integer(start_lookup$start_date - init_date)
  
  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name_actual) |>
    transmute(
      changepoints = as.integer(day),
      vac_values = as.numeric(booster_vac_value)
    ) |>
    filter(changepoints >= start_day) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster input rows found for scenario: ", scenario_name_actual)
  }
  
  scenario_daily <- bind_rows(
    tibble(changepoints = 0L, vac_values = 0),
    scenario_daily
  ) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
  
  active_rows <- scenario_daily |>
    filter(vac_values > 0)
  
  if (nrow(active_rows) == 0) {
    stop("All booster values are zero for scenario: ", scenario_name_actual)
  }
  
  print(
    active_rows |>
      summarise(
        scenario = scenario_name_actual,
        intended_start_date = start_lookup$start_date,
        first_booster_day = min(changepoints),
        first_booster_date = init_date + first_booster_day,
        last_booster_day = max(changepoints),
        last_booster_date = init_date + last_booster_day,
        total_boosters = sum(vac_values),
        max_daily_boosters = max(vac_values)
      )
  )
  
  list(
    changepoints = scenario_daily$changepoints,
    vac_values = scenario_daily$vac_values
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
  stopifnot(length(booster_days) == length(booster_vals))
  stopifnot(all(diff(booster_days) >= 0))
  
  spec <- mp_optimized_spec(calibrator, "modified")
  spec$integers$booster_vac_changepoints <- as.integer(booster_days)
  spec$default$booster_vac_values <- as.numeric(booster_vals)
  spec
}

update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  
  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(
      default = list(timeseries_beta = new_beta)
    )
}

spec_actual <- generate_scenario_boosters("Actual Boosters") |>
  with(set_booster_vax(changepoints, vac_values))

spec_early <- generate_scenario_boosters(early_scenario_name) |>
  with(set_booster_vax(changepoints, vac_values))

spec_faster <- generate_scenario_boosters("Faster Rate") |>
  with(set_booster_vax(changepoints, vac_values))

spec_slower <- generate_scenario_boosters("Slower Rate") |>
  with(set_booster_vax(changepoints, vac_values))

spec_natl <- generate_scenario_boosters("National Rate") |>
  with(set_booster_vax(changepoints, vac_values)) |>
  mp_tmb_insert(
    default = list(V2_init = new_V2_init)
  )

vax_scenario_specs <- list(
  "Actual Boosters" = spec_actual,
  "Slower Rate" = spec_slower,
  "Faster Rate" = spec_faster,
  "National Rate" = spec_natl
)

vax_scenario_specs[[early_scenario_name]] <- spec_early

policy_list <- list(
  "No-ALS Only" = padded_beta_noals,
  "Actual ALS" = padded_beta_base,
  "ALS-4 Only" = padded_beta_als4
)

policy_order <- c(
  "No-ALS Only",
  "Actual ALS",
  "ALS-4 Only"
)

vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  early_plot_label,
  "Faster rate"
)

severity_fill_cols <- c(
  "National rate" = "springgreen",
  "Slower rate" = "#CC79A7",
  "Actual boosters" = "#4D4D4D",
  setNames("#D55E00", early_plot_label),
  "Faster rate" = "#0072B2"
)

vax_short_labs <- c(
  "National rate" = "National\nrate",
  "Slower rate" = "Slower\nboosters",
  "Actual boosters" = "Actual\nboosters",
  setNames("Earlier\nboosters", early_plot_label),
  "Faster rate" = "Faster\nboosters"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

model_outputs <- c(
  "incS",
  "incv2",
  "incv3",
  "inc_leq2dose",
  "inc_booster",
  "cum_inc_leq2dose",
  "cum_inc_booster",
  "cum_inc_total",
  "hosp_leq2dose",
  "hosp_booster",
  "hosp_total",
  "death_leq2dose",
  "death_booster",
  "death_total",
  "cum_hosp_leq2dose",
  "cum_hosp_booster",
  "cum_hosp_total",
  "cum_death_leq2dose",
  "cum_death_booster",
  "cum_death_total"
)

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  
  spec_i <- update_vax_spec_beta(
    vax_scenario_specs[[sim_grid$vax_name[i]]],
    policy_list[[sim_grid$policy_name[i]]]
  )
  
  all_sim_results[[i]] <- spec_i |>
    mp_simulator(
      time_steps = sim_horizon,
      outputs = model_outputs
    ) |>
    mp_trajectory() |>
    mutate(
      date = analysis_start + (time - offset),
      policy = sim_grid$policy_name[i],
      vax_scenario = sim_grid$vax_name[i]
    )
}

severity_df <- bind_rows(all_sim_results) |>
  mutate(
    date = as.Date(date),
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      vax_scenario,
      "National Rate" = "National rate",
      "Slower Rate" = "Slower rate",
      "Actual Boosters" = "Actual boosters",
      "Faster Rate" = "Faster rate",
      .default = vax_scenario
    ),
    vax_scenario = if_else(
      tolower(vax_scenario) == tolower(early_scenario_name),
      early_plot_label,
      vax_scenario
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  ) |>
  filter(
    date >= analysis_start,
    date <= analysis_end
  )

diagnostic_summary <- severity_df |>
  group_by(policy, vax_scenario) |>
  summarise(
    inc_leq2dose = max(
      value[matrix == "cum_inc_leq2dose"],
      na.rm = TRUE
    ),
    inc_booster = max(
      value[matrix == "cum_inc_booster"],
      na.rm = TRUE
    ),
    hospitalizations = max(
      value[matrix == "cum_hosp_total"],
      na.rm = TRUE
    ),
    deaths = max(
      value[matrix == "cum_death_total"],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(diagnostic_summary, n = Inf)

figure3_df <- bind_rows(
  severity_df |>
    filter(matrix == "cum_hosp_total") |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Hospitalizations"),
  
  severity_df |>
    filter(matrix == "cum_death_total") |>
    group_by(policy, vax_scenario) |>
    summarise(
      severe_outcome = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outcome = "Deaths")
) |>
  mutate(
    policy_label = recode(
      as.character(policy),
      "No-ALS Only" = "No NPIs",
      "Actual ALS" = "Actual NPIs",
      "ALS-4 Only" = "Strict NPIs"
    ),
    policy_label = factor(
      policy_label,
      levels = c("No NPIs", "Actual NPIs", "Strict NPIs")
    ),
    outcome = factor(
      outcome,
      levels = c("Hospitalizations", "Deaths")
    ),
    vax_scenario = factor(
      vax_scenario,
      levels = vax_order
    )
  )

x_order_df <- figure3_df |>
  filter(outcome == "Hospitalizations") |>
  group_by(policy_label) |>
  arrange(desc(severe_outcome), .by_group = TRUE) |>
  mutate(x_pos = row_number()) |>
  ungroup() |>
  select(
    policy_label,
    vax_scenario,
    x_pos
  )

figure3_df <- figure3_df |>
  left_join(
    x_order_df,
    by = c("policy_label", "vax_scenario")
  ) |>
  mutate(
    x_lab = vax_short_labs[as.character(vax_scenario)],
    value_label = comma(round(severe_outcome, 0))
  )

figure3 <- ggplot(
  figure3_df,
  aes(
    x = x_pos,
    y = severe_outcome,
    fill = vax_scenario
  )
) +
  geom_col(
    colour = "black",
    width = 0.72,
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = value_label),
    vjust = -0.35,
    size = 10,
    colour = "black"
  ) +
  facet_grid(
    outcome ~ policy_label,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = severity_fill_cols,
    breaks = vax_order,
    labels = c(
      "National rate" = "National rate",
      "Slower rate" = "Slower boosters",
      "Actual boosters" = "Actual boosters",
      setNames("Earlier boosters", early_plot_label),
      "Faster rate" = "Faster boosters"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:5,
    labels = NULL,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Severe outcomes",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 35),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35, colour = "black"),
    strip.text = element_text(size = 35, face = "plain"),
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 35),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.45, "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm"),
    plot.margin = margin(t = 10, r = 20, b = 55, l = 20),
    plot.background = element_blank()
  ) +
  geom_text(
    data = figure3_df |> filter(outcome == "Deaths"),
    aes(
      y = -Inf,
      label = x_lab
    ),
    vjust = 1.15,
    size = 8,
    lineheight = 0.85,
    inherit.aes = TRUE
  ) +
  coord_cartesian(clip = "off")

save_png(
  file.path(fig_dir, "figure2.png"),
  figure3,
  width = 10000,
  height = 5000
)