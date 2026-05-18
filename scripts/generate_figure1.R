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
#   library(patchwork)
#   library(lubridate)
#   library(forcats)
# })
# 
# set.seed(2025)
# options(macpan2_log_dir = ".")
# options(dplyr.summarise.inform = FALSE)
# options(scipen = 9999)
# 
# loadEnvironments()
# 
# init_date <- as.Date("2021-07-24")
# analysis_start <- as.Date("2021-12-15")
# analysis_end <- as.Date("2022-05-22")
# release_date <- as.Date("2022-03-17")
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
# padded_beta_noals <- rdsRead("policy_neutral")
# padded_beta_base <- rdsRead("base_beta")
# padded_beta_als4 <- rdsRead("beta_als4")
# 
# stopifnot(exists("booster_inputs_df"))
# 
# generate_scenario_boosters <- function(scenario_name) {
#   scenario_daily <- booster_inputs_df |>
#     filter(scenario == scenario_name) |>
#     arrange(day) |>
#     transmute(
#       changepoints = as.integer(day),
#       vac_values = as.numeric(booster_vac_value)
#     )
#   
#   if (nrow(scenario_daily) == 0) {
#     stop("No booster_inputs_df rows found for scenario: ", scenario_name)
#   }
#   
#   scenario_daily <- scenario_daily |>
#     group_by(changepoints) |>
#     summarise(vac_values = sum(vac_values, na.rm = TRUE), .groups = "drop") |>
#     arrange(changepoints)
#   
#   if (!any(scenario_daily$changepoints == 0)) {
#     scenario_daily <- bind_rows(
#       tibble(changepoints = 0L, vac_values = 0),
#       scenario_daily
#     ) |>
#       arrange(changepoints)
#   }
#   
#   list(
#     changepoints = scenario_daily$changepoints,
#     vac_values = scenario_daily$vac_values
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
# spec_early <- generate_scenario_boosters("Early Boosters (-35 Days)") |>
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
#   "Early Boosters (-35 Days)" = spec_early,
#   "National Rate" = spec_natl
# )
# 
# policy_list <- list(
#   "Actual ALS" = padded_beta_base,
#   "ALS-4 Only" = padded_beta_als4,
#   "No-ALS Only" = padded_beta_noals
# )
# 
# policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")
# 
# vax_order <- c(
#   "National rate",
#   "Slower rate",
#   "Actual boosters",
#   "Early boosters (-35 Days)",
#   "Faster rate"
# )
# 
# dose_fill_cols <- c(
#   "Booster dose" = "#0072B2",
#   "\u22642 dose" = "goldenrod"
# )
# 
# policy_line_cols <- c(
#   "No NPIs" = "#009E73",
#   "Strict NPIs" = "#0072B2"
# )
# 
# speed_line_cols <- c(
#   "Slower boosters" = "purple",
#   "Faster boosters" = "blue"
# )
# 
# timing_line_cols <- c(
#   "Faster boosters" = "blue",
#   "Early boosters" = "red"
# )
# 
# sim_grid <- expand.grid(
#   policy_name = names(policy_list),
#   vax_name = names(vax_scenario_specs),
#   stringsAsFactors = FALSE
# )
# 
# all_sim_results <- vector("list", nrow(sim_grid))
# 
# for (i in seq_len(nrow(sim_grid))) {
#   policy_name <- sim_grid$policy_name[i]
#   vax_name <- sim_grid$vax_name[i]
#   
#   spec_i <- update_vax_spec_beta(
#     vax_scenario_specs[[vax_name]],
#     policy_list[[policy_name]]
#   )
#   
#   all_sim_results[[i]] <- spec_i |>
#     mp_simulator(
#       sim_horizon,
#       outputs = c(
#         "cum_booster",
#         "inc",
#         "inc_leq2dose",
#         "inc_booster",
#         "cum_inc_total",
#         "cum_inc_leq2dose",
#         "cum_inc_booster"
#       )
#     ) |>
#     mp_trajectory() |>
#     mutate(
#       date = init_date + time,
#       policy = policy_name,
#       vax_scenario = vax_name
#     )
# }
# 
# results_df <- bind_rows(all_sim_results) |>
#   mutate(
#     value = as.numeric(value),
#     policy = factor(policy, levels = policy_order),
#     vax_scenario = recode(
#       vax_scenario,
#       "National Rate" = "National rate",
#       "Slower Rate" = "Slower rate",
#       "Actual Boosters" = "Actual boosters",
#       "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
#       "Faster Rate" = "Faster rate"
#     ),
#     vax_scenario = factor(vax_scenario, levels = vax_order)
#   )
# 
# incidence_df <- results_df |>
#   filter(date >= analysis_start)
# 
# cumulative_df <- results_df |>
#   filter(date >= analysis_start)
# 
# ab_ymax <- incidence_df |>
#   filter(
#     vax_scenario == "Actual boosters",
#     policy %in% c("No-ALS Only", "ALS-4 Only"),
#     matrix %in% c("inc_leq2dose", "inc_booster")
#   ) |>
#   group_by(policy, date) |>
#   summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
#   summarise(ymax = max(total_inc, na.rm = TRUE)) |>
#   pull(ymax)
# 
# make_text_key <- function(
#     x,
#     y1,
#     y2,
#     label1 = "Solid: \u22642 dose",
#     label2 = "Dashed: Booster dose"
# ) {
#   list(
#     annotate(
#       "text",
#       x = x,
#       y = y1,
#       label = label1,
#       hjust = 0,
#       size = 9,
#       colour = "black"
#     ),
#     annotate(
#       "text",
#       x = x,
#       y = y2,
#       label = label2,
#       hjust = 0,
#       size = 9,
#       colour = "black"
#     )
#   )
# }
# 
# make_shaded_incidence_plot <- function(
#     df,
#     scenario_name,
#     policy_name,
#     panel_tag,
#     panel_title,
#     legend_position = c(0.03, 0.96),
#     y_limit = NULL
# ) {
#   plot_df <- df |>
#     filter(
#       policy == policy_name,
#       vax_scenario == scenario_name,
#       matrix %in% c("inc_leq2dose", "inc_booster")
#     ) |>
#     select(date, matrix, value) |>
#     pivot_wider(
#       names_from = matrix,
#       values_from = value,
#       values_fill = 0
#     ) |>
#     mutate(total_inc = inc_leq2dose + inc_booster)
#   
#   y_upper <- if (is.null(y_limit)) {
#     max(plot_df$total_inc, na.rm = TRUE)
#   } else {
#     y_limit
#   }
#   
#   ggplot(plot_df, aes(x = date)) +
#     geom_area(aes(y = total_inc, fill = "\u22642 dose"), alpha = 1) +
#     geom_area(aes(y = inc_booster, fill = "Booster dose"), alpha = 1) +
#     geom_line(aes(y = total_inc), colour = "black", linewidth = 2) +
#     geom_line(aes(y = inc_booster), colour = "black", linewidth = 2) +
#     geom_vline(
#       xintercept = release_date,
#       linetype = "dashed",
#       linewidth = 2,
#       colour = "black"
#     ) +
#     scale_fill_manual(
#       values = dose_fill_cols,
#       breaks = c("Booster dose", "\u22642 dose")
#     ) +
#     scale_x_date(
#       date_breaks = "1 month",
#       date_labels = "%d %b %Y",
#       expand = c(0.01, 0)
#     ) +
#     scale_y_continuous(
#       limits = c(0, y_upper * 1.06),
#       labels = comma,
#       expand = expansion(mult = c(0, 0.04))
#     ) +
#     labs(
#       title = panel_title,
#       x = NULL,
#       y = "Daily incidence",
#       fill = NULL,
#       tag = panel_tag
#     ) +
#     theme_clean() +
#     theme(
#       axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
#       axis.text.y = element_text(size = 30),
#       axis.title.y = element_text(size = 30, colour = "black"),
#       legend.position = legend_position,
#       legend.justification = "center",
#       legend.direction = "vertical",
#       legend.title = element_blank(),
#       legend.text = element_text(size = 25),
#       legend.background = element_blank(),
#       plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
#       panel.grid.major.x = element_blank(),
#       plot.tag = element_text(size = 35, face = "bold"),
#       plot.background = element_blank()
#     )
# }
# 
# make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
#   plot_df <- df |>
#     filter(
#       vax_scenario == scenario_name,
#       policy %in% c("No-ALS Only", "ALS-4 Only"),
#       matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
#     ) |>
#     mutate(
#       dose_group = recode(
#         matrix,
#         "cum_inc_leq2dose" = "\u22642 dose",
#         "cum_inc_booster" = "Booster dose"
#       ),
#       policy_label = recode(
#         as.character(policy),
#         "No-ALS Only" = "No NPIs",
#         "ALS-4 Only" = "Strict NPIs"
#       ),
#       policy_label = factor(policy_label, levels = c("No NPIs", "Strict NPIs")),
#       dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
#     )
#   
#   ymax <- max(plot_df$value, na.rm = TRUE)
#   xmax <- max(plot_df$date, na.rm = TRUE)
#   
#   label_df <- tibble(
#     policy_label = factor(c("No NPIs", "Strict NPIs"), levels = c("No NPIs", "Strict NPIs")),
#     label = c("No NPIs", "Strict NPIs"),
#     label_x = xmax - 75,
#     label_y = c(ymax * 0.20, ymax * 0.12)
#   )
#   
#   ggplot(
#     plot_df,
#     aes(
#       x = date,
#       y = value,
#       colour = policy_label,
#       linetype = dose_group,
#       group = interaction(policy_label, dose_group)
#     )
#   ) +
#     geom_line(linewidth = 3) +
#     geom_text(
#       data = label_df,
#       aes(x = label_x, y = label_y, label = label, colour = policy_label),
#       inherit.aes = FALSE,
#       size = 9,
#       fontface = "bold",
#       hjust = 0
#     ) +
#     make_text_key(
#       x = as.Date("2022-06-02"),
#       y1 = ymax * 0.18,
#       y2 = ymax * 0.10
#     ) +
#     scale_colour_manual(values = policy_line_cols) +
#     scale_linetype_manual(
#       values = c("\u22642 dose" = "solid", "Booster dose" = "twodash")
#     ) +
#     guides(colour = "none", linetype = "none") +
#     scale_x_date(
#       date_breaks = "1 month",
#       date_labels = "%d %b %Y",
#       expand = c(0.01, 0)
#     ) +
#     scale_y_continuous(
#       labels = comma,
#       expand = expansion(mult = c(0, 0.12))
#     ) +
#     labs(
#       title = NULL,
#       x = NULL,
#       y = "Cumulative infections",
#       colour = NULL,
#       linetype = NULL,
#       tag = panel_tag
#     ) +
#     theme_clean() +
#     theme(
#       axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
#       axis.text.y = element_text(size = 30),
#       axis.title.y = element_text(size = 30, colour = "black"),
#       legend.position = "none",
#       plot.title = element_blank(),
#       panel.grid.major.x = element_blank(),
#       plot.tag = element_text(size = 35, face = "bold"),
#       plot.background = element_blank()
#     )
# }
# 
# make_cumulative_speed_plot <- function(df, panel_tag) {
#   plot_df <- df |>
#     filter(
#       policy == "ALS-4 Only",
#       vax_scenario %in% c("Slower rate", "Faster rate"),
#       matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
#     ) |>
#     mutate(
#       dose_group = recode(
#         matrix,
#         "cum_inc_leq2dose" = "\u22642 dose",
#         "cum_inc_booster" = "Booster dose"
#       ),
#       rollout_speed = recode(
#         as.character(vax_scenario),
#         "Slower rate" = "Slower boosters",
#         "Faster rate" = "Faster boosters"
#       ),
#       rollout_speed = factor(rollout_speed, levels = c("Slower boosters", "Faster boosters")),
#       dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
#     )
#   
#   ymax <- max(plot_df$value, na.rm = TRUE)
#   xmax <- max(plot_df$date, na.rm = TRUE)
#   
#   label_df <- tibble(
#     rollout_speed = factor(c("Slower boosters", "Faster boosters"),
#                            levels = c("Slower boosters", "Faster boosters")),
#     label = c("Slower boosters", "Faster boosters"),
#     label_x = xmax - 90,
#     label_y = c(ymax * 0.20, ymax * 0.12)
#   )
#   
#   ggplot(
#     plot_df,
#     aes(
#       x = date,
#       y = value,
#       colour = rollout_speed,
#       linetype = dose_group,
#       group = interaction(rollout_speed, dose_group)
#     )
#   ) +
#     geom_line(linewidth = 3) +
#     geom_text(
#       data = label_df,
#       aes(x = label_x, y = label_y, label = label, colour = rollout_speed),
#       inherit.aes = FALSE,
#       size = 9,
#       fontface = "bold",
#       hjust = 0
#     ) +
#     make_text_key(
#       x = as.Date("2021-12-10"),
#       y1 = ymax * 0.92,
#       y2 = ymax * 0.84
#     ) +
#     scale_colour_manual(values = speed_line_cols) +
#     scale_linetype_manual(
#       values = c("\u22642 dose" = "solid", "Booster dose" = "twodash")
#     ) +
#     guides(colour = "none", linetype = "none") +
#     scale_x_date(
#       date_breaks = "1 month",
#       date_labels = "%d %b %Y",
#       expand = c(0.01, 0)
#     ) +
#     scale_y_continuous(
#       labels = comma,
#       expand = expansion(mult = c(0, 0.12))
#     ) +
#     labs(
#       title = NULL,
#       x = NULL,
#       y = "Cumulative infections",
#       colour = NULL,
#       linetype = NULL,
#       tag = panel_tag
#     ) +
#     theme_clean() +
#     theme(
#       axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
#       axis.text.y = element_text(size = 30),
#       axis.title.y = element_text(size = 30, colour = "black"),
#       legend.position = "none",
#       plot.title = element_blank(),
#       panel.grid.major.x = element_blank(),
#       plot.tag = element_text(size = 35, face = "bold"),
#       plot.background = element_blank()
#     )
# }
# 
# p_a <- make_shaded_incidence_plot(
#   incidence_df, "Actual boosters", "No-ALS Only",
#   "(A)", "No NPIs", c(0.50, 0.85), y_limit = ab_ymax
# )
# 
# p_b <- make_shaded_incidence_plot(
#   incidence_df, "Actual boosters", "ALS-4 Only",
#   "(B)", "Strict NPIs", c(0.50, 0.85), y_limit = ab_ymax
# )
# 
# p_c <- make_cumulative_policy_plot(
#   cumulative_df,
#   "Actual boosters",
#   "(C)"#, "Strict NPIs help expand booster coverage prior to infection"
# )
# 
# p_d <- make_shaded_incidence_plot(
#   incidence_df, "Slower rate", "ALS-4 Only",
#   "(D)", "Slower booster doses", c(0.50, 0.85)
# )
# 
# p_e <- make_shaded_incidence_plot(
#   incidence_df, "Faster rate", "ALS-4 Only",
#   "(E)", "Faster booster doses", c(0.12, 0.85)
# )
# 
# p_f <- make_cumulative_speed_plot(
#   cumulative_df,
#   "(F)"#, "Strict NPIs yield most benefit with rapid booster doses"
# )
# 
# p_g <- make_shaded_incidence_plot(
#   incidence_df, "Early boosters (-35 Days)", "ALS-4 Only",
#   "(G)", "Early boosters", c(0.12, 0.85)
# )
# 
# plot_h_df <- incidence_df |>
#   filter(
#     policy == "ALS-4 Only",
#     vax_scenario %in% c("Faster rate", "Early boosters (-35 Days)"),
#     matrix %in% c("inc", "inc_booster")
#   ) |>
#   mutate(
#     incidence_group = recode(
#       matrix,
#       "inc" = "Total incidence",
#       "inc_booster" = "Booster dose"
#     ),
#     scenario_label = recode(
#       as.character(vax_scenario),
#       "Faster rate" = "Faster boosters",
#       "Early boosters (-35 Days)" = "Early boosters"
#     ),
#     scenario_label = factor(scenario_label, levels = c("Faster boosters", "Early boosters")),
#     incidence_group = factor(incidence_group, levels = c("Total incidence", "Booster dose"))
#   )
# 
# h_label_df <- plot_h_df |>
#   filter(matrix == "inc") |>
#   group_by(scenario_label) |>
#   slice_max(value, n = 1, with_ties = FALSE) |>
#   ungroup() |>
#   mutate(
#     label_x = date + 25,
#     label_y = value * 1.06
#   )
# 
# h_ymax <- max(plot_h_df$value, na.rm = TRUE)
# 
# p_h <- ggplot(
#   plot_h_df,
#   aes(
#     x = date,
#     y = value,
#     colour = scenario_label,
#     linetype = incidence_group,
#     group = interaction(scenario_label, incidence_group)
#   )
# ) +
#   geom_line(linewidth = 3) +
#   geom_text(
#     data = h_label_df,
#     aes(x = label_x, y = label_y, label = scenario_label, colour = scenario_label),
#     inherit.aes = FALSE,
#     size = 9,
#     fontface = "bold",
#     hjust = 0
#   ) +
#   make_text_key(
#     x = as.Date("2021-12-10"),
#     y1 = h_ymax * 0.94,
#     y2 = h_ymax * 0.84,
#     label1 = "Solid: Total incidence",
#     label2 = "Dashed: Booster dose"
#   ) +
#   geom_vline(
#     xintercept = release_date,
#     linetype = "dashed",
#     linewidth = 2,
#     colour = "black"
#   ) +
#   scale_colour_manual(values = timing_line_cols) +
#   scale_linetype_manual(
#     values = c("Total incidence" = "solid", "Booster dose" = "twodash")
#   ) +
#   guides(colour = "none", linetype = "none") +
#   scale_x_date(
#     date_breaks = "1 month",
#     date_labels = "%d %b %Y",
#     expand = c(0.01, 0)
#   ) +
#   scale_y_continuous(
#     labels = comma,
#     expand = expansion(mult = c(0, 0.14))
#   ) +
#   labs(
#     title = "Early boosters increase total but not booster peaks",
#     x = NULL,
#     y = "Daily incidence",
#     colour = NULL,
#     linetype = NULL,
#     tag = "(H)"
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
#     axis.text.y = element_text(size = 30),
#     axis.title.y = element_text(size = 30, colour = "black"),
#     legend.position = "none",
#     plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
#     panel.grid.major.x = element_blank(),
#     plot.tag = element_text(size = 35, face = "bold"),
#     plot.background = element_blank()
#   )
# 
# plot_i_df <- incidence_df |>
#   filter(
#     policy == "ALS-4 Only",
#     vax_scenario %in% c("Faster rate", "Early boosters (-35 Days)"),
#     matrix %in% c("inc_leq2dose", "inc_booster")
#   ) |>
#   group_by(vax_scenario, matrix) |>
#   summarise(
#     peak_daily_incidence = max(value, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   mutate(
#     dose_group = recode(
#       matrix,
#       "inc_leq2dose" = "\u22642 dose",
#       "inc_booster" = "Booster dose"
#     ),
#     scenario_label = recode(
#       as.character(vax_scenario),
#       "Faster rate" = "Faster boosters",
#       "Early boosters (-35 Days)" = "Early boosters"
#     ),
#     scenario_label = factor(scenario_label, levels = c("Early boosters", "Faster boosters")),
#     dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
#   )
# 
# p_i <- ggplot(
#   plot_i_df,
#   aes(
#     x = scenario_label,
#     y = peak_daily_incidence,
#     fill = dose_group
#   )
# ) +
#   geom_col(
#     position = position_dodge(width = 0.72),
#     width = 0.68,
#     colour = "black",
#     linewidth = 0.25
#   ) +
#   scale_fill_manual(
#     values = dose_fill_cols,
#     breaks = c("Booster dose", "\u22642 dose")
#   ) +
#   scale_y_continuous(
#     labels = comma,
#     expand = expansion(mult = c(0, 0.1))
#   ) +
#   labs(
#     title = "Peak daily incidence by vaccination group",
#     x = NULL,
#     y = "Peak daily incidence",
#     fill = NULL,
#     tag = "(I)"
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 30, hjust = 0.5, angle = 0),
#     axis.text.y = element_text(size = 30),
#     axis.title.y = element_text(size = 30, colour = "black"),
#     legend.position = c(0.60, 0.85),
#     legend.justification = "center",
#     legend.direction = "vertical",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 25),
#     legend.background = element_blank(),
#     plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
#     panel.grid.major.x = element_blank(),
#     plot.tag = element_text(size = 35, face = "bold"),
#     plot.background = element_blank()
#   )
# 
# figure1_3x3 <- (
#   (p_a | p_b | p_c) /
#     (p_d | p_e | p_f) /
#     (p_g | p_h | p_i)
# ) +
#   plot_layout(
#     heights = c(1, 1, 1),
#     widths = c(1, 1, 1)
#   ) +
#   plot_annotation(
#     theme = theme(
#       plot.title = element_text(size = 50, hjust = 0.5, face = "plain")
#     )
#   )
# 
# save_png(
#   file.path(fig_dir, "figure1.png"),
#   figure1_3x3,
#   width = 13000,
#   height = 7500
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
  library(patchwork)
  library(lubridate)
  library(forcats)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

init_date <- as.Date("2021-07-24")
analysis_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")
release_date <- as.Date("2022-03-17")

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sim_horizon <- 530
new_V2_init <- 277707

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

calibrator <- rdsRead("calibrate")
padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base <- rdsRead("base_beta")
padded_beta_als4 <- rdsRead("beta_als4")

stopifnot(exists("booster_inputs_df"))

generate_scenario_boosters <- function(scenario_name) {
  scenario_daily <- booster_inputs_df |>
    filter(scenario == scenario_name) |>
    arrange(day) |>
    transmute(
      changepoints = as.integer(day),
      vac_values = as.numeric(booster_vac_value)
    )
  
  if (nrow(scenario_daily) == 0) {
    stop("No booster_inputs_df rows found for scenario: ", scenario_name)
  }
  
  scenario_daily <- scenario_daily |>
    group_by(changepoints) |>
    summarise(vac_values = sum(vac_values, na.rm = TRUE), .groups = "drop") |>
    arrange(changepoints)
  
  if (!any(scenario_daily$changepoints == 0)) {
    scenario_daily <- bind_rows(
      tibble(changepoints = 0L, vac_values = 0),
      scenario_daily
    ) |>
      arrange(changepoints)
  }
  
  list(
    changepoints = scenario_daily$changepoints,
    vac_values = scenario_daily$vac_values
  )
}

set_booster_vax <- function(booster_days, booster_vals) {
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
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

spec_actual <- generate_scenario_boosters("Actual Boosters") |>
  with(set_booster_vax(changepoints, vac_values))

spec_early <- generate_scenario_boosters("Early Boosters (-35 Days)") |>
  with(set_booster_vax(changepoints, vac_values))

spec_faster <- generate_scenario_boosters("Faster Rate") |>
  with(set_booster_vax(changepoints, vac_values))

spec_slower <- generate_scenario_boosters("Slower Rate") |>
  with(set_booster_vax(changepoints, vac_values))

spec_natl <- generate_scenario_boosters("National Rate") |>
  with(set_booster_vax(changepoints, vac_values)) |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))

vax_scenario_specs <- list(
  "Actual Boosters" = spec_actual,
  "Slower Rate" = spec_slower,
  "Faster Rate" = spec_faster,
  "Early Boosters (-35 Days)" = spec_early,
  "National Rate" = spec_natl
)

policy_list <- list(
  "Actual ALS" = padded_beta_base,
  "ALS-4 Only" = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "National rate",
  "Slower rate",
  "Actual boosters",
  "Early boosters (-35 Days)",
  "Faster rate"
)

dose_fill_cols <- c(
  "Booster dose" = "maroon",
  "\u22642 dose" = "gold"
)

policy_line_cols <- c(
  "No NPIs" = "red",
  "Strict NPIs" = "green4"
)

speed_line_cols <- c(
  "Slower boosters" = "magenta",
  "Faster boosters" = "#0026ff"
)

timing_line_cols <- c(
  "Faster boosters" = "#0026ff",
  "Early boosters" = "#D55E00"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  policy_name <- sim_grid$policy_name[i]
  vax_name <- sim_grid$vax_name[i]
  
  spec_i <- update_vax_spec_beta(
    vax_scenario_specs[[vax_name]],
    policy_list[[policy_name]]
  )
  
  all_sim_results[[i]] <- spec_i |>
    mp_simulator(
      sim_horizon,
      outputs = c(
        "cum_booster",
        "inc",
        "inc_leq2dose",
        "inc_booster",
        "cum_inc_total",
        "cum_inc_leq2dose",
        "cum_inc_booster"
      )
    ) |>
    mp_trajectory() |>
    mutate(
      date = init_date + time,
      policy = policy_name,
      vax_scenario = vax_name
    )
}

results_df <- bind_rows(all_sim_results) |>
  mutate(
    value = as.numeric(value),
    policy = factor(policy, levels = policy_order),
    vax_scenario = recode(
      vax_scenario,
      "National Rate" = "National rate",
      "Slower Rate" = "Slower rate",
      "Actual Boosters" = "Actual boosters",
      "Early Boosters (-35 Days)" = "Early boosters (-35 Days)",
      "Faster Rate" = "Faster rate"
    ),
    vax_scenario = factor(vax_scenario, levels = vax_order)
  )

incidence_df <- results_df |>
  filter(date >= analysis_start)

cumulative_df <- results_df |>
  filter(date >= analysis_start)

ab_ymax <- incidence_df |>
  filter(
    vax_scenario == "Actual boosters",
    policy %in% c("No-ALS Only", "ALS-4 Only"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(policy, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max(total_inc, na.rm = TRUE)) |>
  pull(ymax)

make_text_key <- function(
    x,
    y1,
    y2,
    label1 = "Solid: \u22642 dose",
    label2 = "Dashed: Booster dose"
) {
  list(
    annotate(
      "text",
      x = x,
      y = y1,
      label = label1,
      hjust = 0,
      size = 9,
      colour = "black"
    ),
    annotate(
      "text",
      x = x,
      y = y2,
      label = label2,
      hjust = 0,
      size = 9,
      colour = "black"
    )
  )
}

make_shaded_incidence_plot <- function(
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    legend_position = c(0.03, 0.96),
    y_limit = NULL
) {
  plot_df <- df |>
    filter(
      policy == policy_name,
      vax_scenario == scenario_name,
      matrix %in% c("inc_leq2dose", "inc_booster")
    ) |>
    select(date, matrix, value) |>
    pivot_wider(
      names_from = matrix,
      values_from = value,
      values_fill = 0
    ) |>
    mutate(total_inc = inc_leq2dose + inc_booster)
  
  y_upper <- if (is.null(y_limit)) {
    max(plot_df$total_inc, na.rm = TRUE)
  } else {
    y_limit
  }
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc, fill = "\u22642 dose"), alpha = 1) +
    geom_area(aes(y = inc_booster, fill = "Booster dose"), alpha = 1) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 2) +
    geom_line(aes(y = inc_booster), colour = "black", linewidth = 2) +
    geom_vline(
      xintercept = release_date,
      linetype = "dashed",
      linewidth = 2,
      colour = "black"
    ) +
    scale_fill_manual(
      values = dose_fill_cols,
      breaks = c("Booster dose", "\u22642 dose")
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_upper * 1.06),
      labels = comma,
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = "Daily incidence",
      fill = NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 30),
      axis.title.y = element_text(size = 30, colour = "black"),
      legend.position = legend_position,
      legend.justification = "center",
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 25),
      legend.background = element_blank(),
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 35, face = "bold"),
      plot.background = element_blank()
    )
}

make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
  plot_df <- df |>
    filter(
      vax_scenario == scenario_name,
      policy %in% c("No-ALS Only", "ALS-4 Only"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "\u22642 dose",
        "cum_inc_booster" = "Booster dose"
      ),
      policy_label = recode(
        as.character(policy),
        "No-ALS Only" = "No NPIs",
        "ALS-4 Only" = "Strict NPIs"
      ),
      policy_label = factor(policy_label, levels = c("No NPIs", "Strict NPIs")),
      dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
    )
  
  ymax <- max(plot_df$value, na.rm = TRUE)
  xmax <- max(plot_df$date, na.rm = TRUE)
  
  label_df <- tibble(
    policy_label = factor(c("No NPIs", "Strict NPIs"), levels = c("No NPIs", "Strict NPIs")),
    label = c("No NPIs", "Strict NPIs"),
    label_x = xmax - 75,
    label_y = c(ymax * 0.20, ymax * 0.12)
  )
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = policy_label,
      linetype = dose_group,
      group = interaction(policy_label, dose_group)
    )
  ) +
    geom_line(linewidth = 3) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = label_y, label = label, colour = policy_label),
      inherit.aes = FALSE,
      size = 9,
      fontface = "bold",
      hjust = 0
    ) +
    make_text_key(
      x = as.Date("2022-06-02"),
      y1 = ymax * 0.18,
      y2 = ymax * 0.10
    ) +
    scale_colour_manual(values = policy_line_cols) +
    scale_linetype_manual(
      values = c("\u22642 dose" = "solid", "Booster dose" = "twodash")
    ) +
    guides(colour = "none", linetype = "none") +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = "Strict NPIs increase booster coverage before infection",
      x = NULL,
      y = "Cumulative infections",
      colour = NULL,
      linetype = NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 30),
      axis.title.y = element_text(size = 30, colour = "black"),
      legend.position = "none",
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 35, face = "bold"),
      plot.background = element_blank()
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  plot_df <- df |>
    filter(
      policy == "ALS-4 Only",
      vax_scenario %in% c("Slower rate", "Faster rate"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "\u22642 dose",
        "cum_inc_booster" = "Booster dose"
      ),
      rollout_speed = recode(
        as.character(vax_scenario),
        "Slower rate" = "Slower boosters",
        "Faster rate" = "Faster boosters"
      ),
      rollout_speed = factor(rollout_speed, levels = c("Slower boosters", "Faster boosters")),
      dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
    )
  
  ymax <- max(plot_df$value, na.rm = TRUE)
  xmax <- max(plot_df$date, na.rm = TRUE)
  
  label_df <- tibble(
    rollout_speed = factor(c("Slower boosters", "Faster boosters"),
                           levels = c("Slower boosters", "Faster boosters")),
    label = c("Slower boosters", "Faster boosters"),
    label_x = xmax - 90,
    label_y = c(ymax * 0.20, ymax * 0.12)
  )
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = rollout_speed,
      linetype = dose_group,
      group = interaction(rollout_speed, dose_group)
    )
  ) +
    geom_line(linewidth = 3) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = label_y, label = label, colour = rollout_speed),
      inherit.aes = FALSE,
      size = 9,
      fontface = "bold",
      hjust = 0
    ) +
    make_text_key(
      x = as.Date("2021-12-10"),
      y1 = ymax * 0.92,
      y2 = ymax * 0.84
    ) +
    scale_colour_manual(values = speed_line_cols) +
    scale_linetype_manual(
      values = c("\u22642 dose" = "solid", "Booster dose" = "twodash")
    ) +
    guides(colour = "none", linetype = "none") +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = "Strict NPIs yield greatest benefit with rapid boosters",
      x = NULL,
      y = "Cumulative infections",
      colour = NULL,
      linetype = NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 30),
      axis.title.y = element_text(size = 30, colour = "black"),
      legend.position = "none",
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 35, face = "bold"),
      plot.background = element_blank()
    )
}

p_a <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "No-ALS Only",
  "(A)", "No NPIs", c(0.50, 0.85), y_limit = ab_ymax
)

p_b <- make_shaded_incidence_plot(
  incidence_df, "Actual boosters", "ALS-4 Only",
  "(B)", "Strict NPIs", c(0.50, 0.85), y_limit = ab_ymax
)

p_c <- make_cumulative_policy_plot(
  cumulative_df,
  "Actual boosters",
  "(C)"
)

p_d <- make_shaded_incidence_plot(
  incidence_df, "Slower rate", "ALS-4 Only",
  "(D)", "Slower booster doses", c(0.50, 0.85)
)

p_e <- make_shaded_incidence_plot(
  incidence_df, "Faster rate", "ALS-4 Only",
  "(E)", "Faster booster doses", c(0.12, 0.85)
)

p_f <- make_cumulative_speed_plot(
  cumulative_df,
  "(F)"
)

p_g <- make_shaded_incidence_plot(
  incidence_df, "Early boosters (-35 Days)", "ALS-4 Only",
  "(G)", "Early boosters", c(0.12, 0.85)
)

plot_h_df <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Faster rate", "Early boosters (-35 Days)"),
    matrix %in% c("inc", "inc_booster")
  ) |>
  mutate(
    incidence_group = recode(
      matrix,
      "inc" = "Total incidence",
      "inc_booster" = "Booster dose"
    ),
    scenario_label = recode(
      as.character(vax_scenario),
      "Faster rate" = "Faster boosters",
      "Early boosters (-35 Days)" = "Early boosters"
    ),
    scenario_label = factor(scenario_label, levels = c("Faster boosters", "Early boosters")),
    incidence_group = factor(incidence_group, levels = c("Total incidence", "Booster dose"))
  )

h_label_df <- plot_h_df |>
  filter(matrix == "inc") |>
  group_by(scenario_label) |>
  slice_max(value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    label_x = date + 25,
    label_y = value * 1.06
  )

h_ymax <- max(plot_h_df$value, na.rm = TRUE)

p_h <- ggplot(
  plot_h_df,
  aes(
    x = date,
    y = value,
    colour = scenario_label,
    linetype = incidence_group,
    group = interaction(scenario_label, incidence_group)
  )
) +
  geom_line(linewidth = 3) +
  geom_text(
    data = h_label_df,
    aes(x = label_x, y = label_y, label = scenario_label, colour = scenario_label),
    inherit.aes = FALSE,
    size = 9,
    fontface = "bold",
    hjust = 0
  ) +
  make_text_key(
    x = as.Date("2021-12-10"),
    y1 = h_ymax * 0.94,
    y2 = h_ymax * 0.84,
    label1 = "Solid: Total incidence",
    label2 = "Dashed: Booster dose"
  ) +
  geom_vline(
    xintercept = release_date,
    linetype = "dashed",
    linewidth = 2,
    colour = "black"
  ) +
  scale_colour_manual(values = timing_line_cols) +
  scale_linetype_manual(
    values = c("Total incidence" = "solid", "Booster dose" = "twodash")
  ) +
  guides(colour = "none", linetype = "none") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%d %b %Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = "Early boosters increase total but not booster peaks",
    x = NULL,
    y = "Daily incidence",
    colour = NULL,
    linetype = NULL,
    tag = "(H)"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.background = element_blank()
  )

plot_i_df <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Faster rate", "Early boosters (-35 Days)"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(
    peak_daily_incidence = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    dose_group = recode(
      matrix,
      "inc_leq2dose" = "\u22642 dose",
      "inc_booster" = "Booster dose"
    ),
    scenario_label = recode(
      as.character(vax_scenario),
      "Faster rate" = "Faster boosters",
      "Early boosters (-35 Days)" = "Early boosters"
    ),
    scenario_label = factor(scenario_label, levels = c("Early boosters", "Faster boosters")),
    dose_group = factor(dose_group, levels = c("\u22642 dose", "Booster dose"))
  )

p_i <- ggplot(
  plot_i_df,
  aes(
    x = scenario_label,
    y = peak_daily_incidence,
    fill = dose_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = dose_fill_cols,
    breaks = c("Booster dose", "\u22642 dose")
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Peak daily incidence by vaccination group",
    x = NULL,
    y = "Peak daily incidence",
    fill = NULL,
    tag = "(I)"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, colour = "black"),
    legend.position = c(0.60, 0.85),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.background = element_blank()
  )

figure1_3x3 <- (
  (p_a | p_b | p_c) /
    (p_d | p_e | p_f) /
    (p_g | p_h | p_i)
) +
  plot_layout(
    heights = c(1, 1, 1),
    widths = c(1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 50, hjust = 0.5, face = "plain")
    )
  )

save_png(
  file.path(fig_dir, "figure1.png"),
  figure1_3x3,
  width = 13000,
  height = 7500
)