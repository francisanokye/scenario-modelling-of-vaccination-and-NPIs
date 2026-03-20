# # ============================================================================#
# #  FAST policy-duration counterfactuals (peak daily infections threshold)
# #  Goal: For each ALS setting in isolation (No-ALS, ALS-2, ALS-4),
# #        find the MIN extra days beyond March 14, 2022 such that
# #        peak daily infections (Dec 17, 2021–May 22, 2022) <= 1261,
# #        given observed vaccination pace (100%).
# # ============================================================================#
# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(purrr)
#   library(tibble)
#   library(ggplot2)
#   library(ggthemes)
#   library(scales)
#   library(shellpipes)
#   library(macpan2)
# })
# 
# set.seed(2025)
# options(macpan2_log_dir = ".")
# options(dplyr.summarise.inform = FALSE)
# options(scipen = 9999)
# loadEnvironments()
# 
# # ------------------------------- Config ------------------------------------#
# omicron_start_date    <- as.Date("2021-12-17")
# study_end_date        <- as.Date("2022-05-22")
# offset_days           <- 150L
# 
# measures_lift_date    <- as.Date("2022-03-14")
# infection_peak_target <- 1261
# 
# vaccination_scale_observed <- 1.0
# 
# # Search controls (fast)
# extra_days_max     <- 120L            # maximum extension searched
# coarse_step_days   <- 7L              # coarse scan step to bracket solution
# bisection_max_iter <- 12L             # ~2^12 = 4096 resolution, more than enough for days
# 
# # Outputs / figures
# out_dir <- "../output"
# fig_dir <- "../figures/policy_duration_infection_target_fast"
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# 
# csv_sum <- file.path(out_dir, "min_policy_duration_infections1261_summary_fast.csv")
# csv_det <- file.path(out_dir, "min_policy_duration_infections1261_detail_fast.csv")
# 
# save_png <- function(path, p, width, height, dpi = 300) {
#   ggsave(filename = path, plot = p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
#   invisible(path)
# }
# 
# # ------------------------------- Inputs ------------------------------------#
# calibrator  <- rdsRead("calibrate")
# fitserodata <- rdsRead("fitsero.rds")
# stopifnot("time" %in% names(fitserodata))
# time_steps <- max(fitserodata$time, na.rm = TRUE)
# stopifnot(is.finite(time_steps), time_steps > 0)
# 
# beta_no_als <- rdsRead("policy_neutral")
# beta_als2   <- rdsRead("beta_als2")
# beta_als4   <- rdsRead("beta_als4")
# 
# check_beta <- function(x, nm) {
#   if (is.null(x) || length(x) == 0) stop("Missing beta series: ", nm)
#   if (length(x) != time_steps) stop("Beta series `", nm, "` length (", length(x), ") != time_steps (", time_steps, ").")
#   if (!all(is.finite(x))) stop("Beta series `", nm, "` contains non-finite values.")
#   invisible(TRUE)
# }
# check_beta(beta_no_als, "policy_neutral")
# check_beta(beta_als2,   "beta_als2")
# check_beta(beta_als4,   "beta_als4")
# 
# policy_beta_sets <- list(
#   "No-ALS Only" = beta_no_als,
#   "ALS-2 Only"  = beta_als2,
#   "ALS-4 Only"  = beta_als4
# )
# 
# # --------------------------- Base spec / caching ----------------------------#
# spec_base <- mp_optimized_spec(calibrator, "modified")
# spec_base$default$time_steps <- time_steps
# 
# # Cache pristine vaccination series once
# if (is.null(spec_base$default$double_vac_values_base))
#   spec_base$default$double_vac_values_base  <- spec_base$default$double_vac_values
# if (is.null(spec_base$default$booster_vac_values_base))
#   spec_base$default$booster_vac_values_base <- spec_base$default$booster_vac_values
# stopifnot(!is.null(spec_base$default$double_vac_values_base),
#           !is.null(spec_base$default$booster_vac_values_base))
# 
# set_vaccination_observed <- function(spec) {
#   spec$default$double_vac_values  <- vaccination_scale_observed * spec$default$double_vac_values_base
#   spec$default$booster_vac_values <- vaccination_scale_observed * spec$default$booster_vac_values_base
#   spec
# }
# 
# # Insert timeseries_beta ONCE; afterwards we only overwrite spec$default$timeseries_beta
# spec_with_beta_slot <- function(spec, beta_vec) {
#   stopifnot(length(beta_vec) == spec$default$time_steps)
#   spec2 <- spec
#   
#   # Add the beta slot once (delete/insert only once)
#   spec2 <- spec2 |>
#     mp_tmb_delete(phase = "before", at = 24L) |>
#     mp_tmb_delete(phase = "before", at = 24L) |>
#     mp_tmb_insert(default = list(timeseries_beta = beta_vec))
#   
#   spec2
# }
# 
# # Map time indices to dates (for finding lift index)
# time_date_map <- tibble(
#   time = seq_len(time_steps),
#   date = omicron_start_date + (time - offset_days)
# )
# 
# lift_idx <- which(time_date_map$date == measures_lift_date)
# if (length(lift_idx) != 1L) {
#   lift_idx <- which.min(abs(as.numeric(time_date_map$date - measures_lift_date)))
#   warning("measures_lift_date not uniquely found; using nearest mapped date: ",
#           as.character(time_date_map$date[lift_idx]), " (index ", lift_idx, ").")
# }
# 
# extend_policy_duration <- function(beta_vec, extra_days) {
#   beta_ext <- beta_vec
#   restricted_level <- beta_vec[max(1L, lift_idx - 1L)]
#   end_hold_idx <- min(length(beta_vec), lift_idx + as.integer(extra_days))
#   if (end_hold_idx >= lift_idx) beta_ext[lift_idx:end_hold_idx] <- restricted_level
#   beta_ext
# }
# 
# simulate_peak_daily_infections <- function(spec) {
#   df <- spec |>
#     mp_simulator(time_steps, outputs = c("inc")) |>
#     mp_trajectory() |>
#     filter(matrix == "inc") |>
#     transmute(
#       date = omicron_start_date + (as.integer(time) - offset_days),
#       daily_infections = pmax(0, as.numeric(value))
#     ) |>
#     filter(date >= omicron_start_date, date <= study_end_date)
#   
#   max(df$daily_infections, na.rm = TRUE)
# }
# 
# # -------------------------- Fast search (coarse + bisection) ----------------#
# # Returns: list(min_extra_days, detail_df)
# find_min_extra_days_fast <- function(spec_template, beta_base, als_label) {
#   
#   # Coarse scan to bracket the first feasible day
#   coarse_days <- unique(pmin(extra_days_max, seq(0L, extra_days_max, by = coarse_step_days)))
#   coarse_res <- map_dfr(coarse_days, function(d) {
#     beta_cf <- extend_policy_duration(beta_base, d)
#     spec_template$default$timeseries_beta <- beta_cf
#     peak <- simulate_peak_daily_infections(spec_template)
#     tibble(extra_policy_days = d, peak_daily_infections = peak)
#   }) |>
#     mutate(meets_target = peak_daily_infections <= infection_peak_target,
#            als_policy = als_label)
#   
#   if (any(coarse_res$meets_target)) {
#     hi <- min(coarse_res$extra_policy_days[coarse_res$meets_target])
#     lo <- max(0L, hi - coarse_step_days)
#     
#     # Bisection on integer days within [lo, hi]
#     iter <- 0L
#     while (lo < hi && iter < bisection_max_iter) {
#       mid <- as.integer(floor((lo + hi) / 2))
#       beta_cf <- extend_policy_duration(beta_base, mid)
#       spec_template$default$timeseries_beta <- beta_cf
#       peak <- simulate_peak_daily_infections(spec_template)
#       
#       coarse_res <- bind_rows(coarse_res, tibble(
#         extra_policy_days = mid,
#         peak_daily_infections = peak,
#         meets_target = peak <= infection_peak_target,
#         als_policy = als_label
#       ))
#       
#       if (peak <= infection_peak_target) hi <- mid else lo <- mid + 1L
#       iter <- iter + 1L
#     }
#     
#     min_day <- hi
#   } else {
#     min_day <- NA_integer_
#   }
#   
#   detail <- coarse_res |>
#     distinct(als_policy, extra_policy_days, .keep_all = TRUE) |>
#     arrange(extra_policy_days)
#   
#   list(min_extra_days = min_day, detail = detail)
# }
# 
# # ------------------------------- Run ---------------------------------------#
# results <- imap(policy_beta_sets, function(beta_vec, als_label) {
#   
#   # Build template spec once per ALS
#   spec_template <- spec_base |>
#     set_vaccination_observed() |>
#     spec_with_beta_slot(beta_vec)
#   
#   find_min_extra_days_fast(spec_template, beta_vec, als_label)
# })
# 
# min_duration_summary <- tibble(
#   als_policy = factor(names(policy_beta_sets), levels = names(policy_beta_sets)),
#   min_extra_days = map_int(results, ~ ifelse(is.na(.x$min_extra_days), NA_integer_, .x$min_extra_days)),
#   note = map_chr(results, ~ ifelse(is.na(.x$min_extra_days),
#                                    "Target not achieved within search range.",
#                                    ""))
# )
# 
# policy_duration_detail <- bind_rows(map(results, "detail")) |>
#   mutate(als_policy = factor(als_policy, levels = names(policy_beta_sets))) |>
#   arrange(als_policy, extra_policy_days)
# 
# print(min_duration_summary)
# write.csv(min_duration_summary, csv_sum, row.names = FALSE)
# write.csv(policy_duration_detail, csv_det, row.names = FALSE)
# 
# # ------------------------------- Plots -------------------------------------#
# p_peak_curve <- ggplot(policy_duration_detail, aes(extra_policy_days, peak_daily_infections)) +
#   geom_hline(yintercept = infection_peak_target, linetype = 4, linewidth = 0.9) +
#   geom_line(linewidth = 1.0) +
#   facet_wrap(~als_policy, nrow = 1, scales = "free_y") +
#   scale_x_continuous(breaks = pretty_breaks(7)) +
#   scale_y_continuous(labels = comma) +
#   labs(
#     title = paste0("Policy duration required to keep peak daily infections \u2264 ", comma(infection_peak_target)),
#     x = "Extra days maintaining policy beyond March 14, 2022",
#     y = "Peak daily infections (Dec 17, 2021–May 22, 2022)"
#   ) +
#   theme_clean() +
#   theme(
#     plot.title = element_text(size = 16, hjust = 0.5, face = "plain"),
#     strip.text = element_text(size = 13, face = "plain"),
#     axis.text  = element_text(size = 11),
#     axis.title = element_text(size = 12)
#   )
# 
# save_png(file.path(fig_dir, "peak_infections_vs_extra_days.png"), p_peak_curve, 12, 4.2)
# 
# p_min_days <- ggplot(min_duration_summary, aes(als_policy, min_extra_days)) +
#   geom_col(width = 0.7, na.rm = TRUE) +
#   geom_text(aes(label = ifelse(is.na(min_extra_days), "Not achieved", paste0(min_extra_days, " days"))),
#             vjust = -0.4, size = 4) +
#   scale_y_continuous(breaks = pretty_breaks(7), expand = expansion(mult = c(0, 0.15))) +
#   labs(
#     title = paste0("Minimum extra days needed (target: peak daily infections \u2264 ", comma(infection_peak_target), ")"),
#     x = NULL, y = "Minimum extra days"
#   ) +
#   theme_clean() +
#   theme(
#     plot.title  = element_text(size = 16, hjust = 0.5, face = "plain"),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 11),
#     axis.title.y = element_text(size = 12)
#   )
# 
# save_png(file.path(fig_dir, "min_extra_days_by_policy.png"), p_min_days, 7.5, 4.5)
# 
# message("Done. Wrote:\n- ", csv_sum, "\n- ", csv_det,
#         "\nSaved figures to:\n- ", fig_dir)


# ============================== FIXED + INFORMATIVE PLOTS ==============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggthemes)
  library(scales)
})

plot_dir <- "../figures/duration_infection_target"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, p, width, height, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

# --- core constants used in your run (edit ONLY if different) ---
extra_days_max <- max(policy_duration_detail$extra_policy_days, na.rm = TRUE)
infection_peak_target <- 1261

policy_levels <- names(policy_beta_sets)

detail_df <- policy_duration_detail |>
  mutate(
    als_policy = factor(als_policy, levels = policy_levels),
    meets_target = as.logical(meets_target)
  )

# --- recompute minimum extra days from the detail (robust; avoids NA-only plots) ---
min_days_df <- detail_df |>
  group_by(als_policy) |>
  summarise(
    min_extra_days = if (any(meets_target, na.rm = TRUE)) min(extra_policy_days[meets_target], na.rm = TRUE) else NA_integer_,
    .groups = "drop"
  ) |>
  mutate(
    achieved = !is.na(min_extra_days),
    # value used only for plotting bars (so "not achieved" still shows)
    bar_height_days = ifelse(achieved, min_extra_days, extra_days_max),
    bar_label = ifelse(achieved, paste0(min_extra_days, " days"), paste0("> ", extra_days_max, " days")),
    als_policy = factor(als_policy, levels = policy_levels)
  )

# =======================================================================================
# 1) Curve plot: peak infections vs extra days + target + feasible points + min-day line
# =======================================================================================
vline_df <- min_days_df |>
  filter(achieved) |>
  transmute(als_policy, xint = min_extra_days)

not_achieved_df <- min_days_df |>
  filter(!achieved) |>
  transmute(
    als_policy,
    x = extra_days_max * 0.55,
    y = infection_peak_target * 1.05,
    lab = paste0("Not achieved\nwithin ", extra_days_max, " days")
  )

p_curve <- ggplot(detail_df, aes(extra_policy_days, peak_daily_infections)) +
  geom_hline(yintercept = infection_peak_target, linetype = 4, linewidth = 0.9) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(shape = meets_target), size = 2.0, alpha = 0.9) +
  geom_vline(data = vline_df, aes(xintercept = xint), linetype = 3, linewidth = 0.9) +
  geom_text(
    data = not_achieved_df,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~als_policy, nrow = 1, scales = "free_y") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = NULL,
                     labels = c("Does not meet target", "Meets target")) +
  scale_x_continuous(breaks = pretty_breaks(7)) +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste0("Policy duration required to keep peak daily infections \u2264 ", comma(infection_peak_target)),
    x = "Extra days maintaining policy beyond March 14, 2022",
    y = "Peak daily infections (Dec 17, 2021–May 22, 2022)"
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 13, face = "plain"),
    axis.text  = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "top",
    plot.margin = margin(10, 18, 10, 12)
  )

save_png(file.path(plot_dir, "peak_infections_vs_extra_days_fixed.png"), p_curve, width = 13, height = 4.6)

# =======================================================================================
# 2) Bar plot: minimum extra days (shows “> max” when not achieved)
# =======================================================================================
p_min_days <- ggplot(min_days_df, aes(als_policy, bar_height_days, fill = achieved)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = bar_label), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(`TRUE` = "grey40", `FALSE` = "grey80"), guide = "none") +
  scale_y_continuous(
    breaks = pretty_breaks(7),
    limits = c(0, extra_days_max * 1.15),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = paste0("Minimum extra days needed (target: peak daily infections \u2264 ", comma(infection_peak_target), ")"),
    x = NULL,
    y = "Minimum extra days"
  ) +
  theme_clean() +
  theme(
    plot.title  = element_text(size = 16, hjust = 0.5, face = "plain"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    plot.margin = margin(10, 18, 10, 12)
  )

save_png(file.path(plot_dir, "min_extra_days_by_policy_fixed.png"), p_min_days, width = 7.8, height = 4.8)

# =======================================================================================
# 3) Bar plot: peak infections at 0 days vs at minimum required extension
# =======================================================================================
peak_at_zero <- detail_df |>
  filter(extra_policy_days == 0) |>
  select(als_policy, peak0 = peak_daily_infections)

peak_at_min <- detail_df |>
  inner_join(min_days_df |> filter(achieved) |> select(als_policy, min_extra_days),
             by = "als_policy") |>
  filter(extra_policy_days == min_extra_days) |>
  select(als_policy, peak_min = peak_daily_infections)

peak_compare <- min_days_df |>
  left_join(peak_at_zero, by = "als_policy") |>
  left_join(peak_at_min,  by = "als_policy") |>
  mutate(
    peak_min = ifelse(is.na(peak_min), NA_real_, peak_min)
  ) |>
  transmute(
    als_policy,
    `Peak at 0 extra days` = peak0,
    `Peak at minimum extra days` = peak_min
  ) |>
  pivot_longer(-als_policy, names_to = "scenario", values_to = "peak") |>
  mutate(
    scenario = factor(scenario, levels = c("Peak at 0 extra days", "Peak at minimum extra days"))
  )

p_peak_bars <- ggplot(peak_compare, aes(als_policy, peak, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = infection_peak_target, linetype = 4, linewidth = 0.9) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Peak daily infections before and after the minimum required extension",
    x = NULL,
    y = "Peak daily infections",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "plain"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    plot.margin = margin(10, 18, 10, 12)
  )

save_png(file.path(plot_dir, "peak_infections_before_after_min_extension.png"), p_peak_bars, width = 10.5, height = 4.8)

message("Saved improved figures to: ", plot_dir)
