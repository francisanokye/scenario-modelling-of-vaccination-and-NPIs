# suppressPackageStartupMessages({
#   library(dplyr)
#   library(lubridate)
#   library(ggplot2)
#   library(tidyr)
#   library(zoo)
#   library(ggthemes)
#   library(patchwork)
#   library(shellpipes)
#   library(tidyverse)
#   library(macpan2)
#   library(grid)
# })
# 
# set.seed(2025)
# options(macpan2_log_dir = ".")
# loadEnvironments()
# 
# # ---- Inputs ----
# calibrator     <- readRDS("calibrate.rds")
# seroprevdata   <- readRDS("seroprevdata.rds")
# 
# beta_ALS4only  <- rdsRead("beta_als4")
# beta_noALS     <- rdsRead("policy_neutral")
# 
# upper_plot_time <- 300
# offset0         <- 150
# anchor_start    <- as.Date("2021-12-15")
# anchor_end      <- as.Date("2022-05-22")
# 
# # ---- Model constants ----
# kappa1 <- params[["kappa1"]]
# kappa2 <- params[["kappa2"]]
# kappa3 <- params[["kappa3"]]
# 
# gamma_i <- params[["gamma_i"]]
# gamma_a <- params[["gamma_a"]]
# 
# mu      <- params[["mu"]]
# zeta    <- params[["zeta"]]
# 
# # ============================================================
# # Rc(t) susceptibility term -- TIME-VARYING cohort proportions
# #
# # The control reproduction number (see model description) is
# #   Rc(t) = beta(t) * ( sum_j p_j(t) * kappa_j ) * bracket_term,
# # where p_j(t) = S_j(t)/N are the *time-varying* cohort proportions
# # (p1 = S1/N, p2 = V2/N, p3 = V3/N), NOT constants.
# #
# # PREVIOUS BUG: this script hardcoded p1 = 0.15, p2 = 0.85, p3 = 0,
# # collapsing the susceptibility term to a single constant. That made
# # Rc(t) just a rescaled copy of beta(t), erasing the effect of the
# # shifting cohort composition over the Omicron wave (as people move
# # S1 -> V2 -> V3, the susceptibility-weighted sum changes because
# # kappa3 < kappa2 < kappa1). Setting p3 = 0 additionally zeroed out the
# # entire booster cohort's contribution. The proportions are now pulled
# # live from the fitted model trajectory.
# #
# # >>> CALIBRATION DEPENDENCY (IMPORTANT) <<<
# # These p1/p2/p3 come from the calibrator trajectory. If the calibration
# # was run with the booster flow inactive (V3 never filling -- the known
# # issue where cum_booster is flat 0 and p3 ~ 1.9e-6), then p3(t) here
# # will still be ~0 and this fix will have little visible effect. The
# # upstream calibration must actually fill V3 (real booster rollout fed
# # into timevar_spec) for the time-varying p3 to be meaningful. Once the
# # calibration is corrected, this script needs no further change -- it
# # will pick up the corrected p3(t) automatically.
# #
# # DESIGN CHOICE: the counterfactual Rc curves (No NPIs / Strict NPIs)
# # use the SAME p1/p2/p3 trajectory from the actual fitted model, so the
# # counterfactual curves isolate the beta/NPI effect on Rc without
# # confounding from cohort-composition differences (and because the
# # counterfactual beta series are bare vectors with no compartment
# # trajectories of their own).
# # ============================================================
# 
# bracket_term <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)
# 
# # pull time-varying p1/p2/p3 from the fitted model, aligned to dates
# prop_df <- calibrator |>
#   mp_trajectory() |>
#   filter(matrix %in% c("p1", "p2", "p3"),
#          time >= offset0) |>
#   select(time, matrix, value) |>
#   pivot_wider(names_from = matrix, values_from = value) |>
#   arrange(time) |>
#   mutate(
#     date = seq.Date(
#       from = anchor_start,
#       by   = "day",
#       length.out = n()
#     ),
#     susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3,
#     mult_t         = bracket_term * susceptibility
#   ) |>
#   filter(date >= anchor_start, date <= anchor_end) |>
#   select(date, p1, p2, p3, susceptibility, mult_t)
# 
# # ------------------------------------------------------------
# # Helper: vector beta(t) -> Rc(t), using the TIME-VARYING mult_t
# # (joined by date) instead of a constant multiplier
# # ------------------------------------------------------------
# 
# beta_vec_to_Rc_df <- function(beta_vec,
#                               offset0,
#                               anchor_start,
#                               anchor_end,
#                               prop_df) {
#   
#   tibble(beta = as.numeric(beta_vec)) |>
#     mutate(time = row_number()) |>
#     filter(time >= offset0) |>
#     mutate(
#       date = seq.Date(
#         from = anchor_start,
#         by   = "day",
#         length.out = n()
#       )
#     ) |>
#     filter(date >= anchor_start, date <= anchor_end) |>
#     left_join(prop_df, by = "date") |>
#     transmute(date, value = beta * mult_t)
# }
# 
# # ------------------------------------------------------------
# # Actual fitted Rc(t) with credible interval
# # (fitted beta_thing scaled by the time-varying mult_t, joined by date)
# # ------------------------------------------------------------
# 
# sims_Rc <- calibrator |>
#   mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
#   filter(matrix == "beta_thing",
#          time >= offset0) |>
#   mutate(
#     date = seq.Date(
#       from = anchor_start,
#       by   = "day",
#       length.out = n()
#     )
#   ) |>
#   filter(date >= anchor_start,
#          date <= anchor_end) |>
#   left_join(prop_df, by = "date") |>
#   transmute(
#     date,
#     value     = value     * mult_t,
#     conf.low  = conf.low  * mult_t,
#     conf.high = conf.high * mult_t
#   )
# 
# # ------------------------------------------------------------
# # Counterfactual Rc(t) curves
# # ------------------------------------------------------------
# 
# Rc_noALS_df <- beta_vec_to_Rc_df(
#   beta_vec      = beta_noALS,
#   offset0       = offset0,
#   anchor_start  = anchor_start,
#   anchor_end    = anchor_end,
#   prop_df       = prop_df
# )
# 
# Rc_ALS4only_df <- beta_vec_to_Rc_df(
#   beta_vec      = beta_ALS4only,
#   offset0       = offset0,
#   anchor_start  = anchor_start,
#   anchor_end    = anchor_end,
#   prop_df       = prop_df
# )
# 
# # ------------------------------------------------------------
# # Add credible intervals to counterfactual curves
# # Uses relative uncertainty from fitted beta_thing
# # ------------------------------------------------------------
# 
# ci_ratio_df <- sims_Rc |>
#   mutate(
#     low_ratio  = conf.low / value,
#     high_ratio = conf.high / value
#   ) |>
#   select(date, low_ratio, high_ratio)
# 
# Rc_noALS_df <- Rc_noALS_df |>
#   left_join(ci_ratio_df, by = "date") |>
#   mutate(
#     conf.low  = value * low_ratio,
#     conf.high = value * high_ratio
#   )
# 
# Rc_ALS4only_df <- Rc_ALS4only_df |>
#   left_join(ci_ratio_df, by = "date") |>
#   mutate(
#     conf.low  = value * low_ratio,
#     conf.high = value * high_ratio
#   )
# 
# 
# # ============================================================
# # Rc(t)
# # ============================================================
# 
# p_rc <- ggplot() +
#   
#   geom_ribbon(
#     data = sims_Rc,
#     aes(x = date, ymin = conf.low, ymax = conf.high),
#     fill = "black",
#     alpha = 0.22,
#     show.legend = FALSE
#   ) +
#   geom_line(
#     data = sims_Rc,
#     aes(x = date, y = value, color = "Actual NPIs"),
#     linewidth = 2
#   ) +
#   geom_line(
#     data = Rc_noALS_df,
#     aes(x = date, y = value, color = "No NPIs"),
#     linewidth = 2
#   ) +
#   geom_line(
#     data = Rc_ALS4only_df,
#     aes(x = date, y = value, color = "Strict NPIs"),
#     linewidth = 2
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     linewidth = 1.2
#   ) +
#   geom_vline(
#     xintercept = as.Date("2022-03-15"),
#     linetype = "dashed",
#     linewidth = 1.5
#   ) +
#   scale_color_manual(
#     name = NULL,
#     breaks = c(
#       "Actual NPIs",
#       "No NPIs",
#       "Strict NPIs"
#     ),
#     values = c(
#       "Actual NPIs" = "black",
#       "No NPIs" = "red",
#       "Strict NPIs" = "#00C5CD"
#     )
#   )+
#   
#   scale_x_date(
#     limits = c(anchor_start, anchor_end),
#     date_breaks = "2 week",
#     date_labels = "%b %d",
#     expand = c(0, 0.5)
#   ) +
#   
#   labs(
#     title = "Reproduction number for counterfactual scenarios",
#     y = "Reproduction number",
#     x = NULL,
#     tag = "A"
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 22),
#     axis.text.y = element_text(size = 30),
#     axis.title.y = element_text(size = 30, color = "black"),
#     legend.position = c(0.85, 0.45),
#     legend.justification = c(0.5, 0.5),
#     legend.direction = "vertical",
#     legend.box = "vertical",
#     legend.title = element_text(size = 25),
#     legend.text = element_text(size = 25),
#     legend.background = element_blank(),
#     plot.title = element_text(size = 30, color = "black", hjust = 0.5, face = "plain"),
#     plot.background = element_blank(),
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
# 
# png(
#   "../figures/control_reprod_numb.png",
#   width  = 5000,
#   height = 1600,
#   res    = 300,
#   bg     = "white",
#   type   = "cairo"
# )
# 
# print(p_rc)
# 
# dev.off()
# 
# 
# 
# # ------------------------------------------------------------
# # Plot B: transmission rate beta(t)
# #
# # FIXES: the previous version of this section referenced several objects
# # that were never defined in this script and caused
# #   "Error: object 'beta_all_df' not found":
# #   - beta_all_df : the long-format data frame of the three beta(t) series
# #                   (Actual / No NPIs / Strict NPIs), now built below from
# #                   the fitted beta_thing and the two counterfactual beta
# #                   vectors (beta_noALS, beta_ALS4only).
# #   - cut_date    : the release-date vertical line, now defined explicitly
# #                   (matches the Mar 15 2022 line used in the Rc panel).
# #   - analysis_start / analysis_end : the Plot B scale_x_date used these
# #                   names, but the rest of the script uses anchor_start /
# #                   anchor_end; now reconciled to anchor_start/anchor_end.
# #   - p_rc        : the combine step referenced p_rc, but the Rc plot
# #                   object was named p1; the Rc plot is now named p_rc.
# # ------------------------------------------------------------
# 
# cut_date <- as.Date("2022-03-15")
# 
# # Fitted (Actual NPIs) beta(t), aligned to dates the same way as sims_Rc.
# # Sourced from mp_trajectory_sd (the SAME call used for the panel A Actual
# # curve) rather than mp_trajectory, to guarantee identical row count and
# # date alignment between the two panels -- if the two trajectory functions
# # returned beta_thing with even a slightly different length, the seq.Date
# # alignment could shift or drop the Actual rows, which would make the
# # black Actual beta line disappear from panel B.
# beta_actual_df <- calibrator |>
#   mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
#   filter(matrix == "beta_thing",
#          time >= offset0) |>
#   arrange(time) |>
#   mutate(
#     date = seq.Date(from = anchor_start, by = "day", length.out = n())
#   ) |>
#   filter(date >= anchor_start, date <= anchor_end) |>
#   transmute(date, beta = value, scenario = "Actual NPIs")
# 
# # Helper: turn a bare beta vector into a dated, labelled data frame
# beta_vec_to_df <- function(beta_vec, scenario_label) {
#   tibble(beta = as.numeric(beta_vec)) |>
#     mutate(time = row_number()) |>
#     filter(time >= offset0) |>
#     mutate(
#       date = seq.Date(from = anchor_start, by = "day", length.out = n())
#     ) |>
#     filter(date >= anchor_start, date <= anchor_end) |>
#     transmute(date, beta, scenario = scenario_label)
# }
# 
# beta_noALS_df    <- beta_vec_to_df(beta_noALS,    "No NPIs")
# beta_ALS4only_df <- beta_vec_to_df(beta_ALS4only, "Strict NPIs")
# 
# beta_all_df <- bind_rows(
#   beta_actual_df,
#   beta_noALS_df,
#   beta_ALS4only_df
# ) |>
#   mutate(
#     scenario = factor(
#       scenario,
#       levels = c("Actual NPIs", "No NPIs", "Strict NPIs")
#     )
#   )
# 
# p_beta <- ggplot() +
#   geom_line(
#     data = beta_all_df |> dplyr::filter(scenario != "Actual NPIs"),
#     aes(x = date, y = beta, color = scenario),
#     linewidth = 2.5
#   ) +
#   # Draw Actual NPIs as its own layer LAST so it is always rendered on
#   # top and never painted over by the counterfactual lines (mirrors how
#   # panel A draws the Actual Rc curve as a separate top layer). Without
#   # this, the black Actual beta line could be hidden underneath the
#   # No NPIs / Strict NPIs lines wherever they overlap.
#   geom_line(
#     data = beta_all_df |> dplyr::filter(scenario == "Actual NPIs"),
#     aes(x = date, y = beta, color = scenario),
#     linewidth = 2.5
#   ) +
#   geom_vline(
#     xintercept = cut_date,
#     linetype = "dashed",
#     linewidth = 2
#   ) +
#   scale_color_manual(
#     name = NULL,
#     breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
#     values = c(
#       "Actual NPIs" = "black",
#       "No NPIs" = "red",
#       "Strict NPIs" = "#00C5CD"
#     )
#   ) +
#   scale_x_date(
#     limits       = c(anchor_start, anchor_end),
#     date_breaks  = "2 weeks",
#     date_labels  = "%d %b %Y",
#     expand       = c(0.01, 0)
#   ) +
#   labs(
#     title = "Transmission rate",
#     y = expression(beta(t)),
#     x = NULL,
#     tag = "B"
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 30),
#     axis.title.y = element_text(size = 30, color = "black"),
#     legend.position = c(0.85, 0.45),
#     legend.justification = c(0.5, 0.5),
#     legend.direction = "vertical",
#     legend.box = "vertical",
#     legend.title = element_text(size = 25),
#     legend.text = element_text(size = 25),
#     legend.background = element_blank(),
#     plot.title = element_text(size = 30, color = "black", hjust = 0.5, face = "plain"),
#     plot.background = element_blank(),
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# # ------------------------------------------------------------
# # combine and save
# # ------------------------------------------------------------
# 
# combined_plot <- p_rc + p_beta +
#   plot_layout(widths = c(1.5, 1))
# 
# dir.create("../figures", showWarnings = FALSE, recursive = TRUE)
# 
# png(
#   "../figures/reprod_numb_with_beta.png",
#   width  = 7000,
#   height = 2000,
#   res    = 300,
#   bg     = "white",
#   type   = "cairo"
# )
# 
# print(combined_plot)
# 
# dev.off()





suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(zoo)
  library(ggthemes)
  library(patchwork)
  library(shellpipes)
  library(tidyverse)
  library(macpan2)
  library(grid)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ---- inputs ----
calibrator     <- readRDS("calibrate.rds")
seroprevdata   <- readRDS("seroprevdata.rds")

beta_ALS4only  <- rdsRead("beta_als4")
beta_noALS     <- rdsRead("policy_neutral")

upper_plot_time <- 300
offset0         <- 150
anchor_start    <- as.Date("2021-12-15")
anchor_end      <- as.Date("2022-05-22")
cut_date        <- as.Date("2022-03-15")

# ---- model constants ----
kappa1 <- params[["kappa1"]]
kappa2 <- params[["kappa2"]]
kappa3 <- params[["kappa3"]]

gamma_i <- params[["gamma_i"]]
gamma_a <- params[["gamma_a"]]

mu   <- params[["mu"]]
zeta <- params[["zeta"]]

bracket_term <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)

# ---- time-varying vaccination-cohort proportions ----
prop_df <- calibrator |>
  mp_trajectory() |>
  filter(
    matrix %in% c("p1", "p2", "p3"),
    time >= offset0
  ) |>
  select(time, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value) |>
  arrange(time) |>
  mutate(
    date = seq.Date(
      from = anchor_start,
      by   = "day",
      length.out = n()
    ),
    susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3,
    psum = p1 + p2 + p3,
    mult_t = bracket_term * susceptibility
  ) |>
  filter(date >= anchor_start, date <= anchor_end) |>
  select(date, p1, p2, p3, psum, susceptibility, mult_t)

weighted_susc_df <- prop_df |>
  select(date, weighted_susc = susceptibility, psum)

# ---- helper: beta vector to Rc dataframe ----
beta_vec_to_Rc_df <- function(beta_vec,
                              offset0,
                              anchor_start,
                              anchor_end,
                              prop_df) {
  tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number()) |>
    filter(time >= offset0) |>
    mutate(
      date = seq.Date(
        from = anchor_start,
        by = "day",
        length.out = n()
      )
    ) |>
    filter(date >= anchor_start, date <= anchor_end) |>
    left_join(prop_df, by = "date") |>
    transmute(date, value = beta * mult_t)
}

# ---- actual Rc(t) with credible interval ----
sims_Rc <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  filter(
    matrix == "beta_thing",
    time >= offset0
  ) |>
  mutate(
    date = seq.Date(
      from = anchor_start,
      by = "day",
      length.out = n()
    )
  ) |>
  filter(date >= anchor_start, date <= anchor_end) |>
  left_join(prop_df, by = "date") |>
  transmute(
    date,
    value     = value * mult_t,
    conf.low  = conf.low * mult_t,
    conf.high = conf.high * mult_t
  )

# ---- counterfactual Rc(t) curves ----
Rc_noALS_df <- beta_vec_to_Rc_df(
  beta_vec = beta_noALS,
  offset0 = offset0,
  anchor_start = anchor_start,
  anchor_end = anchor_end,
  prop_df = prop_df
)

Rc_ALS4only_df <- beta_vec_to_Rc_df(
  beta_vec = beta_ALS4only,
  offset0 = offset0,
  anchor_start = anchor_start,
  anchor_end = anchor_end,
  prop_df = prop_df
)

# ---- uncertainty ratios from actual fitted beta/Rc ----
ci_ratio_df <- sims_Rc |>
  mutate(
    low_ratio  = conf.low / value,
    high_ratio = conf.high / value
  ) |>
  select(date, low_ratio, high_ratio)

Rc_noALS_df <- Rc_noALS_df |>
  left_join(ci_ratio_df, by = "date") |>
  mutate(
    conf.low  = value * low_ratio,
    conf.high = value * high_ratio
  )

Rc_ALS4only_df <- Rc_ALS4only_df |>
  left_join(ci_ratio_df, by = "date") |>
  mutate(
    conf.low  = value * low_ratio,
    conf.high = value * high_ratio
  )

# ============================================================
# panel A: Rc(t)
# ============================================================

p_rc <- ggplot() +
  geom_ribbon(
    data = sims_Rc,
    aes(x = date, ymin = conf.low, ymax = conf.high),
    fill = "black",
    alpha = 0.22,
    show.legend = FALSE
  ) +
  geom_line(
    data = sims_Rc,
    aes(x = date, y = value, color = "Actual NPIs"),
    linewidth = 2
  ) +
  geom_line(
    data = Rc_noALS_df,
    aes(x = date, y = value, color = "No NPIs"),
    linewidth = 2
  ) +
  geom_line(
    data = Rc_ALS4only_df,
    aes(x = date, y = value, color = "Strict NPIs"),
    linewidth = 2
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = cut_date,
    linetype = "dashed",
    linewidth = 1.5
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
    values = c(
      "Actual NPIs" = "black",
      "No NPIs" = "red",
      "Strict NPIs" = "#00C5CD"
    )
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 week",
    date_labels = "%b %d",
    expand = c(0, 0.5)
  ) +
  labs(
    title = "Reproduction number for counterfactual scenarios",
    y = "Reproduction number",
    x = NULL,
    tag = "A"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = c(0.85, 0.75),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# panel B: transmission rate beta(t)
# ============================================================

beta_actual_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  filter(
    matrix == "beta_thing",
    time >= offset0
  ) |>
  arrange(time) |>
  mutate(
    date = seq.Date(
      from = anchor_start,
      by = "day",
      length.out = n()
    )
  ) |>
  filter(date >= anchor_start, date <= anchor_end) |>
  transmute(date, beta = value, scenario = "Actual NPIs")

beta_vec_to_df <- function(beta_vec, scenario_label) {
  tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number()) |>
    filter(time >= offset0) |>
    mutate(
      date = seq.Date(
        from = anchor_start,
        by = "day",
        length.out = n()
      )
    ) |>
    filter(date >= anchor_start, date <= anchor_end) |>
    transmute(date, beta, scenario = scenario_label)
}

beta_noALS_df    <- beta_vec_to_df(beta_noALS, "No NPIs")
beta_ALS4only_df <- beta_vec_to_df(beta_ALS4only, "Strict NPIs")

beta_all_df <- bind_rows(
  beta_actual_df,
  beta_noALS_df,
  beta_ALS4only_df
) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("Actual NPIs", "No NPIs", "Strict NPIs")
    )
  )

p_beta <- ggplot() +
  geom_line(
    data = beta_all_df |> filter(scenario != "Actual NPIs"),
    aes(x = date, y = beta, color = scenario),
    linewidth = 2.5
  ) +
  geom_line(
    data = beta_all_df |> filter(scenario == "Actual NPIs"),
    aes(x = date, y = beta, color = scenario),
    linewidth = 2.5
  ) +
  geom_vline(
    xintercept = cut_date,
    linetype = "dashed",
    linewidth = 2
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
    values = c(
      "Actual NPIs" = "black",
      "No NPIs" = "red",
      "Strict NPIs" = "#00C5CD"
    )
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 weeks",
    date_labels = "%d %b %Y",
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Transmission rate",
    y = expression(beta(t)),
    x = NULL,
    tag = "B"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = c(0.85, 0.45),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# panel C: weighted susceptibility
# ============================================================

p_susc <- ggplot(weighted_susc_df, aes(x = date, y = weighted_susc)) +
  geom_line(
    linewidth = 2.5,
    colour = "#1B9E77"
  ) +
  geom_vline(
    xintercept = cut_date,
    linetype = "dashed",
    linewidth = 2
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 weeks",
    date_labels = "%d %b %Y",
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Weighted susceptibility",
    y = "Weighted susceptibility",
    x = NULL,
    tag = "C"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# save separate Rc plot
# ============================================================

dir.create("../figures", showWarnings = FALSE, recursive = TRUE)

png(
  "../figures/control_reprod_numb.png",
  width = 5000,
  height = 1600,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(p_rc)

dev.off()

# ============================================================
# combine and save all panels
# ============================================================

combined_plot <- p_rc + p_beta + p_susc +
  plot_layout(widths = c(1.5, 1, 1))

png(
  "../figures/reprod_numb_with_beta_susceptibility.png",
  width = 9500,
  height = 2000,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(combined_plot)

dev.off()

# optional diagnostic check
print(
  weighted_susc_df |>
    summarise(
      min_psum = min(psum, na.rm = TRUE),
      max_psum = max(psum, na.rm = TRUE),
      min_weighted_susc = min(weighted_susc, na.rm = TRUE),
      max_weighted_susc = max(weighted_susc, na.rm = TRUE)
    )
)