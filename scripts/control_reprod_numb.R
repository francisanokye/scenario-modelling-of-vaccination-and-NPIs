
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
  library(grid)   # unit()
})
#

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ---- Inputs / plotting window ----
calibrator     <- readRDS("calibrate.rds")
seroprevdata   <- readRDS("seroprevdata.rds")

beta_ALS4only  <- rdsRead("beta_als4")         
beta_noALS     <- rdsRead("policy_neutral")    

upper_plot_time <- 300
offset0         <- 150
anchor_start    <- as.Date("2021-12-15")
anchor_end      <- as.Date("2022-05-22")

# ---- Model constants used for Rc(t) scaling ----
# Rc(t) = beta(t) * [ (mu/gamma_i) + ((1-mu)*zeta/gamma_a) ] * [p1*kappa1 + p2*kappa2 + p3*kappa3]
kappa1 <- 1; kappa2 <- 0.91; kappa3 <- 0.3
gamma_i <- 1/7; gamma_a <- 1/10; mu <- 0.678; zeta <- 0.75
p1 <- 0.15; p2 <- 0.85; p3 <- 0

bracket_term   <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)
susceptibility <- p1 * kappa1 + p2 * kappa2 + p3 * kappa3
mult_const     <- bracket_term * susceptibility

# ------------------------------------------------------------
# numeric beta(t) vector -> Rc(t) tibble + 95% CI bands
# - CI bands are computed from Rc values using quantiles (0.025, 0.975)
# ------------------------------------------------------------
beta_vec_to_Rc_df <- function(beta_vec,
                              offset0, upper_plot_time,
                              anchor_start, anchor_end,
                              mult_const, z = qnorm(0.975)) {
  
  df <- tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number()) |>
    filter(time >= offset0) |>
    mutate(date = seq.Date(from = anchor_start, by = "day", length.out = n())) |>
    mutate(Rc = beta * mult_const, 
           beta_sd = sd(Rc, na.rm = TRUE),
           beta_lower  = Rc - z * beta_sd,
           beta_upper  = Rc + z * beta_sd
           
           ) |>
    transmute(
      date,
      value = Rc,
      lower = beta_lower,
      upper = beta_upper
    ) |>
    filter(date >= anchor_start, date <= anchor_end) 
}

# ---- Convert No-ALS and ALS-4-only beta(t) -> Rc(t) with bands ----
Rc_noALS_df    <- beta_vec_to_Rc_df(beta_noALS,    offset0, upper_plot_time, anchor_start, anchor_end, mult_const)
Rc_ALS4only_df <- beta_vec_to_Rc_df(beta_ALS4only, offset0, upper_plot_time, anchor_start, anchor_end, mult_const)

# ------------------------------------------------------------
# Calibrated Rc(t) from fitted beta_thing (time-varying CI)
# ------------------------------------------------------------
sims_Rc <- calibrator %>%
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) %>%
  filter(matrix == "beta_thing", time >= offset0) %>%
  mutate(date = seq.Date(from = anchor_start, by = "day", length.out = n())) %>%
  filter(date >= anchor_start, date <= anchor_end) %>%
  transmute(
    date,
    value     = value     * mult_const,
    conf.low  = conf.low  * mult_const,
    conf.high = conf.high * mult_const
  )

# ------------------------------------------------------------
# ALS-4 shading only (remove all other shaded phases)
# NL ALS-4 period in your earlier shading table: 2022-01-08 to 2022-02-07
# ------------------------------------------------------------
als4_shade <- tibble(
  xmin = as.Date("2022-01-08"),
  xmax = as.Date("2022-02-07")
)

als_fill_colors <- c("ALS-4" = "#FFD580")

# ---- K–12 timeline panel 
bracket_df <- tibble(
  xmin      = as.Date(c("2021-12-15","2021-12-20","2022-01-26")),
  xmax      = as.Date(c("2021-12-19","2022-01-25","2022-05-22")),
  sch_label = c("","K-12 Schools Closed","K-12 Schools Open"),
  seg_col   = c("navy","red","navy"),
  lab_col   = c("navy","red","navy")
) %>%
  mutate(x_label = as.Date((as.numeric(xmin) + as.numeric(xmax))/2, origin = "1970-01-01"))

# ---- booster expansion timeline panel 
booster_expansion_df <- tibble(
  xmin      = as.Date(c("2021-12-15")),
  xmax      = as.Date(c("2022-05-22")),
  sch_label = c("Booster Eligibility Expansion For All Adults"),
  seg_col   = c("darkgreen"),
  lab_col   = c("darkgreen")
) %>%
  mutate(x_label = as.Date((as.numeric(xmin) + as.numeric(xmax))/2, origin = "1970-01-01"))


# ============================================================
# Plot 1: Rc(t) curves + CI + ALS-4 shading only
# ============================================================
p1 <- ggplot() +
  # Calibrated Rc(t) ribbon (time-varying CI)
  geom_ribbon(
    data = sims_Rc,
    aes(x = date, ymin = conf.low, ymax = conf.high),
    fill = "black",
    alpha = 0.25,
    show.legend = FALSE
  ) +
  # No-ALS ribbon (quantile band AFTER Rc conversion)
  geom_ribbon(
    data = Rc_noALS_df,
    aes(x = date, ymin = lower, ymax = upper),
    fill = "red",
    alpha = 0.20,
    show.legend = FALSE
  ) +
  # ALS-4 only ribbon (quantile band AFTER Rc conversion)
  geom_ribbon(
    data = Rc_ALS4only_df,
    aes(x = date, ymin = lower, ymax = upper),
    fill = "blue",
    alpha = 0.20,
    show.legend = FALSE
  ) +
  # Lines
  geom_line(
    data = sims_Rc,
    aes(x = date, y = value, color = "Actual ALS"),
    linewidth = 2
  ) +
  geom_line(
    data = Rc_noALS_df,
    aes(x = date, y = value, color = "No-ALS"),
    linewidth = 2
  ) +
  geom_line(
    data = Rc_ALS4only_df,
    aes(x = date, y = value, color = "ALS-4 Only"),
    linewidth = 2
  ) +
  # Threshold line
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2) +
  # Scales / labels
  scale_color_manual(
    name = NULL,
    values = c(
      "Actual ALS"  = "black",
      "No-ALS"      = "red",
      "ALS-4 Only"  = "blue"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = als_fill_colors,
    breaks = "ALS-4",
    guide  = guide_legend(override.aes = list(alpha = 0.6))
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 week",
    date_labels = "%b %d",
    expand = c(0, 0.5)
  ) +
  geom_vline(
             aes(xintercept = as.Date("2022-03-15")),
             linetype = "dashed",
             linewidth = 1.5) +
  labs(
    title = expression("Time-varying control reproduction number for counterfactual scenarios"),
    y     = "Control reproduction number",
    x     = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 22),
    #axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 30),
    axis.title.y     = element_text(size = 30, color = "black"),
    legend.position  = c(0.78,0.85),
    legend.justification = c(0.5, 0.5),
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 25),
    legend.text      = element_text(size = 25),
    legend.background = element_blank(),
    plot.title       = element_text(size = 30, color = "black", hjust = 0.5),
    plot.background  = element_blank()
  ) +
  labs(tag = "A") +
    theme(
      plot.tag = element_text(size = 35, face = "bold", colour = "black"),
      plot.tag.position = c(0.01, 0.98)   # top-left inside panel
    )

# ============================================================
# Plot 2: K–12 bracket timeline
# ============================================================

p2 <- ggplot() +
  geom_segment(
    data = bracket_df,
    aes(x = xmin, xend = xmax, y = 1, yend = 1),
    colour = c("navy", "brown", "navy"),
    linewidth = 1.5,
    arrow = arrow(angle = 90, ends = "both", length = unit(0.5, "cm"))
  ) +
  geom_text(
    data = bracket_df,
    aes(x = x_label, y = 0.25, label = sch_label),
    colour = c("navy", "brown", "navy"),
    size = 9
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 week",
    date_labels = "%b",
    expand = c(0, 0.5)
  ) +
  ylim(-1, 1) +
  theme_void() +
  theme(
    #axis.text.x = element_text(size = 30),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# ============================================================
# Plot 3: Booster Eligibility Expansion 
# ============================================================


# p3 <- ggplot() +
#   geom_segment(
#     data = booster_expansion_df,
#     aes(x = xmin, xend = xmax, y = 1, yend = 1),
#     colour = c("darkgreen"),
#     linewidth = 2,
#     arrow = arrow(angle = 90, ends = "both", length = unit(0.5, "cm"))
#   ) +
#   geom_text(
#     data = booster_expansion_df,
#     aes(x = x_label, y = -0.15, label = sch_label),
#     colour = c("darkgreen"),
#     size = 7
#   ) +
#   scale_x_date(
#     limits = c(anchor_start, anchor_end),
#     date_breaks = "2 week",
#     date_labels = "%b",
#     expand = c(0, 0.5)
#   ) +
#   ylim(-1, 1) +
#   theme_void() +
#   theme(
#     axis.text.x      = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin  = margin(t = 0, r = 0, b = 0, l = 0)
#   )


# ---- Combine ----
gg <- (p1 / p2 ) + plot_layout(heights = c(3.0, 0.3))

# ---- Save (publication settings) ----
png("../figures/control_reprod_numb.png", width = 5000, height = 2500, res = 300, bg = "white", type = "cairo")
gg
dev.off()
