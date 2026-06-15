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

# ---- Inputs ----
calibrator     <- readRDS("calibrate.rds")
seroprevdata   <- readRDS("seroprevdata.rds")

beta_ALS4only  <- rdsRead("beta_als4")
beta_noALS     <- rdsRead("policy_neutral")

upper_plot_time <- 300
offset0         <- 150
anchor_start    <- as.Date("2021-12-15")
anchor_end      <- as.Date("2022-05-22")

# ---- Model constants ----
kappa1 <- params[["kappa1"]]
kappa2 <- params[["kappa2"]]
kappa3 <- params[["kappa3"]]

gamma_i <- params[["gamma_i"]]
gamma_a <- params[["gamma_a"]]

mu      <- params[["mu"]]
zeta    <- params[["zeta"]]

p1 <- 0.15
p2 <- 0.85
p3 <- 0

bracket_term   <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)
susceptibility <- p1 * kappa1 + p2 * kappa2 + p3 * kappa3
mult_const     <- bracket_term * susceptibility

# ------------------------------------------------------------
# Helper: vector beta(t) -> Rc(t)
# ------------------------------------------------------------

beta_vec_to_Rc_df <- function(beta_vec,
                              offset0,
                              anchor_start,
                              anchor_end,
                              mult_const) {
  
  tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number()) |>
    filter(time >= offset0) |>
    mutate(
      date = seq.Date(
        from = anchor_start,
        by   = "day",
        length.out = n()
      ),
      value = beta * mult_const
    ) |>
    filter(date >= anchor_start,
           date <= anchor_end) |>
    transmute(date, value)
}

# ------------------------------------------------------------
# Actual fitted Rc(t) with credible interval
# ------------------------------------------------------------

sims_Rc <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  filter(matrix == "beta_thing",
         time >= offset0) |>
  mutate(
    date = seq.Date(
      from = anchor_start,
      by   = "day",
      length.out = n()
    )
  ) |>
  filter(date >= anchor_start,
         date <= anchor_end) |>
  transmute(
    date,
    value     = value * mult_const,
    conf.low  = conf.low * mult_const,
    conf.high = conf.high * mult_const
  )

# ------------------------------------------------------------
# Counterfactual Rc(t) curves
# ------------------------------------------------------------

Rc_noALS_df <- beta_vec_to_Rc_df(
  beta_vec      = beta_noALS,
  offset0       = offset0,
  anchor_start  = anchor_start,
  anchor_end    = anchor_end,
  mult_const    = mult_const
)

Rc_ALS4only_df <- beta_vec_to_Rc_df(
  beta_vec      = beta_ALS4only,
  offset0       = offset0,
  anchor_start  = anchor_start,
  anchor_end    = anchor_end,
  mult_const    = mult_const
)

# ------------------------------------------------------------
# Add credible intervals to counterfactual curves
# Uses relative uncertainty from fitted beta_thing
# ------------------------------------------------------------

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
# Plot 1: Rc(t)
# ============================================================

p1 <- ggplot() +
  
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
    xintercept = as.Date("2022-03-15"),
    linetype = "dashed",
    linewidth = 1.5
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c(
      "Actual NPIs",
      "No NPIs",
      "Strict NPIs"
    ),
    values = c(
      "Actual NPIs" = "black",
      "No NPIs" = "red",
      "Strict NPIs" = "#00C5CD"
    )
  )+
  # scale_color_manual(
  #   name = NULL,
  #   values = c(
  #     "Actual NPIs" = "black",
  #     "Strict NPIs" = "#00C5CD",
  #     "No NPIs" = "red"
  #   )
  # ) +
  
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
    legend.position = c(0.25, 0.85),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 30, color = "black", hjust = 0.5, face = "plain"),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

dir.create("../figures", showWarnings = FALSE, recursive = TRUE)

png(
  "../figures/control_reprod_numb.png",
  width  = 5000,
  height = 1600,
  res    = 300,
  bg     = "white",
  type   = "cairo"
)

print(p1)

dev.off()
