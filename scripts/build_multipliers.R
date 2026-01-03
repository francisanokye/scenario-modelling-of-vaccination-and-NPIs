
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2); library(ggh4x);
  library(scales); library(tibble); library(tidyverse);
  library(broom.mixed); library(zoo); library(patchwork)
})

set.seed(2025)
try(loadEnvironments(), silent = TRUE)

date0    <- as.Date("2021-12-15")
date_end <- as.Date("2022-05-22")
offset0  <- 150
upper_plot_time <- 300
start_date <- date0
basis_cols <- 5
tau_days <- 2.0

# load saved inputs
calibrator   <- rdsRead("calibrate")
seroprevdata <- rdsRead("seroprevdata")
time_steps   <- max(seroprevdata$time)

# initialize calibrated model spec with optimized parameters
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# extract β̂(t) from calibrator
beta_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  dplyr::filter(matrix == "beta_thing", time >= offset0, time <= upper_plot_time) |>
  mutate(date = seq(start_date, by = "day", length.out = length(value)), beta_hat = value) |>
  filter(date >= date0, date <= date_end) |>
  select(date, beta_hat)

# ALS scenarios required
als_cols <- c("als2_only", "als4_only", "no_als_only")
als_fac  <- factor(als_cols, levels = als_cols)

####################################### Counterfactual Scenarios #######################################

in_int  <- function(d, a, b) d >= as.Date(a) & d <= as.Date(b)

base_df  <- tibble(date = beta_df$date, t = as.numeric(date - start_date) + 1)

policy_df <- base_df |>
  mutate(
    ALS = case_when(
      in_int(date,"2021-12-15","2021-12-24") ~ "ALS2",
      in_int(date,"2021-12-25","2022-01-03") ~ "ALS3",
      in_int(date,"2022-01-04","2022-02-07") ~ "ALS4",
      in_int(date,"2022-02-08","2022-03-14") ~ "ALS3",
      in_int(date,"2022-03-15","2022-05-22") ~ "NoALS",
      TRUE ~ "ALS2"
    ),
    school = if_else(in_int(date,"2021-12-20","2022-01-25"), "closed", "open")
  )

# multipliers (rebased so no-als = 1)
# R_ALS4 <- 1.52; R_ALS3 <- 2.12; R_NoALS <- 2.32; R_open <- 2.03; R_closed <- 1.44; R_ALS2 <- sqrt(R_ALS3 * R_NoALS)
R_ALS4 <- 1.91; R_ALS3 <- 2.03; R_NoALS <- 2.26; R_open <- 2.17; R_closed <- 1.49; R_ALS2 <- sqrt(R_ALS3 * R_NoALS)
m_als    <- c(ALS4 = R_ALS4/R_NoALS, ALS3 = R_ALS3/R_NoALS, ALS2 = R_ALS2/R_NoALS, NoALS = 1.00)
m_school <- c(open = 1.0,  closed = R_closed/R_open)

# smooth categorical timelines to series
smooth_step <- function(t, t0, tau) plogis((as.numeric(t - t0))/tau)

make_smooth_series <- function(date, lev, map, tau = tau_days) {
  lev <- as.character(lev)
  seg_idx    <- c(1L, which(lev[-1L] != lev[-length(lev)]) + 1L)
  seg_dates  <- date[seg_idx]
  seg_levels <- lev[seg_idx]
  m <- rep(unname(map[seg_levels[1L]]), length(date))
  if (length(seg_idx) >= 2L) {
    for (j in 2:length(seg_idx)) {
      t0 <- seg_dates[j]
      m_target <- unname(map[seg_levels[j]])
      s <- smooth_step(date, t0, tau)
      m <- m*(1 - s) + m_target*s
    }
  }
  m
}

z <- qnorm(0.975)

df_est <- policy_df %>%
  left_join(beta_df, by = "date") %>%
  mutate(
    m_ALS_hat    = make_smooth_series(date, ALS,    map = m_als,    tau = tau_days),
    m_school_hat = make_smooth_series(date, school, map = m_school, tau = tau_days),
    m_obs_hat    = m_ALS_hat * m_school_hat,
    beta0_hat    = beta_hat / pmax(m_obs_hat, 1e-8),
    beta_recon   = beta0_hat * m_obs_hat,
    beta_hat_sd  = sd(beta_hat, na.rm = TRUE),
    beta0_hat_sd = sd(beta0_hat, na.rm = TRUE),
    beta_hat_l   = beta_hat  - z * beta_hat_sd,
    beta_hat_u   = beta_hat  + z * beta_hat_sd,
    beta0_hat_l  = beta0_hat - z * beta0_hat_sd,
    beta0_hat_u  = beta0_hat + z * beta0_hat_sd
  )

# ------------------------- define scenarios ---------------------------------
stopifnot(!any(is.na(df_est$beta0_hat)))

# # counterfactual betas (schools as observed; ALS constant throughout) 
k_ALS <- c(ALS2 = unname(m_als["ALS2"]),
           ALS4 = unname(m_als["ALS4"]),
           NoALS = unname(m_als["NoALS"]))

cut_school <- as.Date("2022-03-14")

m_closed <- unname(m_school["closed"])
m_open   <- unname(m_school["open"])

# # find first date with No-ALS in your policy schedule and enforce across when noals starts
noals_start <- df_est %>%
  dplyr::filter(ALS == "NoALS") %>%
  summarise(start = min(date)+1) %>%
  pull(start)
# 
s_school <- plogis(as.numeric(df_est$date - cut_school) / tau_days) # logistic step: ~0 before cut, ~1 after cut
m_s2_school <- m_closed * (1 - s_school) + m_open * s_school # smoothly blend from closed -> open

df_est <- df_est %>%
  dplyr::mutate(
    m_scen_ALS2  = if_else(date < noals_start + 1, k_ALS["ALS2"] * m_school_hat, m_school_hat),
    m_scen_ALS4  = if_else(date < noals_start,k_ALS["ALS4"] * m_school_hat, m_school_hat),
    m_scen_NoALS = k_ALS["NoALS"] * m_school_hat,

    beta_cf_ALS2  = beta0_hat * m_scen_ALS2,
    beta_cf_ALS4  = beta0_hat * m_scen_ALS4,
    beta_cf_NoALS = beta0_hat * m_scen_NoALS,
    
    # new betas for schools open/closed (counterfactuals) 
    beta_cf_Obs_ALS_Sch = beta_hat,
    beta_cf_Obs_ALS_Sch_Clo_March  = beta0_hat * m_ALS_hat * m_s2_school, # apply observed ALS
    beta_cf_Obs_ALS_Sch_Clos_Tru = beta0_hat * m_ALS_hat * m_closed,
    beta_cf_Obs_ALS_Sch_Open_Tru = beta0_hat * m_ALS_hat * m_open,
    beta_cf_NoALS_Sch_Clos_Tru = beta0_hat * m_closed
  )

print(df_est)

# ------------------ school counterfactual multipliers -----------------------

# observed ALS is already m_ALS_hat
# school counterfactuals:
df_est <- df_est %>%
  mutate(
    m_sch_obs = m_school_hat,
    m_sch_close_to_mar14 = m_s2_school,
    m_sch_closed_all = m_closed,
    m_sch_open_all = m_open,
    
    # combined multipliers (ALS observed unless explicitly removed)
    m_comb_obs                = m_ALS_hat * m_sch_obs,
    m_comb_close_to_mar14      = m_ALS_hat * m_sch_close_to_mar14,
    m_comb_closed_all          = m_ALS_hat * m_sch_closed_all,
    m_comb_open_all            = m_ALS_hat * m_sch_open_all,
    m_comb_noals_closed_all    = 1.0      * m_sch_closed_all
  )


als_change_dates <- df_est |>
  mutate(chg = ALS != lag(ALS)) |>
  dplyr::filter(!is.na(chg) & chg) |>
  distinct(date)

school_change_dates <- df_est |>
  mutate(chg = school != lag(school)) |>
  dplyr::filter(!is.na(chg) & chg) |>
  distinct(date)

# methods figure: (a) multipliers, (b) β0 vs β̂
mult_long <- df_est |>
  select(date, `ALS`=m_ALS_hat, `School`=m_school_hat, `Combined`=m_obs_hat) |>
  pivot_longer(-date, names_to = "series", values_to = "value")

pB_mult <- ggplot(mult_long, aes(date, value, color = series)) +
  geom_line(linewidth = 0.8, alpha = 1.0) +
  scale_color_manual(values = c("ALS"="#D01125","School"="#40B0A6","Combined"="#E69F00")) +
  labs(title="Policy effect multipliers", y="Multiplier", x=NULL, color=NULL) +
  scale_x_date(date_breaks="2 weeks", date_labels="%b %d", expand = c(0,0)) +
  geom_vline(data = als_change_dates,    aes(xintercept = date), linetype="dashed",  linewidth=1) +
  geom_vline(data = school_change_dates, aes(xintercept = date), linetype="dotdash", linewidth=1, color="brown")+
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.85),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black", hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(color = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank(),
    legend.position = c(0.85, 0.5)
  )

p_beta_overlay <- ggplot(df_est, aes(date)) +
  geom_line(aes(y = beta_hat,  color = "Estimated (Implemented)"),      linewidth = 1.5) +
  geom_line(aes(y = beta0_hat, color = "No Policy"),   linewidth = 1) +
  geom_line(aes(y = beta_cf_ALS2, color = "Beta2"), linewidth = 1) +
  geom_line(aes(y = beta_cf_ALS4, color = "Beta4"),   linewidth = 1) +
  labs(title = "Estimated transmission rates vs No policy transmission",x = NULL, y = "Transmission rate") +
  scale_color_manual(values = c("Estimated (Implemented)" = "red","No Policy" = "#0072B2", "Beta2" = "#E69F00", "Beta4" = "#009E73"), name = NULL) +
  scale_fill_manual(values = c("Estimated (Implemented)" = "red", "No Policy" = "blue"),guide = "none") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0, 0)) +
  geom_vline(data = als_change_dates,    aes(xintercept = date), linetype="dashed",  linewidth=1) +
  geom_vline(data = school_change_dates, aes(xintercept = date), linetype="dotdash", linewidth=1, color="brown")+
  theme_clean() +
  theme(
    axis.text.x  = element_text(size = 12, angle = 0, hjust = 0.85),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title   = element_text(size = 12, color = "black", hjust = 0.5, face = "plain"),
    strip.text   = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 11),
    legend.background = element_rect(color = NA),
    legend.margin     = margin(0, 0, 0, 0),
    plot.background   = element_blank(),
    legend.position   = c(0.9, 0.25),
    plot.margin  = margin(0, 6, -10, 6)
  )

## ALS
als_df <- tibble(
  xmin      = as.Date(c("2021-12-15","2021-12-25","2022-01-05", "2022-02-11", "2022-03-20")),
  xmax      = as.Date(c("2021-12-25","2022-01-05","2022-02-11", "2022-03-20", "2022-05-22")),
  sch_label = c("ALS-2","ALS-3","ALS-4","ALS-3","No-ALS"),
  seg_col   = c("navy","navy","navy","navy","navy"),
  lab_col   = c("navy","navy","navy","navy","navy")
) %>%
  mutate(x_label = as.Date((as.numeric(xmin) + as.numeric(xmax))/2, origin = "1970-01-01"))

als_imp <- ggplot() +
  geom_segment(
    data = als_df,
    aes(x = xmin, xend = xmax, y = 1, yend = 1),
    colour = c("navy","navy","navy","navy","navy"),
    linewidth = 1.5,
    arrow = arrow(angle = 90, ends = "both", length = unit(0.5, "cm"))
  ) +
  geom_text(
    data = als_df,
    aes(x = x_label, y = 0.50, label = sch_label),
    colour = c("navy","navy","navy","navy","navy"),
    size = 4
  ) +
  scale_x_date(
    limits = c(as.Date("2021-12-15"), as.Date("2022-05-22")),
    expand = c(0, 0),
    breaks = NULL
  ) +
  ylim(-1, 1) +
  theme_void() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(0, 0, 0, 0)
  )

sch_op_df <- tibble(
  xmin      = as.Date(c("2021-12-15","2021-12-20","2022-01-29")),
  xmax      = as.Date(c("2021-12-19","2022-01-28","2022-05-22")),
  sch_label = c("K-12\nSch.\nOpen","K-12 Schools Closed","K-12 Schools Open"),
  seg_col   = c("navy","red","navy"),
  lab_col   = c("navy","red","navy")
) %>%
  mutate(x_label = as.Date((as.numeric(xmin) + as.numeric(xmax))/2, origin = "1970-01-01"))

sch_operation <- ggplot() +
  geom_segment(
    data = sch_op_df,
    aes(x = xmin, xend = xmax, y = 1, yend = 1),
    colour = c("navy","brown","navy"),
    linewidth = 1.5,
    arrow = arrow(angle = 90, ends = "both", length = unit(0.5, "cm"))
  ) +
  geom_text(
    data = sch_op_df,
    aes(x = x_label, y = 0.35, label = sch_label),
    colour = c("navy","red","navy"),
    size = 4
  ) +
  scale_x_date(
    limits = c(as.Date("2021-12-15"), as.Date("2022-05-22")),
    expand = c(0, 0),
    breaks = NULL
  ) +
  ylim(-1, 1) +
  theme_void() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(0, 0, 0, 0)
  )

fig_methods <- pB_mult / p_beta_overlay / als_imp /sch_operation + plot_layout(heights = c(1, 1, 0.25, 0.25))

print(fig_methods)

# time-varying multiplier heatmap
heatdf <- df_est %>%
  select(date, m_ALS_hat, m_school_hat, m_obs_hat) %>%
  pivot_longer(-date, names_to = "component", values_to = "multiplier") %>%
  mutate(component = recode(component, m_ALS_hat = "ALS",m_school_hat = "School", m_obs_hat = "Combined"))

p_heat_time <- ggplot(heatdf, aes(date, component, fill = multiplier)) +
  geom_tile(height = 0.8) +
  scale_fill_viridis_c(name = "", option = "C") +
  labs(title = "Policy effect multipliers over time", x = NULL, y = NULL) +
  geom_vline(data=als_change_dates, aes(xintercept = as.numeric(date)), linetype = "dashed", linewidth = 1) +
  geom_vline(data=school_change_dates, aes(xintercept = as.numeric(date)), linetype = "dotdash", linewidth = 1, color = "brown") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0,0)) +
  theme_clean()+
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.85),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black", hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 12, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(color = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank(),
    legend.position = "right"
  )

print(p_heat_time)

png("../figures/multipliers.png", width = 4000, height = 2000, res = 300, bg = "white", type = "cairo")
fig_methods
dev.off()

png("../figures/multiplier_heatmap.png", width = 3000, height = 2500, res = 300, bg = "white", type = "cairo")
p_heat_time
dev.off()



# ------------------ beta overlays: school counterfactuals -------------------

beta_school_long <- df_est %>%
  select(
    date,
    beta_hat,
    beta0_hat,
    beta_cf_Obs_ALS_Sch,
    beta_cf_Obs_ALS_Sch_Clo_March,
    beta_cf_Obs_ALS_Sch_Clos_Tru,
    beta_cf_Obs_ALS_Sch_Open_Tru,
    beta_cf_NoALS_Sch_Clos_Tru
  ) %>%
  pivot_longer(-date, names_to = "series", values_to = "beta") %>%
  mutate(
    series = recode(
      series,
      beta_hat = "Estimated (implemented)",
      beta0_hat = "No policy baseline",
      beta_cf_Obs_ALS_Sch = "Observed ALS + observed school",
      beta_cf_Obs_ALS_Sch_Clo_March = "Observed ALS + school closed → Mar 14",
      beta_cf_Obs_ALS_Sch_Clos_Tru = "Observed ALS + school closed throughout",
      beta_cf_Obs_ALS_Sch_Open_Tru = "Observed ALS + school open throughout",
      beta_cf_NoALS_Sch_Clos_Tru = "No ALS + school closed throughout"
    ),
    series = factor(series, levels = c(
      "Estimated (implemented)",
      "No policy baseline",
      "Observed ALS + observed school",
      "Observed ALS + school closed → Mar 14",
      "Observed ALS + school closed throughout",
      "Observed ALS + school open throughout",
      "No ALS + school closed throughout"
    ))
  )

p_beta_school <- ggplot(beta_school_long, aes(date, beta, color = series)) +
  geom_line(linewidth = 1.05) +
  labs(
    title = "Transmission rate under school closure counterfactuals",
    x = NULL,
    y = "Transmission rate",
    color = NULL
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0, 0)) +
  geom_vline(data = als_change_dates, aes(xintercept = date), linetype = "dashed", linewidth = 1) +
  geom_vline(data = school_change_dates, aes(xintercept = date), linetype = "dotdash", linewidth = 1, color = "brown") +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, color = "black", hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(color = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank(),
    legend.position = c(0.8,0.85)
  )

print(p_beta_school)
png("../figures/beta_school_counterfactuals.png", width = 4000, height = 2000, res = 300, bg = "white", type = "cairo")
print(p_beta_school)
dev.off()

# ------------------ multiplier lines: school counterfactuals ----------------

mult_school_long <- df_est %>%
  select(
    date,
    `Observed ALS` = m_ALS_hat,
    `School (observed)` = m_sch_obs,
    `School (closed → Mar 14)` = m_sch_close_to_mar14,
    `School (closed throughout)` = m_sch_closed_all,
    `School (open throughout)` = m_sch_open_all
  ) %>%
  pivot_longer(-date, names_to = "series", values_to = "value")

p_mult_school <- ggplot(mult_school_long, aes(date, value, color = series)) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "School closure counterfactual multipliers (ALS held observed)",
    x = NULL,
    y = "Multiplier",
    color = NULL
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0, 0)) +
  geom_vline(data = als_change_dates, aes(xintercept = date), linetype = "dashed", linewidth = 1) +
  geom_vline(data = school_change_dates, aes(xintercept = date), linetype = "dotdash", linewidth = 1, color = "brown") +
  theme_clean()+
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, color = "black", hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(color = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank(),
    legend.position = c(0.8,0.6)
  )


print(p_mult_school)

png("../figures/mult_school_counterfactuals.png", width = 4000, height = 2000, res = 300, bg = "white", type = "cairo")
print(p_mult_school)
dev.off()


# ------------------ heatmap: school counterfactuals -------------------------

heat_school <- df_est %>%
  transmute(
    date,
    beta0_hat,
    `Observed ALS + observed school`            = beta_cf_Obs_ALS_Sch,
    `Observed ALS + school closed → Mar 14`     = beta_cf_Obs_ALS_Sch_Clo_March,
    `Observed ALS + school closed throughout`   = beta_cf_Obs_ALS_Sch_Clos_Tru,
    `Observed ALS + school open throughout`     = beta_cf_Obs_ALS_Sch_Open_Tru,
    `No ALS + school closed throughout`         = beta_cf_NoALS_Sch_Clos_Tru
  ) %>%
  pivot_longer(-c(date, beta0_hat), names_to = "scenario", values_to = "beta") %>%
  mutate(
    beta_rel = beta / beta0_hat,
    scenario = factor(scenario, levels = c(
      "Observed ALS + observed school",
      "Observed ALS + school closed → Mar 14",
      "Observed ALS + school closed throughout",
      "Observed ALS + school open throughout",
      "No ALS + school closed throughout"
    ))
  )

p_heat_school <- ggplot(heat_school, aes(date, scenario, fill = beta_rel)) +
  geom_tile(height = 0.9) +
  scale_fill_viridis_c(option = "C", name = "Relative\ntransmission") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0, 0)) +
  labs(
    title = "Relative transmission under school counterfactual scenarios",
    x = NULL,
    y = NULL
  ) +
  theme_clean()+
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, color = "black", hjust = 0.5, face = "plain"),
    strip.text = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(color = NA),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank(),
    legend.position = "top"
  )


print(p_heat_school)

png("../figures/heat_school_counterfactuals.png", width = 6000, height = 4200, res = 300, bg = "white", type = "cairo")
print(p_heat_school)
dev.off()

fig_school_methods <- p_mult_school / p_beta_school + plot_layout(heights = c(1, 1))

png("../figures/school_counterfactuals_methods.png", width = 6000, height = 4200, res = 300, bg = "white", type = "cairo")
print(fig_school_methods)
dev.off()



rdsSave(df_est)

