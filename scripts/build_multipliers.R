
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2); library(ggh4x)
  library(scales); library(tibble); library(tidyverse)
  library(broom.mixed); library(zoo); library(patchwork)
})

set.seed(2025)
try(loadEnvironments(), silent = TRUE)

date0    <- as.Date("2021-12-15")
date_end <- as.Date("2022-05-22")
offset0  <- 150
upper_plot_time <- 308
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
  mutate(date = seq(start_date, by = "day", length.out = length(value)),
         beta_hat = value) |>
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
#R_ALS4 <- 1.91; R_ALS3 <- 2.03; R_NoALS <- 2.26; R_open <- 2.17; R_closed <- 1.49; R_ALS2 <- sqrt(R_ALS3 * R_NoALS)
R_ALS4 <- 2.23; R_ALS3 <- 2.74; R_NoALS <- 3.0; R_open <- 2.71; R_closed <- 1.98; R_ALS2 <- sqrt(R_ALS3 * R_NoALS)
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

# counterfactual betas (schools as observed; ALS constant throughout)
k_ALS <- c(
  ALS4  = unname(m_als["ALS4"]),
  NoALS = unname(m_als["NoALS"])
)

noals_start <- df_est %>%
  dplyr::filter(ALS == "NoALS") %>%
  summarise(start = min(date) + 1) %>%
  pull(start)

df_est <- df_est %>%
  mutate(
    # Create the weight for the transition (0 to 1)
    s_trans = plogis(as.numeric(date - noals_start) / tau_days),
    
    # Smoothly transition from ALS multiplier to 1.0
    m_scen_ALS4  = (k_ALS["ALS4"] * (1 - s_trans) + 1.0 * s_trans) * m_school_hat,
    m_scen_NoALS = k_ALS["NoALS"] * m_school_hat,
    
    # Now calculate your betas
    beta_cf_ALS4  = beta0_hat * m_scen_ALS4,
    beta_cf_NoALS = beta0_hat * m_scen_NoALS
  )


# ===================== ALS-4 DAYS ONLY COUNTERFACTUALS =====================
# ALS-4 starts on Dec 15, 2021 and lasts exactly d days, then NoALS.
# Built on full model grid so length(beta_ALS4_*days_only) == time_steps

full_df <- tibble(
  time = seq_len(time_steps),
  date = start_date + (seq_len(time_steps) - 1L)
)

df_est_full <- full_df %>%
  left_join(
    df_est %>% select(date, beta0_hat, m_school_hat),
    by = "date"
  ) %>%
  tidyr::fill(beta0_hat, m_school_hat, .direction = "downup")

als4_start_only <- as.Date("2021-12-15")
als4_only_days  <- c(0, 7L, 14L, 28L, 40L, 60L, 90L, 120L, 150L, 180L)

make_ALS4_days_only_then_NoALS <- function(date, als4_days,
                                           als4_start = as.Date("2021-12-15"),
                                           end_date   = max(date)) {
  als4_end_new <- als4_start + (as.integer(als4_days) - 1L)
  
  case_when(
    date >= als4_start & date <= als4_end_new ~ "ALS4",
    date >= (als4_end_new + 1L) & date <= end_date ~ "NoALS",
    TRUE ~ "ALS4"
  )
}

for (d in als4_only_days) {
  
  nm_ALS   <- paste0("m_ALS_hat_ALS4_", d, "days_only")
  nm_obs   <- paste0("m_obs_hat_ALS4_", d, "days_only")
  nm_beta  <- paste0("beta_ALS4_", d, "days_only")
  nm_sd    <- paste0("beta_ALS4_", d, "days_only_sd")
  nm_low   <- paste0("beta_ALS4_", d, "days_only_l")
  nm_high  <- paste0("beta_ALS4_", d, "days_only_u")
  
  ALS_tmp_full <- make_ALS4_days_only_then_NoALS(
    date       = df_est_full$date,
    als4_days  = d,
    als4_start = als4_start_only,
    end_date   = max(df_est_full$date)
  )
  
  df_est_full[[nm_ALS]] <- make_smooth_series(
    date = df_est_full$date,
    lev  = ALS_tmp_full,
    map  = m_als,
    tau  = tau_days
  )
  
  df_est_full[[nm_obs]]  <- df_est_full[[nm_ALS]] * df_est_full$m_school_hat
  df_est_full[[nm_beta]] <- df_est_full$beta0_hat * df_est_full[[nm_obs]]
  
  df_est_full[[nm_sd]]   <- sd(df_est_full[[nm_beta]], na.rm = TRUE)
  df_est_full[[nm_low]]  <- df_est_full[[nm_beta]] - z * df_est_full[[nm_sd]]
  df_est_full[[nm_high]] <- df_est_full[[nm_beta]] + z * df_est_full[[nm_sd]]
  
  stopifnot(length(df_est_full[[nm_beta]]) == time_steps)
}

# copy plotting-window values back into df_est using same names
keep_cols_days_only <- unlist(lapply(als4_only_days, function(d) {
  c(
    paste0("m_ALS_hat_ALS4_", d, "days_only"),
    paste0("m_obs_hat_ALS4_", d, "days_only"),
    paste0("beta_ALS4_", d, "days_only"),
    paste0("beta_ALS4_", d, "days_only_sd"),
    paste0("beta_ALS4_", d, "days_only_l"),
    paste0("beta_ALS4_", d, "days_only_u")
  )
}))

df_est <- df_est %>%
  left_join(
    df_est_full %>% select(date, all_of(keep_cols_days_only)),
    by = "date"
  )

# optional check
for (d in als4_only_days) {
  als4_end_new <- als4_start_only + (as.integer(d) - 1L)
  noals_start  <- als4_end_new + 1L
  message(
    "ALS-4 ", d, " days only: ",
    format(als4_start_only), " to ", format(als4_end_new),
    " | NoALS from ", format(noals_start)
  )
}


# ------------------ change dates for vertical lines -----------------------
als_change_dates <- df_est |>
  mutate(chg = ALS != lag(ALS)) |>
  dplyr::filter(!is.na(chg) & chg) |>
  distinct(date)

school_change_dates <- df_est |>
  mutate(chg = school != lag(school)) |>
  dplyr::filter(!is.na(chg) & chg) |>
  distinct(date)

# ------------------ plot: ALS-4 days-only counterfactuals -----------------
cols_new <- c(
  "Estimated (Implemented)" = "black",
  "No ALS"                  = "red",
  "ALS-4 Only"              = "blue",
  "ALS-4 0 days only"      = "#1b9e77",
  "ALS-4 7 days only"      = "#1b9e77",
  "ALS-4 14 days only"      = "#1b9e77",
  "ALS-4 28 days only"      = "#d95f02",
  "ALS-4 40 days only"      = "#7570b3",
  "ALS-4 60 days only"      = "#e7298a",
  "ALS-4 90 days only"      = "#66a61e",
  "ALS-4 120 days only"      = "brown",
  "ALS-4 150 days only"      = "magenta",
  "ALS-4 180 days only"      = "gray60"
)

p_beta_ALS4_days_only <- ggplot(df_est, aes(date)) +
  geom_ribbon(aes(ymin = beta_hat_l, ymax = beta_hat_u,
                  fill = "Estimated (Implemented)"), alpha = 0.18) +
  geom_ribbon(aes(ymin = beta0_hat_l, ymax = beta0_hat_u,
                  fill = "No ALS"), alpha = 0.18) +
  geom_ribbon(aes(ymin = beta_ALS4_14days_only_l, ymax = beta_ALS4_0days_only_u,
                  fill = "ALS-4 0 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_14days_only_l, ymax = beta_ALS4_7days_only_u,
                  fill = "ALS-4 7 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_14days_only_l, ymax = beta_ALS4_14days_only_u,
                  fill = "ALS-4 14 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_28days_only_l, ymax = beta_ALS4_28days_only_u,
                  fill = "ALS-4 28 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_40days_only_l, ymax = beta_ALS4_40days_only_u,
                  fill = "ALS-4 40 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_60days_only_l, ymax = beta_ALS4_60days_only_u,
                  fill = "ALS-4 60 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_90days_only_l, ymax = beta_ALS4_90days_only_u,
                  fill = "ALS-4 90 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_120days_only_l, ymax = beta_ALS4_120days_only_u,
                  fill = "ALS-4 120 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_150days_only_l, ymax = beta_ALS4_150days_only_u,
                  fill = "ALS-4 150 days only"), alpha = 0.10) +
  geom_ribbon(aes(ymin = beta_ALS4_180days_only_l, ymax = beta_ALS4_180days_only_u,
                  fill = "ALS-4 180 days only"), alpha = 0.10) +
  geom_line(aes(y = beta_hat, color = "Estimated (Implemented)"), linewidth = 2.0) +
  geom_line(aes(y = beta0_hat, color = "No ALS"), linewidth = 2.0) +
  geom_line(aes(y = beta_ALS4_14days_only, color = "ALS-4 0 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_14days_only, color = "ALS-4 7 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_14days_only, color = "ALS-4 14 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_28days_only, color = "ALS-4 28 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_40days_only, color = "ALS-4 40 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_60days_only, color = "ALS-4 60 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_90days_only, color = "ALS-4 90 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_120days_only, color = "ALS-4 120 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_150days_only, color = "ALS-4 150 days only"), linewidth = 1.5) +
  geom_line(aes(y = beta_ALS4_180days_only, color = "ALS-4 180 days only"), linewidth = 1.5) +
  
  scale_color_manual(values = cols_new, name = NULL) +
  scale_fill_manual(values = cols_new, guide = "none") +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d",
    expand = c(0, 0)
  ) +
  geom_vline(data = als_change_dates,
             aes(xintercept = date),
             linetype = "dashed",
             linewidth = 1.0) +
  geom_vline(data = school_change_dates,
             aes(xintercept = date),
             linetype = "dotdash",
             linewidth = 1.0,
             color = "brown") +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 22, hjust = 0.5, face = "plain"),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    strip.background = element_blank(),
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Estimated transmission rate and ALS-4 duration-only counterfactuals",
    x = NULL,
    y = "Transmission rate"
  )

p_beta_only <- ggplot(df_est, aes(date)) +
  geom_line(aes(y = beta_hat, color = "Estimated (Implemented)"), linewidth = 2.0) +
  geom_line(aes(y = beta0_hat, color = "No ALS"), linewidth = 2.0) +
  geom_line(aes(y = beta_cf_ALS4, color = "ALS-4 Only"), linewidth = 2.0) +
  scale_color_manual(values = cols_new, name = NULL) +
  scale_fill_manual(values = cols_new, guide = "none") +
  scale_x_date(date_breaks = "2 weeks",date_labels = "%b %d",expand = c(0, 0)) +
  geom_vline(data = als_change_dates,
             aes(xintercept = date),
             linetype = "dashed",
             linewidth = 1.0) +
  geom_vline(data = school_change_dates,
             aes(xintercept = date),
             linetype = "dotdash",
             linewidth = 1.0,
             color = "brown") +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 22, hjust = 0.5, face = "plain"),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    strip.background = element_blank(),
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Estimated transmission rate and ALS-4 duration-only counterfactuals",
    x = NULL,
    y = "Transmission rate"
  )

# png("../figures/p_beta_only.png",width = 4200, height = 2400, res = 300, bg = "white", type = "cairo")
# print(p_beta_only)
# dev.off()

selected_columns <- c("date","beta_hat", "beta0_hat", "beta_cf_ALS4",  
                      "beta_ALS4_0days_only","beta_ALS4_7days_only",
                      "beta_ALS4_14days_only","beta_ALS4_28days_only",
                      "beta_ALS4_40days_only","beta_ALS4_60days_only", 
                      "beta_ALS4_90days_only", "beta_ALS4_120days_only",
                      "beta_ALS4_150days_only", "beta_ALS4_180days_only")

df_est <- df_est |> 
  select(selected_columns)

print(df_est)

# rdsSave(df_est)
saveVars(df_est)
