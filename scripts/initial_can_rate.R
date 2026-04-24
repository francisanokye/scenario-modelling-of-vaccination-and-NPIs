
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(ggthemes)
  library(shellpipes)
  library(grid)
  library(patchwork)
})

loadEnvironments()

# ============================================================
# CONSTANTS
# ============================================================
NL_POP      <- 510550
CA_POP      <- 36991981
NL_ELIG_POP <- 426198

init_date      <- as.Date("2021-07-24")
booster_start  <- as.Date("2021-11-08")
backlog_end    <- as.Date("2022-01-11")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
study_end      <- as.Date("2022-05-22")

clip01 <- function(x) pmin(pmax(x, 0), 1)

# ============================================================
# LOAD RAW DATA
# ============================================================
dat <- csvRead()

# ============================================================
# CANADIAN BOOSTER COUNTS
# ============================================================
ca_booster_raw <- dat |>
  filter(grepl("^canada$", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional), 0, numtotal_2nd_additional),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  select(week_end, numtotal_additional) |>
  rename(value = numtotal_additional) |>
  filter(week_end >= booster_start, week_end <= as.Date("2022-06-19")) |>
  arrange(week_end)

# ============================================================
# CANADIAN FULLY VACCINATED COUNTS
# ============================================================
ca_fully <- dat |>
  filter(grepl("^canada$", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  transmute(
    week_end,
    fully_count = as.numeric(numtotal_fully)
  ) |>
  arrange(week_end)

# ============================================================
# REDISTRIBUTE POST-FLAT JUMPS
# ============================================================
redistribute_postflat_jumps <- function(df, backlog_end) {
  df <- df |>
    arrange(week_end)
  
  n <- nrow(df)
  i <- 2
  
  while (i <= n) {
    if (df$week_end[i] > backlog_end && df$corrected_value[i] == df$corrected_value[i - 1]) {
      plateau_start <- i - 1
      j <- i
      
      while (j <= n && df$corrected_value[j] == df$corrected_value[plateau_start]) {
        j <- j + 1
      }
      
      if (j <= n) {
        start_val <- df$corrected_value[plateau_start]
        end_val   <- df$corrected_value[j]
        
        if (end_val > start_val) {
          idx <- plateau_start:j
          df$corrected_value[idx] <- seq(start_val, end_val, length.out = length(idx))
        }
        
        i <- j
      } else {
        i <- n + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  df
}

# ============================================================
# BUILD CORRECTED CANADIAN WEEKLY SERIES
# IMPORTANT:
# - use actual reporting dates
# - add one explicit zero row at booster_start
# ============================================================
ca_booster_corr <- ca_booster_raw |>
  transmute(
    week_end,
    corrected_value = value
  ) |>
  bind_rows(
    tibble(
      week_end = booster_start,
      corrected_value = 0
    )
  ) |>
  distinct(week_end, .keep_all = TRUE) |>
  arrange(week_end) |>
  redistribute_postflat_jumps(backlog_end = backlog_end)

# ============================================================
# JOIN TO FULLY VACCINATED COUNTS
# If booster_start is not present in ca_fully, interpolate fully_count
# ============================================================
fully_fun <- approxfun(
  x = as.numeric(ca_fully$week_end),
  y = ca_fully$fully_count,
  method = "linear",
  rule = 2,
  ties = "ordered"
)

ca_weekly <- ca_booster_corr |>
  mutate(
    fully_count = fully_fun(as.numeric(week_end))
  ) |>
  mutate(
    booster_uptake_prop = clip01(corrected_value / pmax(fully_count, 1)),
    fully_prop_totalpop = clip01(fully_count / CA_POP)
  ) |>
  arrange(week_end)

# ============================================================
# MAP CANADIAN BOOSTER UPTAKE PROPORTION TO NL ELIGIBLE POP
# ============================================================
ca_weekly_nl_scaled <- ca_weekly |>
  mutate(
    nl_scaled_cum  = booster_uptake_prop * NL_ELIG_POP,
    nl_scaled_prop = clip01(nl_scaled_cum / NL_ELIG_POP)
  )

# ============================================================
# DAILY INTERPOLATION
# Use approxfun, not splinefun, for robustness
# ============================================================
daily_dates <- seq(booster_start, study_end, by = "1 day")

prop_fun <- approxfun(
  x = as.numeric(ca_weekly_nl_scaled$week_end),
  y = ca_weekly_nl_scaled$booster_uptake_prop,
  method = "linear",
  rule = 2,
  ties = "ordered"
)

ca_daily_nl_scaled <- tibble(date = daily_dates) |>
  mutate(
    uptake_prop = clip01(prop_fun(as.numeric(date))),
    nl_scaled_cum = uptake_prop * NL_ELIG_POP,
    nl_scaled_cum = cummax(nl_scaled_cum),
    nl_scaled_cum = pmin(nl_scaled_cum, NL_ELIG_POP),
    daily_rate = c(0, diff(nl_scaled_cum)),
    daily_rate = pmax(daily_rate, 0),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1),
    days_since_init = as.integer(date - init_date)
  ) |>
  group_by(week_end) |>
  mutate(
    weekly_rate = sum(daily_rate, na.rm = TRUE),
    weekly_cum = max(nl_scaled_cum, na.rm = TRUE),
    weekly_prop = max(nl_scaled_cum / NL_ELIG_POP, na.rm = TRUE)
  ) |>
  ungroup()

# ============================================================
# INITIAL VALUES AT ANALYSIS START
# ============================================================
fully_prop_ts <- approxfun(
  x = as.numeric(ca_weekly$week_end),
  y = ca_weekly$fully_prop_totalpop,
  method = "linear",
  rule = 2,
  ties = "ordered"
)(as.numeric(init_date))

booster_uptake_prop_ts <- approxfun(
  x = as.numeric(ca_weekly$week_end),
  y = ca_weekly$booster_uptake_prop,
  method = "linear",
  rule = 2,
  ties = "ordered"
)(as.numeric(booster_start))

nl_full_vax_total_ts <- fully_prop_ts * NL_POP
nl_boosted_ts <- nl_full_vax_total_ts * booster_uptake_prop_ts
nl_two_dose_no_booster_ts <- pmax(nl_full_vax_total_ts - nl_boosted_ts, 0)

can_init_vax_for_nl <- tibble(
  date = init_date,
  ca_fully_prop_totalpop = fully_prop_ts,
  ca_booster_uptake_prop = booster_uptake_prop_ts,
  nl_full_vax_total = nl_full_vax_total_ts,
  nl_two_dose_no_booster = nl_two_dose_no_booster_ts,
  nl_boosted = nl_boosted_ts
)

# ============================================================
# MODEL INPUTS
# ============================================================
can_booster_inputs_df <- ca_daily_nl_scaled |>
  transmute(
    scenario = "Canada proportions mapped to NL",
    date,
    day = days_since_init,
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  filter(day >= 0) |>
  arrange(day)

can_weekly_nl_inputs <- ca_daily_nl_scaled |>
  group_by(week_end) |>
  summarise(
    weekly_rate = max(weekly_rate, na.rm = TRUE),
    nl_scaled_cum = max(nl_scaled_cum, na.rm = TRUE),
    nl_scaled_prop = max(weekly_prop, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# PLOTS
# ============================================================
theme_pub <- theme_clean() +
  theme(
    axis.text.x = element_text(size = 14, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18, color = "black"),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.background = element_blank()
  )

ca_weekly_pts <- ca_weekly_nl_scaled |>
  transmute(
    week_end,
    nl_scaled_prop,
    nl_scaled_cum
  )

p_prop <- ggplot(ca_daily_nl_scaled, aes(x = date, y = nl_scaled_cum / NL_ELIG_POP)) +
  geom_line(linewidth = 1.1, colour = "blue") +
  geom_point(
    data = ca_weekly_pts,
    aes(x = week_end, y = nl_scaled_prop),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.8
  ) +
  scale_x_date(
    limits = c(booster_start, analysis_end),
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Canadian booster uptake proportion mapped to NL",
    y = "Proportion of NL eligible population",
    x = NULL
  ) +
  theme_pub

p_cum <- ggplot(ca_daily_nl_scaled, aes(x = date, y = nl_scaled_cum)) +
  geom_line(linewidth = 1.1, colour = "darkgreen") +
  geom_point(
    data = ca_weekly_pts,
    aes(x = week_end, y = nl_scaled_cum),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.8
  ) +
  scale_x_date(
    limits = c(booster_start, analysis_end),
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Canadian booster cumulative doses scaled to NL",
    y = "Cumulative booster doses",
    x = NULL
  ) +
  theme_pub

p_all <- p_prop / p_cum
print(p_all)

# png("../figures/canada_scaled_to_nl_boosters_corrected.png",width = 5000, height = 5200, res = 300, bg = "white")
# print(p_all)
# dev.off()

# save
saveVars(can_init_vax_for_nl, can_weekly_nl_inputs, ca_weekly_pts, ca_daily_nl_scaled)
# saveVars(canada_init_vax_nl, canada_weekly_vax_nl, ca_weekly_pts, canada_daily_vax_nl)
