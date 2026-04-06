suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
  library(scales)
  library(ggthemes)
  library(shellpipes)
  library(grid)
  library(patchwork)
})

loadEnvironments()

# ============================================================
# Constants
# ============================================================
NL_POP           <- 510550
CA_POP           <- 36991981
eligible_pop <- 426198

init_date      <- as.Date("2021-07-24")   
booster_start  <- as.Date("2021-11-08")
backlog_end    <- as.Date("2022-01-11")
cutoff_date    <- as.Date("2022-02-07")
analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")
study_end      <- as.Date("2022-05-22")   
earlier_vax_date <- as.Date("2021-10-04")
booster_offset <- 150

target_date       <- as.Date("2022-01-11")
target_cumulative <- 106730
shape_power       <- 1.8

day0 <- as.Date("2021-12-17") - lubridate::days(booster_offset)

clip01 <- function(x) pmin(pmax(x, 0), 1)

# ============================================================
# Load raw data
# ============================================================
dat <- csvRead()

# load and plot canada booster porportion scaled to NL
can_lower_init <- ca_weekly_pts |>
  select(week_end, nl_scaled_prop) 

# bind this the other scenarios after they are created in line 
canada_weekly_vax_nl <- can_weekly_nl_inputs |> select(week_end, weekly_rate)|> 
  mutate(days = as.numeric(week_end - day0), scenario = "National Rate")|> 
  select(scenario, week_end, weekly_rate)


# ============================================================
# Raw NL booster cumulative counts
# ============================================================
df_raw <- dat |>
  filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional), 0, numtotal_2nd_additional),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  select(week_end, numtotal_additional) |>
  rename(value = numtotal_additional) |>
  filter(week_end >= earlier_vax_date, week_end <= as.Date("2022-06-19")) |>
  arrange(week_end)

# ============================================================
# Correct cumulative booster series
# ============================================================

first_obs <- df_raw |>
  filter(!is.na(value)) |>
  arrange(week_end) |>
  slice(1)

first_obs_date  <- first_obs$week_end
first_obs_value <- first_obs$value

dates_a <- seq.Date(booster_start, target_date, by = "week")
t_idx   <- seq(0, 1, length.out = length(dates_a))

seg_a <- tibble(
  week_end = dates_a,
  corrected_value = target_cumulative * (t_idx ^ shape_power),
  point_type = "interpolated"
)

if (first_obs_date > target_date) {
  dates_b <- seq.Date(target_date, first_obs_date, by = "week")
  seg_b <- tibble(
    week_end = dates_b,
    corrected_value = seq(target_cumulative, first_obs_value, length.out = length(dates_b)),
    point_type = "interpolated"
  ) |>
    filter(week_end > target_date)
} else {
  seg_b <- tibble(
    week_end = as.Date(character()),
    corrected_value = numeric(),
    point_type = character()
  )
}

df_imputed <- bind_rows(seg_a, seg_b) |>
  filter(week_end < first_obs_date) |>
  distinct(week_end, .keep_all = TRUE) |>
  arrange(week_end)

df_true <- df_raw |>
  filter(week_end >= first_obs_date) |>
  transmute(
    week_end,
    corrected_value = value,
    point_type = "observed"
  )

# ============================================================
# Redistribute post-flat jumps backward
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
          df$point_type[idx[idx < j]] <- "interpolated"
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
# Process observed booster data
# ============================================================
nl_actual <- bind_rows(df_imputed, df_true) |>
  arrange(week_end) |>
  redistribute_postflat_jumps(backlog_end = backlog_end) |>
  mutate(
    cum_prop        = corrected_value / eligible_pop,
    daily_inc_star  = (corrected_value - lag(corrected_value, default = 0)) / 7,
    weekly_inc_star = corrected_value - lag(corrected_value, default = 0)
  )

# ============================================================
# Build baseline daily cumulative rollout
# ============================================================
daily_dates <- seq(booster_start, study_end, by = "1 day")

fit_func <- splinefun(
  x = nl_actual$week_end,
  y = nl_actual$corrected_value,
  method = "hyman"
)

last_obs_date <- max(nl_actual$week_end)
last_obs_cum  <- nl_actual |>
  filter(week_end == last_obs_date) |>
  pull(corrected_value)

daily_base <- tibble(date = daily_dates) |>
  mutate(
    cum_value_raw = fit_func(date),
    cum_value_raw = pmax(cum_value_raw, 0),
    cum_value = if_else(date <= last_obs_date, cum_value_raw, last_obs_cum),
    cum_value = pmin(cum_value, eligible_pop),
    cum_value = cummax(cum_value),
    daily_inc = c(0, diff(cum_value)),
    daily_inc = pmax(daily_inc, 0),
    week_end  = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(week_end) |>
  mutate(
    weekly_rate = sum(daily_inc, na.rm = TRUE)
  ) |>
  ungroup()

# ============================================================
# Time-warp cumulative rollout
# speed > 1  => faster
# speed < 1  => slower
# shift_days < 0 => earlier
# ============================================================
make_rollout_scenario <- function(df, speed = 1, shift_days = 0, label) {
  base_dates <- df$date
  base_num   <- as.numeric(base_dates)
  base_cum   <- df$cum_value
  
  cum_fun <- approxfun(
    x = base_num,
    y = base_cum,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  t0 <- min(base_num)
  scen_num <- as.numeric(base_dates)
  
  mapped_t <- t0 + speed * ((scen_num - shift_days) - t0)
  
  scen_cum <- cum_fun(mapped_t)
  scen_cum <- pmax(scen_cum, 0)
  scen_cum <- pmin(scen_cum, eligible_pop)
  scen_cum <- cummax(scen_cum)
  
  tibble(
    date    = base_dates,
    cum_val = scen_cum
  ) |>
    mutate(
      daily_rate = c(0, diff(cum_val)),
      daily_rate = pmax(daily_rate, 0),
      cum_prop   = clip01(cum_val / eligible_pop),
      week_end   = ceiling_date(date, unit = "week", week_start = 7) - days(1),
      scenario   = label
    ) |>
    group_by(scenario, week_end) |>
    mutate(
      weekly_rate = sum(daily_rate, na.rm = TRUE)
    ) |>
    ungroup()
}

# ============================================================
# Vaccination scenarios
# ============================================================
daily_scenarios_df <- bind_rows(
  make_rollout_scenario(daily_base, speed = 1.00,   shift_days =   0, label = "Actual Boosters"),
  make_rollout_scenario(daily_base, speed = 1.00,   shift_days = -35, label = "Early Boosters (-35 Days)"),
  make_rollout_scenario(daily_base, speed = 1.50,   shift_days =   0, label = "Faster Rate"),
  make_rollout_scenario(daily_base, speed = 1/1.5,  shift_days =   0, label = "Slower Rate")
) |>
  mutate(
    days_since_init = as.integer(date - init_date)
  ) |>
  filter(days_since_init >= 0)

weekly_scenarios_df <- daily_scenarios_df |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = max(weekly_rate, na.rm = TRUE),
    cum_prop    = max(cum_prop, na.rm = TRUE),
    cum_val     = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  ) 

# ============================================================
# Build model-ready booster inputs (daily)
# one changepoint per day, anchored to init_date
# ============================================================
booster_inputs_df <- daily_scenarios_df |>
  transmute(
    scenario,
    date,
    day = as.integer(date - init_date),
    booster_vac_value = as.numeric(daily_rate)
  ) |>
  arrange(scenario, day)

# ============================================================
# Plotting
# ============================================================
theme_pub <- theme_clean() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, color = "black"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 24, face = "plain", hjust = 0.5),
    plot.background = element_blank()
  )

actual_stars <- nl_actual |>
  filter(week_end <= cutoff_date)

p1 <- ggplot(daily_scenarios_df, aes(x = date, y = daily_rate, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_point(
    data = actual_stars,
    aes(x = week_end, y = daily_inc_star),
    color = "black",
    shape = 8,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_x_date(date_breaks = "month", date_labels = "%Y-%m", expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Daily Booster Delivery Rate Comparison", y = "Daily booster doses", x = NULL) +
  theme_pub


vax_cols <- c(
  "Slower Rate"               = "magenta",
  "Early Boosters (-35 Days)" = "blue",
  "Actual Boosters"           = "black",
  "Faster Rate"               = "orange"
)

p2 <- ggplot(daily_scenarios_df |> filter(date <= as.Date("2022-05-22")),aes(x = date, y = cum_prop, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(
    data = weekly_scenarios_df |> filter(week_end <= as.Date("2022-05-22"),!(scenario == "Actual Boosters" & week_end <= cutoff_date)),
    aes(x = week_end, y = cum_prop, color = scenario),
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = actual_stars,
    aes(x = week_end, y = cum_prop),
    color = "black",
    shape = 8,
    size = 8,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = can_lower_init[!is.na(can_lower_init$nl_scaled_prop), ],
    aes(x = week_end, y = nl_scaled_prop, color = "National Rate"),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  geom_point(
    data = can_lower_init,
    aes(x = week_end, y = nl_scaled_prop, color = "National Rate"),
    inherit.aes = FALSE,
    size = 5
  ) +
  geom_vline(
    xintercept = booster_start,
    color = "maroon",
    linetype = 2,
    linewidth = 1.5
  ) +
  scale_colour_manual(values = c(vax_cols,"National Rate" =  "springgreen3")) +
  scale_x_date(
    limits = c(booster_start, as.Date("2022-05-22")),
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0.05,0)
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 0.9)
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.9)) +
  labs(
    title = "Proportions of booster vaccination coverage in NL",
    y = "Proportion of NL eligible population",
    x = NULL
  ) +
  theme_pub+
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(tag = "B") +
  theme(
      plot.tag = element_text(size = 35, face = "bold", colour = "black"),
      plot.tag.position = c(0.01, 0.98)
    )

p3 <- ggplot(weekly_scenarios_df, aes(x = week_end, y = weekly_rate, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_point(
    data = actual_stars,
    aes(x = week_end, y = weekly_inc_star),
    color = "black",
    shape = 8,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Weekly Booster Delivery Rate Comparison", y = "Weekly booster doses", x = NULL) +
  theme_pub

# png("../figures/daily_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
# print(p1)
# dev.off()

png("../figures/cum_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
print(p2)
dev.off()

# png("../figures/weekly_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
# print(p3)
# dev.off()

# ============================================================
# Save outputs for downstream simulation
# ============================================================

# row bind canadaina national average vax to Nl to vaccination scenarios
weekly_scenarios_df <- weekly_scenarios_df|> 
  select(scenario, week_end, weekly_rate)

weekly_scenarios_df <- rbind(weekly_scenarios_df,canada_weekly_vax_nl)

saveVars(weekly_scenarios_df,daily_scenarios_df, booster_inputs_df)

