
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(lubridate)
#   library(readr)
#   library(scales)
#   library(ggthemes)
#   library(shellpipes)
#   library(grid)
#   library(patchwork)
# })
# 
# loadEnvironments()
# 
# # ============================================================
# # Constants
# # ============================================================
# NL_POP <- 510550
# CA_POP <- 36991981
# BOOSTER_ELIG_POP <- 426198
# booster_offset <- 107
# booster_start <- as.Date("2021-11-08")
# backlog_end   <- as.Date("2022-01-11")
# plot_start    <- booster_start
# plot_end      <- as.Date("2022-05-22")
# init_date     <- as.Date("2021-07-24")
# 
# NL_POP <- 510550
# CUTOFF_DATE <- as.Date("2022-02-07")
# STUDY_END <- as.Date("2022-12-15")
# 
# target_cumulative <- 106730
# 
# clip01 <- function(x) pmin(pmax(x, 0), 1)
# 
# # ============================================================
# # Load raw data
# # ============================================================
# dat <- csvRead()
# 
# # ============================================================
# # NL booster raw data
# # ============================================================
# df_raw <- dat |>
#   filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
#   mutate(week_end = as.Date(week_end)) |>
#   select(week_end, contains("numtotal")) |>
#   mutate(
#     numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional), 0, numtotal_2nd_additional),
#     numtotal_additional = numtotal_1additional + numtotal_2nd_additional
#   ) |>
#   select(week_end, numtotal_additional, numtotal_fully) |>
#   pivot_longer(-week_end, names_to = "type", values_to = "value") |>
#   filter(type == "numtotal_additional") |>
#   filter(week_end >= as.Date("2021-11-06") & week_end <= as.Date("2022-06-19")) |>
#   arrange(week_end)
# 
# 
# # ============================================================
# # NL and Canada fully vaccinated counts
# # ============================================================
# nl_fully <- dat |>
#   filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
#   mutate(week_end = as.Date(week_end)) |>
#   transmute(
#     week_end,
#     fully_count = as.numeric(numtotal_fully)
#   )
# 
# ca_fully <- dat |>
#   filter(grepl("^canada$", prename, ignore.case = TRUE)) |>
#   mutate(week_end = as.Date(week_end)) |>
#   transmute(
#     week_end,
#     fully_count = as.numeric(numtotal_fully)
#   )
# 
# # ============================================================
# # Correct actual cumulative booster series
# # ============================================================
# eligible_pop <- BOOSTER_ELIG_POP
# target_date  <- as.Date("2022-01-11")
# target_value <- target_cumulative
# shape_power  <- 1.8
# 
# first_obs <- df_raw |>
#   filter(!is.na(value)) |>
#   arrange(week_end) |>
#   slice(1)
# 
# first_obs_date  <- first_obs$week_end
# first_obs_value <- first_obs$value
# 
# dates_a <- seq.Date(booster_start, target_date, by = "week")
# t_idx <- seq(0, 1, length.out = length(dates_a))
# 
# seg_a <- tibble(
#   week_end = dates_a,
#   corrected_value = target_value * (t_idx ^ shape_power),
#   point_type = "interpolated"
# )
# 
# if (first_obs_date > target_date) {
#   dates_b <- seq.Date(target_date, first_obs_date, by = "week")
#   seg_b <- tibble(
#     week_end = dates_b,
#     corrected_value = seq(target_value, first_obs_value, length.out = length(dates_b)),
#     point_type = "interpolated"
#   ) |>
#     filter(week_end > target_date)
# } else {
#   seg_b <- tibble(
#     week_end = as.Date(character()),
#     corrected_value = numeric(),
#     point_type = character()
#   )
# }
# 
# df_imputed <- bind_rows(seg_a, seg_b) |>
#   filter(week_end < first_obs_date) |>
#   distinct(week_end, .keep_all = TRUE) |>
#   arrange(week_end)
# 
# df_true <- df_raw |>
#   filter(week_end >= first_obs_date) |>
#   transmute(
#     week_end,
#     corrected_value = value,
#     point_type = "observed"
#   )
# 
# # ============================================================
# # Redistribute each post-flat jump backward over the flat stretch
# # ============================================================
# redistribute_postflat_jumps <- function(df, backlog_end) {
#   df <- df |>
#     arrange(week_end)
#   
#   n <- nrow(df)
#   i <- 2
#   
#   while (i <= n) {
#     if (df$week_end[i] > backlog_end && df$corrected_value[i] == df$corrected_value[i - 1]) {
#       plateau_start <- i - 1
#       j <- i
#       
#       while (j <= n && df$corrected_value[j] == df$corrected_value[plateau_start]) {
#         j <- j + 1
#       }
#       
#       if (j <= n) {
#         start_val <- df$corrected_value[plateau_start]
#         end_val   <- df$corrected_value[j]
#         
#         if (end_val > start_val) {
#           idx <- plateau_start:j
#           df$corrected_value[idx] <- seq(start_val, end_val, length.out = length(idx))
#           df$point_type[idx[idx < j]] <- "interpolated"
#         }
#         
#         i <- j
#       } else {
#         i <- n + 1
#       }
#     } else {
#       i <- i + 1
#     }
#   }
#   
#   df
# }
# 
# 
# # Process actual data - Ensure cum_prop is created here
# nl_actual <- bind_rows(df_imputed, df_true) |>
#   arrange(week_end) |>
#   redistribute_postflat_jumps(backlog_end = backlog_end) |>
#   mutate(
#     # Create the column specifically for the stars in Plot B
#     cum_prop = corrected_value / NL_POP,
#     # Create the daily-normalized value for the stars in Plot A
#     daily_inc_star = (corrected_value - lag(corrected_value, default = 0)) / 7,
#     weekly_inc_star = corrected_value - lag(corrected_value, default = 0)
#   )
# 
# # 2. DAILY INTERPOLATION
# daily_dates <- seq(min(nl_actual$week_end), max(nl_actual$week_end), by = "1 day")
# fit_func    <- splinefun(x = nl_actual$week_end, y = nl_actual$corrected_value, method = "hyman")
# 
# daily_base <- data.frame(date = daily_dates) |>
#   mutate(daily_inc = pmax(0, c(0, diff(fit_func(date)))))
# 
# # 3. SCENARIO GENERATION (With Continuous Pace Correction)
# project_scenario <- function(df, rate_mult, date_offset, label) {
#   obs <- df |>
#     mutate(adj_date = date + days(date_offset),
#            adj_inc  = daily_inc * rate_mult)
#   
#   final_pace <- mean(tail(obs$adj_inc, 7))
#   future_dates <- seq(max(obs$adj_date) + 1, STUDY_END, by = "1 day")
#   future_df <- data.frame(adj_date = future_dates, adj_inc = final_pace)
#   
#   bind_rows(obs |> select(adj_date, adj_inc), future_df) |>
#     arrange(adj_date) |>
#     mutate(
#       cum_val = cumsum(adj_inc),
#       adj_inc = ifelse(lag(cum_val, default = 0) >= NL_POP, 0, adj_inc),
#       cum_prop = pmin(1, cumsum(adj_inc) / NL_POP),
#       scenario = label
#     ) |>
#     filter(adj_date <= STUDY_END)
# }
# 
# scenarios_df <- bind_rows(
#   project_scenario(daily_base, 1, 0, "Actual Boosters"),
#   project_scenario(daily_base, 1, -35, "Early Boosters (-35 Days)"),
#   project_scenario(daily_base, 1.5, 0, "150% Faster Rate"),
#   project_scenario(daily_base, 0.66, 0, "150% Slower Rate")
# )
# 
# # WEEKLY RATE FROM DAILY SCENARIOS
# scenarios_df <- scenarios_df |>
#   mutate(week_end = ceiling_date(adj_date, unit = "week", week_start = 7) - days(1)) |>
#   group_by(scenario, week_end) |>
#   summarise(
#     weekly_inc = sum(adj_inc, na.rm = TRUE),
#     cum_prop = max(cum_prop, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Weekly actual stars up to cutoff
# actual_weekly_stars <- nl_actual |>
#   #filter(week_end <= CUTOFF_DATE) |>
#   mutate(
#     weekly_inc_star = corrected_value - lag(corrected_value, default = 0)
#   )
# 
# # 4. PUBLICATION PLOTTING
# theme_pub <- theme_clean() +
#     theme(
#       axis.text.x = element_text(size = 16, hjust = 1, angle = 45),
#       axis.text.y = element_text(size = 20),
#       axis.title.y = element_text(size = 22, color = "black"),
#       legend.position = "bottom",
#       legend.justification = "center",
#       legend.direction = "horizontal",
#       legend.box = "horizontal",
#       legend.box.just = "center",
#       legend.title = element_blank(),
#       legend.text = element_text(size = 20),
#       legend.background = element_blank(),
#       legend.key.width = unit(1.8, "cm"),
#       legend.spacing.x = unit(0.35, "cm"),
#       legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
#       plot.title = element_text(size = 24, face = "plain", hjust = 0.5),
#       plot.background = element_blank()
#     ) 
# 
# # Filter actual points for stars
# actual_stars <- nl_actual |> filter(week_end <= CUTOFF_DATE)
# 
# # Plot A: Daily Incidence
# p1 <- ggplot(scenarios_df, aes(x = adj_date, y = adj_inc, color = scenario)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.5) +
#   geom_point(data = actual_stars, aes(x = week_end, y = daily_inc_star), color = "black", shape = 8, size = 3, inherit.aes = FALSE) +
#   scale_x_date(
#     limits = c(plot_start, plot_end),
#     date_breaks = "month",
#     date_labels = "%Y %m",
#     expand = c(0, 0)
#   ) +
#   scale_color_brewer(palette = "Set1") +
#   labs(title = "Daily Delivery Rate Comparison", tag = "A", y = "Daily Increase (n)", x = NULL) +
#   theme_pub +
#   guides(colour = guide_legend(nrow = 1, byrow = TRUE))
# 
# # Plot B: Cumulative Coverage
# p2 <- ggplot(scenarios_df, aes(x = adj_date, y = cum_prop, color = scenario)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.5) +
#   geom_point(data = actual_stars, aes(x = week_end, y = cum_prop), color = "black", shape = 8, size = 4, inherit.aes = FALSE) +
#   geom_vline(xintercept = as.Date("2021-11-08"), color = "maroon", linetype = 2, linewidth = 2) +
#   scale_x_date(
#     #limits = c(plot_start, plot_end),
#     date_breaks = "2 weeks",
#     date_labels = "%b %Y",
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(labels = percent_format(), limits = c(0, 0.9)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(title = "Proportions of booster vaccination coverage in NL", tag = "B",y = "Proportion of NL eligible population", x = NULL) +
#   theme_pub +
#   guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
#   labs(tag = "B") +
#   theme(
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)   # top-left inside panel
#   )
# 
# # Plot C: Weekly Incidence
# p3 <- ggplot(scenarios_df, aes(x = week_end, y = weekly_inc, color = scenario)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.8) +
#   geom_point(
#     data = actual_weekly_stars,
#     aes(x = week_end, y = weekly_inc_star),
#     color = "black",
#     shape = 8,
#     size = 4,
#     inherit.aes = FALSE
#   ) +
#   scale_x_date(
#     limits = c(plot_start, plot_end),
#     date_breaks = "2 weeks",
#     date_labels = "%b %Y",
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(labels = comma) +
#   scale_color_brewer(palette = "Set1") +
#   labs(
#     title = "Weekly Delivery Rate Comparison",
#     tag = "C",
#     y = "Weekly increase (n)",
#     x = NULL
#   ) +
#   theme_pub +
#   guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
#   theme(
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# # Save figures
# png("../figures/daily_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
# print(p1)
# dev.off()
# 
# png("../figures/cum_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
# print(p2)
# dev.off()
# 
# png("../figures/weekly_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
# print(p3)
# dev.off()


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
NL_POP <- 510550
CA_POP <- 36991981
BOOSTER_ELIG_POP <- 426198
booster_offset <- 107
booster_start <- as.Date("2021-11-08")
backlog_end   <- as.Date("2022-01-11")
plot_start    <- booster_start
plot_end      <- as.Date("2022-05-22")
init_date     <- as.Date("2021-07-24")

NL_POP <- 510550
CUTOFF_DATE <- as.Date("2022-02-07")
STUDY_END <- as.Date("2022-12-15")

target_cumulative <- 106730

clip01 <- function(x) pmin(pmax(x, 0), 1)

# ============================================================
# Load raw data
# ============================================================
dat <- csvRead()

# ============================================================
# NL booster raw data
# ============================================================
df_raw <- dat |>
  filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional), 0, numtotal_2nd_additional),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  select(week_end, numtotal_additional, numtotal_fully) |>
  pivot_longer(-week_end, names_to = "type", values_to = "value") |>
  filter(type == "numtotal_additional") |>
  filter(week_end >= as.Date("2021-11-06") & week_end <= as.Date("2022-06-19")) |>
  arrange(week_end)

# ============================================================
# NL and Canada fully vaccinated counts
# ============================================================
nl_fully <- dat |>
  filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  transmute(
    week_end,
    fully_count = as.numeric(numtotal_fully)
  )

ca_fully <- dat |>
  filter(grepl("^canada$", prename, ignore.case = TRUE)) |>
  mutate(week_end = as.Date(week_end)) |>
  transmute(
    week_end,
    fully_count = as.numeric(numtotal_fully)
  )

# ============================================================
# Correct actual cumulative booster series
# ============================================================
eligible_pop <- BOOSTER_ELIG_POP
target_date  <- as.Date("2022-01-11")
target_value <- target_cumulative
shape_power  <- 1.8

first_obs <- df_raw |>
  filter(!is.na(value)) |>
  arrange(week_end) |>
  slice(1)

first_obs_date  <- first_obs$week_end
first_obs_value <- first_obs$value

dates_a <- seq.Date(booster_start, target_date, by = "week")
t_idx <- seq(0, 1, length.out = length(dates_a))

seg_a <- tibble(
  week_end = dates_a,
  corrected_value = target_value * (t_idx ^ shape_power),
  point_type = "interpolated"
)

if (first_obs_date > target_date) {
  dates_b <- seq.Date(target_date, first_obs_date, by = "week")
  seg_b <- tibble(
    week_end = dates_b,
    corrected_value = seq(target_value, first_obs_value, length.out = length(dates_b)),
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
# Redistribute each post-flat jump backward over the flat stretch
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
# Process actual data
# ============================================================
nl_actual <- bind_rows(df_imputed, df_true) |>
  arrange(week_end) |>
  redistribute_postflat_jumps(backlog_end = backlog_end) |>
  mutate(
    cum_prop = corrected_value / NL_POP,
    daily_inc_star = (corrected_value - lag(corrected_value, default = 0)) / 7,
    weekly_inc_star = corrected_value - lag(corrected_value, default = 0)
  )

# ============================================================
# Daily interpolation: build baseline cumulative + daily series
# ============================================================
daily_dates <- seq(min(nl_actual$week_end), STUDY_END, by = "1 day")

fit_func <- splinefun(
  x = nl_actual$week_end,
  y = nl_actual$corrected_value,
  method = "hyman"
)

daily_base <- tibble(date = daily_dates) |>
  mutate(
    cum_value_raw = fit_func(date),
    cum_value_raw = pmax(cum_value_raw, 0)
  )

# freeze observed cumulative after last observed date
last_obs_date <- max(nl_actual$week_end)
last_obs_cum  <- nl_actual |>
  filter(week_end == last_obs_date) |>
  pull(corrected_value)

daily_base <- daily_base |>
  mutate(
    cum_value = if_else(date <= last_obs_date, cum_value_raw, last_obs_cum),
    cum_value = pmin(cum_value, eligible_pop),
    cum_value = cummax(cum_value),                 # enforce monotone cumulative
    daily_inc = c(0, diff(cum_value)),
    daily_inc = pmax(daily_inc, 0),
    week_end = ceiling_date(date, unit = "week", week_start = 7) - days(1)
  ) |>
  group_by(week_end) |>
  mutate(
    weekly_rate = sum(daily_inc, na.rm = TRUE)
  ) |>
  ungroup()

# ============================================================
# Helper: time-warp cumulative curve, then derive daily/weekly
# speed > 1  => faster
# speed < 1  => slower
# shift_days < 0 => earlier
# ============================================================
make_rollout_scenario <- function(df, speed = 1, shift_days = 0, label) {
  base_dates <- df$date
  base_num   <- as.numeric(base_dates)
  base_cum   <- df$cum_value
  
  # interpolation on cumulative curve
  cum_fun <- approxfun(
    x = base_num,
    y = base_cum,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  t0 <- min(base_num)
  
  # for each scenario date t, evaluate baseline cumulative at:
  # baseline_time = t0 + speed * (t - shift_days - t0)
  scen_num <- as.numeric(base_dates)
  mapped_t <- t0 + speed * ((scen_num - shift_days) - t0)
  
  scen_cum <- cum_fun(mapped_t)
  scen_cum <- pmax(scen_cum, 0)
  scen_cum <- pmin(scen_cum, eligible_pop)
  scen_cum <- cummax(scen_cum)
  
  tibble(
    adj_date = base_dates,
    cum_val  = scen_cum
  ) |>
    mutate(
      adj_inc = c(0, diff(cum_val)),
      adj_inc = pmax(adj_inc, 0),
      cum_prop = clip01(cum_val / eligible_pop),
      week_end = ceiling_date(adj_date, unit = "week", week_start = 7) - days(1),
      scenario = label
    ) |>
    group_by(scenario, week_end) |>
    mutate(
      weekly_rate = sum(adj_inc, na.rm = TRUE)
    ) |>
    ungroup()
}

# ============================================================
# Scenario generation: scientifically consistent versions
# ============================================================
scenarios_df <- bind_rows(
  make_rollout_scenario(daily_base, speed = 1.00, shift_days =   0, label = "Actual Boosters"),
  make_rollout_scenario(daily_base, speed = 1.00, shift_days = -35, label = "Early Boosters (-35 Days)"),
  make_rollout_scenario(daily_base, speed = 1.50, shift_days =   0, label = "150% Faster Rate"),
  make_rollout_scenario(daily_base, speed = 1/1.5, shift_days = 0, label = "150% Slower Rate")
)

# Weekly dataframe with one row per scenario-week
weekly_scenarios_df <- scenarios_df |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = max(weekly_rate, na.rm = TRUE),
    cum_prop = max(cum_prop, na.rm = TRUE),
    cum_val = max(cum_val, na.rm = TRUE),
    .groups = "drop"
  )
# ============================================================
# Publication plotting theme
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
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    plot.title = element_text(size = 24, face = "plain", hjust = 0.5),
    plot.background = element_blank()
  )

# Filter actual points for stars
actual_stars <- nl_actual |>
  filter(week_end <= CUTOFF_DATE)

# ============================================================
# Plot A: Daily Incidence
# ============================================================
p1 <- ggplot(scenarios_df, aes(x = adj_date, y = adj_inc, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_point(
    data = actual_stars,
    aes(x = week_end, y = daily_inc_star),
    color = "black",
    shape = 8,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_x_date(
    #limits = c(plot_start, plot_end),
    date_breaks = "month",
    date_labels = "%Y %m",
    expand = c(0, 0)
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Daily Delivery Rate Comparison",
    tag = "A",
    y = "Daily Increase (n)",
    x = NULL
  ) +
  theme_pub +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

# ============================================================
# Plot B: Cumulative Coverage
# ============================================================
p2 <- ggplot(scenarios_df, aes(x = adj_date, y = cum_prop, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_point(
    data = actual_stars,
    aes(x = week_end, y = cum_prop),
    color = "black",
    shape = 8,
    size = 4,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = as.Date("2021-11-08"), color = "maroon", linetype = 2, linewidth = 2) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.9)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Proportions of booster vaccination coverage in NL",
    tag = "B",
    y = "Proportion of NL eligible population",
    x = NULL
  ) +
  theme_pub +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# Plot C: Weekly Rate
# ============================================================
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
  scale_x_date(
    #limits = c(plot_start, plot_end),
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Weekly Delivery Rate Comparison",
    tag = "C",
    y = "Weekly increase (n)",
    x = NULL
  ) +
  theme_pub +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# Save figures
# ============================================================
png("../figures/daily_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
print(p1)
dev.off()

png("../figures/cum_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
print(p2)
dev.off()

png("../figures/weekly_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
print(p3)
dev.off()

# use 107 days instead of offset = 150 becos we're starting from nov 8, 2021 for boosters
# day0 <- as.Date("2021-11-08") - lubridate::days(booster_offset)
# booster_counts_df <- (booster_counts_df |>
#                         select(week_end, scenario, daily_booster_rate) |>
#                         mutate(days = as.numeric(week_end - day0),  daily_rate = daily_booster_rate)|>
#                         select(days, week_end, scenario, daily_rate) |>
#                         filter(between(days,0,202 + booster_offset))
# )

weekly_scenarios_df <- weekly_scenarios_df |>
  filter(week_end <= plot_end)

print(weekly_scenarios_df)

write.csv(weekly_scenarios_df, "../data/booster_counts_scenarios.csv",row.names = FALSE)

saveVars(weekly_scenarios_df, booster_offset)





#
# ============================================================
# Scenario builders
# ============================================================
# warp_booster_curve <- function(df,
#                                value_col = "booster_prop_corr",
#                                t_ref = booster_start,
#                                time_stretch = 1,
#                                shift_days = 0,
#                                asymptote_mult = 1) {
#   x <- as.numeric(df$week_end)
#   y <- df[[value_col]] * asymptote_mult
#   
#   x_src <- as.numeric(t_ref) + (x - as.numeric(t_ref)) * time_stretch - shift_days
#   
#   y_new <- approx(
#     x = x_src,
#     y = y,
#     xout = x,
#     method = "linear",
#     rule = 2,
#     ties = mean
#   )$y
#   
#   tibble(
#     week_end = df$week_end,
#     booster_prop = clip01(y_new)
#   )
# }
# 
# make_lower_initial_ratio <- function(nl_df, ca_df, init_date, nl_pop, ca_pop) {
#   nl_init <- nl_df |>
#     filter(week_end == init_date) |>
#     pull(fully_count)
#   
#   ca_init <- ca_df |>
#     filter(week_end == init_date) |>
#     pull(fully_count)
#   
#   if (length(nl_init) != 1) stop("NL fully_count missing at init_date.")
#   if (length(ca_init) != 1) stop("Canada fully_count missing at init_date.")
#   
#   ca_rescaled <- ca_init * (nl_pop / ca_pop)
#   ratio <- min(nl_init, ca_rescaled) / nl_init
#   clip01(ratio)
# }
# 
# build_scenario_df <- function(base_df, scenario_name, point_type = "counterfactual") {
#   base_df |>
#     arrange(week_end) |>
#     mutate(
#       scenario = scenario_name,
#       point_type = point_type,
#       booster_count = round(booster_prop * BOOSTER_ELIG_POP)
#     ) |>
#     mutate(
#       weekly_new_boosters = booster_count - lag(booster_count, default = 0),
#       weekly_new_boosters = pmax(weekly_new_boosters, 0),
#       daily_booster_rate = weekly_new_boosters / 7
#     )
# }
# 
# lower_ratio <- make_lower_initial_ratio(
#   nl_df = nl_fully,
#   ca_df = ca_fully,
#   init_date = init_date,
#   nl_pop = NL_POP,
#   ca_pop = CA_POP
# )
# 
# lower_ratio_plot <- min(lower_ratio, 0.90)
# 
# actual_scen <- nl_actual |>
#   transmute(
#     week_end,
#     booster_prop = booster_prop_corr,
#     point_type = point_type
#   ) |>
#   build_scenario_df("Actual boosters", point_type = point_type)
# 
# slower_scen <- warp_booster_curve(
#   df = nl_actual,
#   value_col = "booster_prop_corr",
#   t_ref = booster_start,
#   time_stretch = 1.5,
#   shift_days = 0,
#   asymptote_mult = 1
# ) |>
#   build_scenario_df("Slower rate")
# 
# faster_scen <- warp_booster_curve(
#   df = nl_actual,
#   value_col = "booster_prop_corr",
#   t_ref = booster_start,
#   time_stretch = 1 / 1.5,
#   shift_days = 0,
#   asymptote_mult = 1
# ) |>
#   build_scenario_df("Faster rate")
# 
# earlier_scen <- warp_booster_curve(
#   df = nl_actual,
#   value_col = "booster_prop_corr",
#   t_ref = booster_start,
#   time_stretch = 1,
#   shift_days = 35,
#   asymptote_mult = 1
# ) |>
#   build_scenario_df("Earlier boosters (-35 days)")
# 
# lower_initial_scen <- warp_booster_curve(
#   df = nl_actual,
#   value_col = "booster_prop_corr",
#   t_ref = booster_start,
#   time_stretch = 1,
#   shift_days = 0,
#   asymptote_mult = lower_ratio_plot
# ) |>
#   build_scenario_df("Lower initial 2-dose")
# 
# plot_df <- bind_rows(
#   actual_scen,
#   slower_scen,
#   faster_scen,
#   earlier_scen,
#   lower_initial_scen
# ) |>
#   filter(week_end >= plot_start, week_end <= plot_end)
# 
# scenario_order <- c(
#   "Slower rate",
#   "Earlier boosters (-35 days)",
#   "Faster rate",
#   "Lower initial 2-dose",
#   "Actual boosters"
# )
# 
# plot_df <- plot_df |>
#   mutate(
#     scenario = factor(scenario, levels = scenario_order)
#   ) |>
#   arrange(scenario, week_end)
# 
# # ============================================================
# # Styling
# # ============================================================
# scenario_cols <- c(
#   "Slower rate"                 = "magenta",
#   "Earlier boosters (-35 days)" = "blue",
#   "Faster rate"                 = "goldenrod",
#   "Lower initial 2-dose"        = "#000000",
#   "Actual boosters"             = "red"
# )
# 
# shape_vals <- c(
#   "observed"       = 16,
#   "interpolated"   = 8,
#   "counterfactual" = 16
# )
# 
# actual_df     <- plot_df |> filter(scenario == "Actual boosters")
# actual_obs    <- actual_df |> filter(point_type == "observed")
# actual_interp <- actual_df |> filter(point_type == "interpolated")
# cf_df         <- plot_df |> filter(scenario != "Actual boosters")
# 
# # bridge between last interpolated backlog point and first observed actual point
# actual_bridge <- bind_rows(
#   actual_interp |> arrange(week_end) |> slice_tail(n = 1),
#   actual_obs    |> arrange(week_end) |> slice_head(n = 1)
# )
# 
# # ============================================================
# # Panel A: cumulative counts
# # ============================================================
# p_counts <- ggplot() +
#   geom_line(
#     data = cf_df,
#     aes(week_end, booster_count, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = cf_df,
#     aes(week_end, booster_count, colour = scenario),
#     size = 3
#   ) +
#   geom_line(
#     data = actual_obs,
#     aes(week_end, booster_count, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = actual_obs,
#     aes(week_end, booster_count, colour = scenario, shape = point_type),
#     size = 4
#   ) +
#   geom_line(
#     data = actual_interp,
#     aes(week_end, booster_count, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dotdash"
#   ) +
#   geom_line(
#     data = actual_bridge,
#     aes(week_end, booster_count, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dashed"
#   ) +
#   geom_point(
#     data = actual_interp,
#     aes(week_end, booster_count, colour = scenario, shape = point_type),
#     size = 5
#   ) +
#   scale_color_manual(values = scenario_cols, name = NULL, drop = FALSE) +
#   scale_shape_manual(
#     values = shape_vals,
#     breaks = c("observed", "interpolated"),
#     labels = c("Observed data", "Interpolated backlog"),
#     name = NULL
#   ) +
#   scale_x_date(
#     limits = c(plot_start, plot_end),
#     date_breaks = "2 weeks",
#     date_labels = "%b %d",
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(labels = comma) +
#   labs(
#     title = "Cumulative booster vaccination coverage in NL under different scenarios",
#     x = NULL,
#     y = "Cumulative count",
#     colour = NULL
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 20, hjust = 0.5),
#     axis.text.y = element_text(size = 20),
#     axis.title.y = element_text(size = 22, color = "black"),
#     legend.position = "bottom",
#     legend.justification = "center",
#     legend.direction = "horizontal",
#     legend.box = "vertical",
#     legend.box.just = "center",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 20),
#     legend.background = element_blank(),
#     legend.key.width = unit(1.8, "cm"),
#     legend.spacing.x = unit(0.35, "cm"),
#     legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
#     plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
#     plot.background = element_blank()
#   ) +
#   guides(
#     colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
#     shape = "none"
#     #shape  = guide_legend(nrow = 1, byrow = TRUE, order = 2)
#   ) +
#   labs(tag = "A") +
#   theme(
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# png("../figures/cum_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white", type = "cairo")
# print(p_counts)
# dev.off()
# 
# # ============================================================
# # Panel B: daily rates
# # ============================================================
# p_daily <- ggplot() +
#   geom_line(
#     data = cf_df,
#     aes(week_end, daily_booster_rate, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = cf_df,
#     aes(week_end, daily_booster_rate, colour = scenario),
#     size = 3
#   ) +
#   geom_line(
#     data = actual_obs,
#     aes(week_end, daily_booster_rate, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = actual_obs,
#     aes(week_end, daily_booster_rate, colour = scenario, shape = point_type),
#     size = 3
#   ) +
#   geom_line(
#     data = actual_interp,
#     aes(week_end, daily_booster_rate, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dashed"
#   ) +
#   geom_line(
#     data = actual_bridge,
#     aes(week_end, daily_booster_rate, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dashed"
#   ) +
#   geom_point(
#     data = actual_interp,
#     aes(week_end, daily_booster_rate, colour = scenario, shape = point_type),
#     size = 3
#   ) +
#   scale_color_manual(values = scenario_cols, name = NULL, drop = FALSE) +
#   scale_shape_manual(
#     values = shape_vals,
#     breaks = c("observed", "interpolated"),
#     labels = c("Observed data", "Interpolated backlog"),
#     name = NULL
#   ) +
#   scale_x_date(
#     limits = c(plot_start, plot_end),
#     date_breaks = "2 weeks",
#     date_labels = "%b %d",
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     title = "Daily booster vaccination rate in NL under different scenarios",
#     x = NULL,
#     y = "Daily rate",
#     colour = NULL
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 20, hjust = 0.5),
#     axis.text.y = element_text(size = 20),
#     axis.title.y = element_text(size = 22, color = "black"),
#     legend.position = "bottom",
#     legend.justification = "center",
#     legend.direction = "horizontal",
#     legend.box = "vertical",
#     legend.box.just = "center",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 20),
#     legend.background = element_blank(),
#     legend.key.width = unit(1.8, "cm"),
#     legend.spacing.x = unit(0.35, "cm"),
#     legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
#     plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
#     plot.background = element_blank()
#   ) +
#   guides(
#     colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
#     shape = "none"
#     #shape  = guide_legend(nrow = 1, byrow = TRUE, order = 2)
#   ) +
#   labs(tag = "B") +
#   theme(
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# png("../figures/daily_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white", type = "cairo")
# print(p_daily)
# dev.off()
# 
# y_top <- ceiling(max(cf_df$booster_prop, na.rm = TRUE) / 0.05) * 0.05
# y_top <- max(y_top, 0.75)
# 
# # ============================================================
# # Panel C: cumulative proportion
# # ============================================================
# p_prop <- ggplot() +
#   geom_line(
#     data = cf_df,
#     aes(week_end, booster_prop, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = cf_df,
#     aes(week_end, booster_prop, colour = scenario),
#     size = 3
#   ) +
#   geom_line(
#     data = actual_obs,
#     aes(week_end, booster_prop, colour = scenario, group = scenario),
#     linewidth = 1.5
#   ) +
#   geom_point(
#     data = actual_obs,
#     aes(week_end, booster_prop, colour = scenario, shape = point_type),
#     size = 3
#   ) +
#   geom_line(
#     data = actual_interp,
#     aes(week_end, booster_prop, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dashed"
#   ) +
#   geom_line(
#     data = actual_bridge,
#     aes(week_end, booster_prop, colour = scenario, group = scenario),
#     linewidth = 1.5,
#     linetype = "dashed"
#   ) +
#   geom_point(
#     data = actual_interp,
#     aes(week_end, booster_prop, colour = scenario, shape = point_type),
#     size = 6
#   ) +
#   scale_color_manual(values = scenario_cols, name = NULL, drop = FALSE) +
#   scale_shape_manual(
#     values = shape_vals,
#     breaks = c("observed", "interpolated"),
#     labels = c("Observed data", "Interpolated backlog"),
#     name = NULL
#   ) +
#   scale_x_date(
#     limits = c(plot_start, plot_end),
#     date_breaks = "2 weeks",
#     date_labels = "%b %d",
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(
#     breaks = seq(0, y_top, by = 0.10),
#     labels = percent_format(accuracy = 1),
#     limits = c(0, y_top)
#   ) +
#   labs(
#     title = "Proportions of booster vaccination coverage in NL under different scenarios",
#     x = NULL,
#     y = "Proportion of booster-eligible population (12+)",
#     colour = NULL
#   ) +
#   theme_clean() +
#   theme(
#     axis.text.x = element_text(size = 22, hjust = 0.5),
#     axis.text.y = element_text(size = 25),
#     axis.title.y = element_text(size = 25, color = "black"),
#     legend.position = "bottom",
#     legend.justification = "center",
#     legend.direction = "horizontal",
#     legend.box = "vertical",
#     legend.box.just = "center",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 25),
#     legend.background = element_blank(),
#     legend.key.width = unit(1.8, "cm"),
#     legend.spacing.x = unit(0.35, "cm"),
#     legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
#     plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
#     plot.background = element_blank()
#   ) +
#   guides(
#     colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
#     shape = "none"
#     #shape  = guide_legend(nrow = 1, byrow = TRUE, order = 2)
#   ) +
#   labs(tag = "B") +
#   theme(
#     plot.tag = element_text(size = 35, face = "bold", colour = "black"),
#     plot.tag.position = c(0.01, 0.98)
#   )
# 
# png("../figures/vax_props.png", width = 5000, height = 2800, res = 300, bg = "white", type = "cairo")
# print(p_prop)
# dev.off()

# booster_counts_df <- plot_df |>
#   transmute(
#     week_end,
#     scenario,
#     point_type,
#     booster_prop,
#     booster_count,
#     weekly_new_boosters,
#     daily_booster_rate
#   ) |>
#   arrange(scenario, week_end)
# 
# 
# # use 107 days instead of offset = 150 becos we're starting from nov 8, 2021 for boosters
# day0 <- as.Date("2021-11-08") - lubridate::days(booster_offset)
# booster_counts_df <- (booster_counts_df |>
#                         select(week_end, scenario, daily_booster_rate) |>
#                         mutate(days = as.numeric(week_end - day0),  daily_rate = daily_booster_rate)|>
#                         select(days, week_end, scenario, daily_rate) |>
#                         filter(between(days,0,202 + booster_offset))
#                       )
# 
# write.csv(booster_counts_df, "../data/booster_counts_scenarios.csv",row.names = FALSE)
# 
# saveVars(booster_counts_df, booster_offset)
# 
