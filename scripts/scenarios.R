library(tidyverse)
library(zoo)
library(lubridate)
library(shellpipes)

loadEnvironments()

offset <- 107

analysis_start <- as.Date("2021-11-08")
analysis_end   <- as.Date("2022-05-22")
init_date      <- analysis_start - days(offset)

booster_daily_vac <- boosterframe |>
  arrange(date) |>
  filter(date >= analysis_start, date <= analysis_end) |>
  mutate(day = as.numeric(date - analysis_start))

start_date <- min(booster_daily_vac$date)
max_historical_day <- max(booster_daily_vac$day)

days_earlier     <- 35
speed_up_factor  <- 1.5
slow_down_factor <- 0.5

plot_times_early <- seq(-days_earlier, max_historical_day, by = 1)

plot_times_main <- seq(0, ceiling(max_historical_day / slow_down_factor),by = 1)

make_rate_fun <- function(day, rate) {
  approxfun(
    x = day,
    y = rate,
    method = "linear",
    rule = 1,
    yleft = 0,
    yright = 0
  )
}

early_func <- make_rate_fun(booster_daily_vac$day - days_earlier, booster_daily_vac$daily_rate)

fast_func <- make_rate_fun(booster_daily_vac$day / speed_up_factor, booster_daily_vac$daily_rate * speed_up_factor)

slow_func <- make_rate_fun(booster_daily_vac$day / slow_down_factor, booster_daily_vac$daily_rate * slow_down_factor)

results_list <- list()

#####################################################################
## Actual boosters: keep exactly as original daily data
#####################################################################

results_list[["Actual Boosters"]] <- booster_daily_vac |>
  transmute(
    date,
    Scenario = "Actual Boosters",
    Daily_Rate = daily_rate,
    Cumulative_Total = cumsum(daily_rate)
  )

#####################################################################
## Early boosters: starts before actual rollout
#####################################################################

rates_early <- early_func(plot_times_early)
rates_early[is.na(rates_early)] <- 0
rates_early <- pmax(rates_early, 0)

results_list[["Early Boosters"]] <- tibble(
  date = start_date + plot_times_early,
  Scenario = "Early Boosters",
  Daily_Rate = rates_early,
  Cumulative_Total = cumsum(rates_early)
)

#####################################################################
## Faster boosters: starts same date, finishes earlier
#####################################################################

rates_fast <- fast_func(plot_times_main)
rates_fast[is.na(rates_fast)] <- 0
rates_fast <- pmax(rates_fast, 0)

results_list[["Faster Boosters"]] <- tibble(
  date = start_date + plot_times_main,
  Scenario = "Faster Boosters",
  Daily_Rate = rates_fast,
  Cumulative_Total = cumsum(rates_fast)
)

#####################################################################
## Slower boosters: starts same date, finishes later
#####################################################################

rates_slow <- slow_func(plot_times_main)
rates_slow[is.na(rates_slow)] <- 0
rates_slow <- pmax(rates_slow, 0)

results_list[["Slower Boosters"]] <- tibble(
  date = start_date + plot_times_main,
  Scenario = "Slower Boosters",
  Daily_Rate = rates_slow,
  Cumulative_Total = cumsum(rates_slow)
)

final_plot_df <- bind_rows(results_list)

final_plot_df |>
  group_by(Scenario) |>
  summarise(
    total_doses = max(Cumulative_Total),
    .groups = "drop"
  ) |>
  print()

scen_palette <- c(
  "Actual Boosters" = "#222222",
  "Early Boosters"  = "#1F77B4",
  "Faster Boosters" = "#2CA02C",
  "Slower Boosters" = "#D62728"
)

p1_daily <- ggplot(final_plot_df, aes(date, Daily_Rate, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = scen_palette) +
  scale_x_date(date_breaks = "month", date_labels = "%b %d %Y") +
  labs(
    title = "A) Daily vaccination delivery rates",
    x = "Date",
    y = "Doses per day",
    colour = "Scenario"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p2_cum <- ggplot(final_plot_df, aes(date, Cumulative_Total, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = scen_palette) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title = "B) Cumulative vaccination doses",
    x = "Date",
    y = "Cumulative doses",
    colour = "Scenario"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

print(p1_daily)
print(p2_cum)