# load required libraries
library(tidyverse)
library(ggthemes)
library(coda)
library(lubridate)
library(ggplot2)
library(cowplot)
library(minpack.lm)
library(xtable)
library(mgcv)

# load and filter seroprevalence data
seroprevalence <- read.csv("raw_citf_data.csv") |>
  mutate(week_end = as.Date(week_end)) |>
  filter(geo == "NL", ab_estimate == "Anti-N estimate") |>
  rename(date = week_end) |>
  select(date, pct_mean) |>
  mutate(day = 1:n())

# subset to study period and scale day by 7
est_SERP <- seroprevalence |>
  filter(date >= as.Date("2021-11-01") & date < as.Date("2022-06-09")) |>
  mutate(day = 1:n() * 7)

# create daily prediction frame
daily_seroprev <- data.frame(day = 1:(nrow(est_SERP) * 7))
daily_seroprev$date <- seq.Date(from = as.Date("2021-11-01"), by = "day", length.out = nrow(daily_seroprev))

# define gompertz function
gompertz_func <- function(a, b, c, K, t) {
  a + (c * exp(-exp(-b * (t - K))))
}

# fit gompertz model
gombertz_mod <- nlsLM(
  pct_mean ~ gompertz_func(a, b, c, K, day),
  data = est_SERP,
  start = list(a = 0.5, b = 0.05, c = 1, K = max(est_SERP$pct_mean)),
  control = nls.lm.control(maxiter = 20000)
)

# predict daily seroprevalence
daily_seroprev$daily_serop <- predict(gombertz_mod, newdata = daily_seroprev)

# add confidence intervals
n <- nrow(daily_seroprev)
z <- qnorm(0.975)
ci_width <- z * (sd(daily_seroprev$daily_serop) / sqrt(n))
daily_seroprev <- daily_seroprev |>
  mutate(
    daily_serop_q025 = daily_serop - ci_width,
    daily_serop_q975 = daily_serop + ci_width
  )

# shift by 14 days to account for antibody delay
daily_seroprev <- daily_seroprev |>
  mutate(date = date - days(14)) |>
  filter(date >= as.Date("2021-12-15") & date <= as.Date("2022-06-03")) |>
  mutate(day = 1:n())

# plot gompertz fit
p <- ggplot() +
  geom_point(data = est_SERP, aes(x = date, y = pct_mean, color = "seroprevalence"), size = 2) +
  geom_line(data = daily_seroprev, aes(x = date, y = daily_serop, color = "predicted"), size = 1.5, alpha = 0.95) +
  scale_color_manual(values = c("seroprevalence" = "black", "predicted" = "brown")) +
  labs(x = "date", y = "seroprevalence (%)", title = "gompertz model fit estimating daily seroprevalence") +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    plot.title = element_text(size = 10, face = "bold", color = "black", hjust = 0.5),
    legend.position = c(0.25, 0.75),
    legend.title = element_blank(),
    legend.background = element_rect(color = NA),
    legend.text = element_text(size = 10),
    legend.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )

print(p)

# optional: save output
# write.csv(daily_seroprev, "./daily_se.csv", row.names = FALSE)
