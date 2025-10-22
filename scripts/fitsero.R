library(tidyverse)
library(ggthemes)
library(coda)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(minpack.lm)
library(xtable)
library(mgcv)
library(shellpipes)

loadEnvironments()

# load data and format date column
seroprevalence <- csvRead()
seroprevalence <- seroprevalence |>
  mutate(week_end = as.Date(week_end, format = "%Y-%m-%d")) |>
  # filter antinucleocaspid which indicates only infections
  filter((geo == "NL") & (ab_estimate == "Anti-N estimate")) |>
  rename_at("week_end", ~"date") |>
  # select only the date and pct_mean
  select(date, pct_mean) |>
  # create the count index
  mutate(day = 1:n())

# # Create a copy of data
est_SERP <- seroprevalence

# subset only the study period
est_SERP <- est_SERP[(est_SERP$date >= "2021-11-01") & (est_SERP$date < "2022-06-09"),]
rownames(est_SERP) <- NULL

# create the count index
est_SERP$day <- 1:nrow(est_SERP) * 7
rownames(est_SERP) <- NULL

# create a new dummy dataframe in the daily scale to predict
daily_seroprev <- data.frame(day = 1:(nrow(est_SERP) * 7)) 

# Create a sequence of dates
date <- seq.Date(from = as.Date("2021-11-01"), by = "day", length.out = nrow(daily_seroprev))

# Add dates
daily_seroprev$date <- date

#-----------------------------------------------------------------------------------
# Fit a Gompertz model to convert weekly estimated seroprevalence to daily estimates
#-----------------------------------------------------------------------------------

gompertz_func = function(a,b,c,K, t) {
  a + (c * exp(-exp(-b * (t - K))))
}
#
gombertz_mod <- nlsLM(pct_mean ~ gompertz_func(a,b,c, K, day), data = est_SERP,
                      start = list(a = 0.5, b = 0.05, c = 1, K = max(est_SERP$pct_mean)),
                      control = nls.lm.control(maxiter = 20000))

y_pred_gomb <- predict(gombertz_mod, newdata = daily_seroprev, type = "response")

daily_seroprev$daily_serop <- y_pred_gomb

#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Add confidence interval (CI) to daily seroprevalence estimates
#-----------------------------------------------------------------------------------
n <- nrow(daily_seroprev)
sero_mean_y <- daily_seroprev$daily_serop
sero_sd_y <- sd(daily_seroprev$daily_serop)  # Sample standard deviation

# Z-score for 95% confidence; qnorm(0.975) gives 1.96
z <- qnorm(0.975)
#
# # Calculate confidence interval
sero_ci_width <- z * (sero_sd_y / sqrt(n))
sero_lower_ci <- sero_mean_y - sero_ci_width
sero_upper_ci <- sero_mean_y + sero_ci_width
#
# Add calculated CI to daily seroprevalence estimates
daily_seroprev$daily_serop_q025 <- sero_lower_ci
daily_seroprev$daily_serop_q975 <- sero_upper_ci


A <- ggplot() +
  geom_point(data = est_SERP, aes(x = date, y = pct_mean, color = "Seroprevalence"), size = 2) +
  geom_line(data = daily_seroprev, aes(x = date, y = daily_serop, color = "Predicted"), size = 1.5, alpha = 0.95) +
  scale_color_manual(values = c("Seroprevalence" = "black", "Predicted" = "brown")) +
  labs(x = "Date", y = "Seroprevalence(%)", title = "Gompertz Model Fit Estimating Daily Seroprevalence") +
  theme_clean() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        plot.title = element_text(size = 10, face = "bold", color = "black", hjust = 0.5),
        legend.position = c(0.25, 0.75),
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 10),
        legend.margin = margin(0, 0, 0, 0),
        plot.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

print(A)

# A person may recover from infection within two to three weeks, during which detectable antibodies develop in the blood. 
# Therefore, we account for the delay in detection by adjusting the seroprevalence estimates to align antibody 
# detection with the actual timing of infections. Specifically, the estimates are shifted backward by
# the average recovery period (14 days) to better reflect the infection dynamics rather than the delayed
# immune response

daily_seroprev <- daily_seroprev |> mutate(date = date - days(14))

# select only data from Dec 15, 2021 to June 3, 2022
daily_seroprev <- daily_seroprev[(daily_seroprev$date >= "2021-12-15") & (daily_seroprev$date <= "2022-06-03"),]
daily_seroprev$day <- 1:nrow(daily_seroprev)
rownames(daily_seroprev) <- NULL

serodat <- (daily_seroprev
	|> mutate(sero_inf = daily_serop*params[["N"]]
	)
	|> transmute(time=day
		, date
		, serop = daily_serop
		, sero_inc = round(diff(c(sero_inf,NA)))
	)
	|> pivot_longer(-c(time,date), names_to = "matrix", values_to = "value")
	|> arrange(matrix, time)
	|> mutate(time = time + offset0)
)

## trying one data at a time
#serodat <- filter(serodat, matrix == "sero_inc")

rdsSave(serodat)

