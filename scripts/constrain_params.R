library(zoo)
library(shellpipes)

## TODO: organize the constrain parameters into a meaningful way

## params to constrain fitting
offset <- 150
omicron_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")
burn_in_start <- omicron_start - lubridate::days(offset)
booster_start <- as.Date("2021-11-08")
NL_POP <- 510550

upper_plot_time <- 300
sim_horizon = 530

intervention_removed <- as.Date("2022-03-14")
time_omicron_start <- as.numeric(omicron_start - burn_in_start) + 1
time_analysis_end_date	<- as.numeric(analysis_end - burn_in_start) + 1

saveEnvironment()
