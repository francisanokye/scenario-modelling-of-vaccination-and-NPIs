library(zoo)
library(lubridate)
library(shellpipes)

## TODO: organize the constrain parameters into a meaningful way

## params to constrain fitting
offset <- 150
omicron_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")
burn_in_start <- omicron_start - days(offset)
booster_start <- as.Date("2021-11-08")
early_start_date <- as.Date("2021-10-04")

NL_POP <- 510550
CA_POP <- 36991981
booster_max <- 296186
booster_max <- 297000
upper_plot_time <- 300
sim_horizon = 530
new_V2_init <- NA_real_

intervention_removed <- as.Date("2022-03-14")
time_omicron_start <- as.numeric(omicron_start - burn_in_start) + 1
time_analysis_end_date	<- as.numeric(analysis_end - burn_in_start) + 1

model_outputs <- c(
  "incS", "incv2", "incv3", "inc", "S1", "V2", "V3",
  "E1", "E2", "E3", "A1", "A2", "A3", "I1", "I2", "I3",
  "R1", "R2", "R3", "foi", "beta_thing",
  "inc_leq2dose", "inc_booster",
  "cum_inc_leq2dose", "cum_inc_booster", "cum_inc_total",
  "hosp_leq2dose", "hosp_booster", "hosp_total",
  "death_leq2dose", "death_booster", "death_total",
  "cum_hosp_leq2dose", "cum_hosp_booster", "cum_hosp_total",
  "cum_death_leq2dose", "cum_death_booster", "cum_death_total",
  "cum_double_vac", "cum_booster", "boosterbox"
)

saveEnvironment()
