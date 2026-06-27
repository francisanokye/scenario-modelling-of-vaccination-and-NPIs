library(shellpipes)

## params to constrain fitting
offset <- 150
omicron_start <- as.Date("2021-12-15")
analysis_end <- as.Date("2022-05-22")
burn_in_start <- omicron_start - lubridate::days(offset)
booster_start <- as.Date("2021-11-08")
NL_POP <- 510550

saveEnvironment()
