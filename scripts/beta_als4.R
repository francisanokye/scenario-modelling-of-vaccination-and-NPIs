library(shellpipes)
library(ggplot2)

loadEnvironments()

offset <- 150
omicron_start  <- as.Date("2021-12-15")
intervention_removed <- as.Date("2022-03-14")
burn_in_start  <- omicron_start - lubridate::days(offset)

step_start <- as.numeric(omicron_start - burn_in_start) + 1
step_end   <- as.numeric(intervention_removed - burn_in_start) + 1

als4_constant_value  <- 0.245

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[step_start:step_end] <- als4_constant_value

print(all.equal(base_beta,new_beta))

#print(ggplot() + geom_point(aes(x = seq_along(new_beta[1:length(new_beta)]), y = new_beta)))

rdsSave(new_beta)