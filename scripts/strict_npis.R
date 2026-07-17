library(shellpipes)
library(ggplot2)

loadEnvironments()

als4_constant_value  <- 0.245

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[step_start:step_end] <- als4_constant_value

print(all.equal(base_beta,new_beta))

rdsSave(new_beta)
