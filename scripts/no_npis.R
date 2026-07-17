library(shellpipes)

loadEnvironments()

new_constant_value <- 0.3675

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[step_start:step_end] <- new_constant_value

print(all.equal(base_beta,new_beta))

rdsSave(new_beta)
