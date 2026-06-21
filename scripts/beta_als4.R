library(shellpipes)

loadEnvironments()

stepkeep <- 150
als4_constant_value  <- 0.245

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[stepkeep:242] <- als4_constant_value

print(all.equal(base_beta,new_beta))

rdsSave(new_beta)
