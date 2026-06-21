library(shellpipes)

loadEnvironments()

stepkeep <- 150
new_constant_value <- 0.349

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[stepkeep:242] <- new_constant_value

print(all.equal(base_beta,new_beta))

rdsSave(new_beta)
