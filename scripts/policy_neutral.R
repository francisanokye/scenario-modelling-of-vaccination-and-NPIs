library(shellpipes)

loadEnvironments()

stepkeep <- 150
new_constant_value <- 0.365

base_beta <- rdsRead()

new_beta <- base_beta

new_beta[stepkeep:239] <- new_constant_value

# this is to check if we are changing correctly
print(base_beta[148:239])
print(new_beta[148:239])


print(all.equal(base_beta,new_beta))

rdsSave(new_beta)
