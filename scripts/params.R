library(shellpipes)
set.seed(2024)

loadEnvironments()

# parameter values
params <- c(beta = 0.1
	, kappa1 = 1
	, kappa2 = 0.91
	, kappa3 = 0.3
	, sigma = 1/3
	, gamma_a = 1/10
	, gamma_i = 1/7
	, mu = 0.676 #0.324
	, zeta = 0.75
	, vac2 = 0.09
	, vac3 = 0.0189
	, vac2_const = 0.09
	, vac3_const = 0.0189
	, N = 510550
	, V2_init = V2_init
  
)

# ~95% prior interval
prior_range <- list(kappa2 = c(0.85, 0.95)
	, kappa3 = c(0.2, 0.4)
	, gamma_a = c(1/11, 1/8)
	, gamma_i = c(1/8, 1/6)
	, sigma = c(1/4, 1/2)
	, R1 = c(1, 10)
)

saveEnvironment()

