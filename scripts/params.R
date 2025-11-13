library(shellpipes)
set.seed(2024)

loadEnvironments()

# parameter values
# params <- c(beta = 0.1
# 	, kappa1 = 1
# 	, kappa2 = 0.91
# 	, kappa3 = 0.3
# 	, sigma = 1/3
# 	, gamma_a = 1/10
# 	, gamma_i = 1/7
# 	, mu = 0.324
# 	, zeta = 0.75
# 	, vac2 = 0.09
# 	, vac3 = 0.0189
# 	, vac2_const = 0.09
# 	, vac3_const = 0.0189
# 	, N = 510550
# 	, V2_init = V2_init
# 	, omega_1 = 0 #1/121 # waning after 4 months
# 	, omega_2 = 1/182 # waning after 6 months
# 	, omega_3 = 1/365 # waning after 12 months
#   
# )

params <- c(
  # transmission and disease dynamics
  eps = 1e-9,
  beta      = 0.1,        # base transmission rate
  sigma     = 1/3,        # latent period (3 days)
  gamma_a   = 1/10,       # asymptomatic recovery
  gamma_i   = 1/7,        # symptomatic recovery
  p_sym     = 0.324,      # probability symptomatic
  zeta      = 0.75,       # relative infectiousness of asymptomatic
  
  # susceptibility reduction (leaky protection)
  kappa1    = 1.0,        # unvaccinated susceptibility
  kappa2P   = 0.91,       # vaccinated (dose 2, protected)
  kappa2W   = 0.97,       # vaccinated (dose 2, waned)
  kappa3P   = 0.30,       # boosted (dose 3+, protected)
  kappa3W   = 0.60,       # boosted (dose 3+, waned)
  
  # vaccination & booster rates
  vac2      = 0.09,       # 2-dose scale parameter
  vac3      = 0.0189,     # booster scale parameter
  vac2_const = 0.09,      # constant for time-invariant vaccination
  vac3_const = 0.0189,    # constant for time-invariant booster
  
  # waning rates (per day)
  omega_R1  = 0,          # post-infection waning in cohort 1
  omega_R2  = 1/182,      # post-infection waning in cohort 2
  omega_R3  = 1/365,      # post-infection waning in cohort 3
  omega_V2  = 1/182,      # vaccine-protection waning within cohort 2
  omega_V3  = 1/365,      # vaccine-protection waning within cohort 3
  
  # population
  N         = 510550,
  V2_init   = V2_init     # initial 2-dose coverage
)


# ~95% prior interval
prior_range <- list(
  kappa2P = c(0.85, 0.95)
  , kappa2W = c(0.85, 0.95)
	, kappa3P = c(0.2, 0.4)
	, kappa3W = c(0.2, 0.4)

	, gamma_a = c(1/11, 1/8)
	, gamma_i = c(1/8, 1/6)
	, sigma = c(1/4, 1/2)
	, R1 = c(1, 10)
)

saveEnvironment()
