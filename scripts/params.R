library(shellpipes)
set.seed(2024)

loadEnvironments()

gamma_pmf_discrete <- function(mean, sd, max_lag = 40) {
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  tau <- 0:max_lag
  
  w <- pgamma(tau + 1, shape = shape, scale = scale) -
    pgamma(tau, shape = shape, scale = scale)
  
  w / sum(w)
}

make_delay_state_defaults <- function(prefix, weights) {
  n_lag <- length(weights) - 1
  
  setNames(
    rep(0, n_lag),
    sprintf("%s_%02d", prefix, seq_len(n_lag) - 1)
  )
}

mu_val <- 0.676

# Risks per infection
p_hosp1_inf <- 0.030
p_hosp2_inf <- 0.015
p_hosp3_inf <- 0.006

p_death1_inf <- 0.0040
p_death2_inf <- 0.0015
p_death3_inf <- 0.0005

# Convert to risks conditional on symptomatic infection
p_hosp1_symp <- p_hosp1_inf / mu_val
p_hosp2_symp <- p_hosp2_inf / mu_val
p_hosp3_symp <- p_hosp3_inf / mu_val

p_death1_symp <- p_death1_inf / mu_val
p_death2_symp <- p_death2_inf / mu_val
p_death3_symp <- p_death3_inf / mu_val

w_hosp <- gamma_pmf_discrete(mean = 7, sd = 4, max_lag = 40)
w_death <- gamma_pmf_discrete(mean = 18, sd = 7, max_lag = 40)

severity_delay_states <- c(
  make_delay_state_defaults("H_leq2dose", w_hosp),
  make_delay_state_defaults("H_booster",  w_hosp),
  make_delay_state_defaults("D_leq2dose", w_death),
  make_delay_state_defaults("D_booster",  w_death)
)

params <- c(
  beta = 0.1,
  kappa1 = 1,
  kappa2 = 0.91,
  kappa3 = 0.3,
  sigma = 1 / 3,
  gamma_a = 1 / 10,
  gamma_i = 1 / 7,
  mu = mu_val,
  zeta = 0.75,
  vac2 = 0.09,
  vac3 = 0.0189,
  vac2_const = 0.09,
  vac3_const = 0.0189,
  N = 510550,
  V2_init = V2_init,
  
  double_vac = 0,
  booster_shot = 0,
  
  cum_inc_leq2dose = 0,
  cum_inc_booster = 0,
  cum_inc_total = 0,
  
  p_hosp1 = p_hosp1_symp,
  p_hosp2 = p_hosp2_symp,
  p_hosp3 = p_hosp3_symp,
  
  p_death1 = p_death1_symp,
  p_death2 = p_death2_symp,
  p_death3 = p_death3_symp,
  
  severity_delay_states
)

prior_range <- list(
  kappa2 = c(0.85, 0.95),
  kappa3 = c(0.2, 0.4),
  gamma_a = c(1 / 11, 1 / 8),
  gamma_i = c(1 / 8, 1 / 6),
  sigma = c(1 / 4, 1 / 2),
  R1 = c(1, 10)
)

saveEnvironment()