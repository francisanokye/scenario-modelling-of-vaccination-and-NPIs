
library(shellpipes)

set.seed(2024)

loadEnvironments()

mu_val <- 0.676
# ------------------------------------------------------------
# Severe outcome risks per infection
# ------------------------------------------------------------
## Canadian/Ontario Omicron studies show Omicron had lower hospitalization/death risk than Delta,
# while vaccination—especially a third dose—strongly preserved protection against severe outcomes.
# Ontario data estimated third-dose VE against Omicron severe outcomes around 95%, and Public Health Ontario
# estimated thousands of hospitalizations were prevented by boosters during Omicron.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8855311/

# Because Omicron was more transmissible but less severe than Delta, the average baseline probability of hospitalization 
# for an infected individual was roughly 1.0% to 1.5%, and the baseline probability of death (IFR) was roughly 0.1% to 0.3% across 
# the general population. We obtain how many hospitalizations will occur in a specific vaccination group by multiplying the number of infections 
#  by \(0.045\) (or \(4.5\%\) risk per infection for unvaccinated, fully Vaccinated: by \(0.012\) (or \(1.2\%\) risk per infection and 
# booster vaccinated by \(0.004\) (or \(0.4\%\) risk per infection. For the number of deaths, unvaccinated: multiply infections 
# by \(0.0090\) (or \(0.90\%\) risk per infection), fully vaccinated: multiply infections by \(0.0015\) (or \(0.15\%\) risk per infection) and 
# booster vaccinated: multiply infections by \(0.0003\) (or \(0.03\%\) risk per infection.

# p_hosp1 <- 0.05
# p_hosp2 <- 0.014
# p_hosp3 <- 0.007
# 
# p_death1 <- 0.002
# p_death2 <- 0.001
# p_death3 <- 0.0005

p_hosp1 <- 0.05
p_hosp2 <- 0.015
p_hosp3 <- 0.007

p_death1 <- 0.005
p_death2 <- 0.001
p_death3 <- 0.0007
# ------------------------------------------------------------
# Model parameters
# ------------------------------------------------------------

params <- c(

  # transmission
  beta = 0.1,
  zeta = 0.75,

  # vaccine-modified susceptibility
  kappa1 = 1.00,
  kappa2 = 0.91,
  kappa3 = 0.30,

  # progression
  sigma = 1 / 3,
  gamma_a = 1 / 10,
  gamma_i = 1 / 7,
  mu = mu_val,

  # vaccination
  vac2 = 0.09,
  vac3 = 0.0189,
  vac2_const = 0.09,
  vac3_const = 0.0189,

  # population
  N = 510550,
  V2_init = V2_init,

  # ----------------------------------------------------------
  # vaccination counters
  # ----------------------------------------------------------

  double_vac = 0,
  booster_shot = 0,

  cum_double_vac = 0,
  cum_booster = 0,

  # ----------------------------------------------------------
  # infection incidence
  # ----------------------------------------------------------

  incS = 0,
  incv2 = 0,
  incv3 = 0,

  inc = 0,
  total_inc = 0,

  inc_leq2dose = 0,
  inc_booster = 0,

  # ----------------------------------------------------------
  # symptomatic incidence
  # ----------------------------------------------------------

  inc_symp1 = 0,
  inc_symp2 = 0,
  inc_symp3 = 0,

  inc_symp_leq2dose = 0,
  inc_symp_booster = 0,
  inc_symp_total = 0,

  # ----------------------------------------------------------
  # asymptomatic incidence
  # ----------------------------------------------------------

  inc_asymp1 = 0,
  inc_asymp2 = 0,
  inc_asymp3 = 0,

  inc_asymp_leq2dose = 0,
  inc_asymp_booster = 0,
  inc_asymp_total = 0,

  # ----------------------------------------------------------
  # cumulative infections
  # ----------------------------------------------------------

  cum_inc_leq2dose = 0,
  cum_inc_booster = 0,
  cum_inc_total = 0,

  cum_symp_leq2dose = 0,
  cum_symp_booster = 0,
  cum_symp_total = 0,

  # ----------------------------------------------------------
  # hospitalization incidence
  # ----------------------------------------------------------

  hosp_1dose = 0,
  hosp_2dose = 0,
  hosp_3dose = 0,

  hosp_leq2dose = 0,
  hosp_booster = 0,
  hosp_total = 0,

  # ----------------------------------------------------------
  # death incidence
  # ----------------------------------------------------------

  death_1dose = 0,
  death_2dose = 0,
  death_3dose = 0,

  death_leq2dose = 0,
  death_booster = 0,
  death_total = 0,

  # ----------------------------------------------------------
  # cumulative hospitalizations
  # ----------------------------------------------------------

  cum_hosp_1dose = 0,
  cum_hosp_2dose = 0,
  cum_hosp_3dose = 0,

  cum_hosp_leq2dose = 0,
  cum_hosp_booster = 0,
  cum_hosp_total = 0,

  # ----------------------------------------------------------
  # cumulative deaths
  # ----------------------------------------------------------

  cum_death_1dose = 0,
  cum_death_2dose = 0,
  cum_death_3dose = 0,

  cum_death_leq2dose = 0,
  cum_death_booster = 0,
  cum_death_total = 0,

  # ----------------------------------------------------------
  # serology / reporting
  # ----------------------------------------------------------

  sero_inc_total = 0,
  serop_total = 0,
  newR_total = 0,
  beta_thing = 0,

  # ----------------------------------------------------------
  # severe outcome risks
  # risks are applied per infection
  # ----------------------------------------------------------

  p_hosp1 = p_hosp1,
  p_hosp2 = p_hosp2,
  p_hosp3 = p_hosp3,

  p_death1 = p_death1,
  p_death2 = p_death2,
  p_death3 = p_death3
)

# ------------------------------------------------------------
# Prior ranges
# ------------------------------------------------------------

prior_range <- list(
  kappa2 = c(0.85, 0.95),
  kappa3 = c(0.20, 0.40),
  gamma_a = c(1 / 11, 1 / 8),
  gamma_i = c(1 / 8, 1 / 6),
  sigma = c(1 / 4, 1 / 2),
  R1 = c(1, 10)
)

saveEnvironment()

