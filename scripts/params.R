
library(shellpipes)

set.seed(2024)

loadEnvironments()

# ------------------------------------------------------------
# Severe outcome risks per infection
# ------------------------------------------------------------
## Canadian/Ontario Omicron studies show Omicron had lower hospitalization/death risk than Delta,
# while vaccination—especially a third dose—strongly preserved protection against severe outcomes.
# Ontario data estimated third-dose VE against Omicron severe outcomes around 95%, and Public Health Ontario
# estimated thousands of hospitalizations were prevented by boosters during Omicron.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8855311/


# Hospitalization probability per infection
p_hosp1 <- 0.0186 # 1484
p_hosp2 <- 0.0071
p_hosp3 <- 0.0030

# Death probability per infection
p_death1 <- 0.0050
p_death2 <- 0.0010
p_death3 <- 0.0005
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
  mu = 0.676,

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

