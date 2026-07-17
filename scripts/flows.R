
library(macpan2)
library(shellpipes)

flows <- list(
  foi ~ beta * (zeta * (A1 + A2 + A3) + (I1 + I2 + I3)) / N,
  
  mp_per_capita_flow("S1", "E1", "kappa1 * foi", "incS"),
  mp_per_capita_flow("S1", "V2", "((vac2 * S1) / (vac2 + S1)) / S1", "double_vac"),
  mp_per_capita_flow("E1", "I1", "sigma * mu", "inc_symp1"),
  mp_per_capita_flow("E1", "A1", "sigma * (1 - mu)", "inc_asymp1"),
  mp_per_capita_flow("A1", "R1", "gamma_a", "asymp_recov1"),
  mp_per_capita_flow("I1", "R1", "gamma_i", "symp_recov1"),
  
  mp_per_capita_flow("V2", "E2", "kappa2 * foi", "incv2"),
  mp_per_capita_flow("V2", "V3", "((vac3 * V2) / (vac3 + V2)) / V2", "booster_shot"),
  mp_per_capita_flow("E2", "I2", "sigma * mu", "inc_symp2"),
  mp_per_capita_flow("E2", "A2", "sigma * (1 - mu)", "inc_asymp2"),
  mp_per_capita_flow("A2", "R2", "gamma_a", "asymp_recov2"),
  mp_per_capita_flow("I2", "R2", "gamma_i", "symp_recov2"),
  
  mp_per_capita_flow("V3", "E3", "kappa3 * foi", "incv3"),
  mp_per_capita_flow("E3", "I3", "sigma * mu", "inc_symp3"),
  mp_per_capita_flow("E3", "A3", "sigma * (1 - mu)", "inc_asymp3"),
  mp_per_capita_flow("A3", "R3", "gamma_a", "asymp_recov3"),
  mp_per_capita_flow("I3", "R3", "gamma_i", "symp_recov3"),
  
  S ~ S1 + V2 + V3,
  E ~ E1 + E2 + E3,
  A ~ A1 + A2 + A3,
  I ~ I1 + I2 + I3,
  R ~ R1 + R2 + R3,

  p1 ~ S1/N ,
  p2 ~ V2/N ,
  p3 ~ V3/N ,

  boosterbox ~ V3 + E3 + A3 + I3 + R3,
  
  cum_double_vac ~ cum_double_vac + double_vac,
  cum_booster ~ cum_booster + booster_shot,
  
  inc ~ incS + incv2 + incv3,
  inc_leq2dose ~ incS + incv2,
  inc_booster ~ incv3,
  
  inc_symp_total ~ inc_symp1 + inc_symp2 + inc_symp3,
  inc_symp_leq2dose ~ inc_symp1 + inc_symp2,
  inc_symp_booster ~ inc_symp3,
  
  inc_asymp_total ~ inc_asymp1 + inc_asymp2 + inc_asymp3,
  inc_asymp_leq2dose ~ inc_asymp1 + inc_asymp2,
  inc_asymp_booster ~ inc_asymp3,
  
  cum_inc_leq2dose ~ cum_inc_leq2dose + inc_leq2dose,
  cum_inc_booster ~ cum_inc_booster + inc_booster,
  cum_inc_total ~ cum_inc_total + inc,
  
  cum_symp_leq2dose ~ cum_symp_leq2dose + inc_symp_leq2dose,
  cum_symp_booster ~ cum_symp_booster + inc_symp_booster,
  cum_symp_total ~ cum_symp_total + inc_symp_total,
  
  hosp_leq2dose ~ p_hosp1 * incS + p_hosp2 * incv2,
  hosp_booster ~ p_hosp3 * incv3,
  hosp_total ~ hosp_leq2dose + hosp_booster,
  
  cum_hosp_leq2dose ~ cum_hosp_leq2dose + hosp_leq2dose,
  cum_hosp_booster ~ cum_hosp_booster + hosp_booster,
  cum_hosp_total ~ cum_hosp_total + hosp_total,
  
  death_leq2dose ~ p_death1 * incS + p_death2 * incv2,
  death_booster ~ p_death3 * incv3,
  death_total ~ death_leq2dose + death_booster,
  
  cum_death_leq2dose ~ cum_death_leq2dose + death_leq2dose,
  cum_death_booster ~ cum_death_booster + death_booster,
  cum_death_total ~ cum_death_total + death_total,
  
  newR_total ~ asymp_recov1 + asymp_recov2 + asymp_recov3 + 
               symp_recov1 + symp_recov2 + symp_recov3,
  
#  sero_inc_total ~ inc,
  serop_total ~ R / N,
  beta_thing ~ beta
#  total_inc ~ R1 + R2 + R3
)

saveVars(flows)
