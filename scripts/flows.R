library(macpan2)
library(shellpipes)

# this describes the model flow of how individuals transition across compartments and cohorts
flows <- list(foi ~ beta * (zeta * (A1 + A2 + A3) + (I1 + I2 + I3)) / N
	, mp_per_capita_flow("S1", "E1", "kappa1 * foi", "incS")
	, mp_per_capita_flow("S1", "V2", "((vac2 * S1) / (vac2 + S1)) / S1", "double_vac")
	, mp_per_capita_flow("E1", "I1", "sigma * mu", "inc_symp1")
	, mp_per_capita_flow("E1", "A1", "sigma * (1 - mu)", "inc_asymp1")
	, mp_per_capita_flow("A1", "R1", "gamma_a", "asymp_recov1")
	, mp_per_capita_flow("I1", "R1", "gamma_i", "symp_recov1")
	#, mp_per_capita_flow("R1", "S1", "omega_1", "waning1")
	
	, mp_per_capita_flow("V2", "E2", "kappa2 * foi", "incv2")
	, mp_per_capita_flow("V2", "V3", "((vac3 * V2) / (vac3 + V2)) / V2", "booster_shot")
	, mp_per_capita_flow("E2", "I2", "sigma * mu", "inc_symp2")
	, mp_per_capita_flow("E2", "A2", "sigma * (1 - mu)", "inc_asymp2")
	, mp_per_capita_flow("A2", "R2", "gamma_a", "asymp_recov2")
	, mp_per_capita_flow("I2", "R2", "gamma_i", "symp_recov2")
	#, mp_per_capita_flow("R2", "V2", "omega_2", "waning2")

	, mp_per_capita_flow("V3", "E3", "kappa3 * foi", "incv3")
	, mp_per_capita_flow("E3", "I3", "sigma * mu", "inc_symp3")
	, mp_per_capita_flow("E3", "A3", "sigma * (1 - mu)", "inc_asymp3")
	, mp_per_capita_flow("A3", "R3", "gamma_a", "asymp_recov3")
	, mp_per_capita_flow("I3", "R3", "gamma_i", "symp_recov3")
  	#, mp_per_capita_flow("R3", "V3", "omega_3", "waning3")

 	, S ~ S1 + V2 + V3
	, E ~ E1 + E2 + E3
	, A ~ A1 + A2 + A3
	, I ~ I1 + I2 + I3
	, R ~ R1 + R2 + R3

	, inc ~ incS + incv2 + incv3
	, sero_inc_total ~ inc
	, serop_total ~ R / N
	, beta_thing ~ beta
)

saveVars(flows)
