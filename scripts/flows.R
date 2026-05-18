library(macpan2)
library(shellpipes)

gamma_pmf_discrete <- function(mean, sd, max_lag = 40) {
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  tau <- 0:max_lag
  
  w <- pgamma(tau + 1, shape = shape, scale = scale) -
    pgamma(tau, shape = shape, scale = scale)
  
  w / sum(w)
}

make_delay_queue <- function(prefix, input_expr, weights, output_name) {
  n_lag <- length(weights) - 1
  q <- sprintf("%s_%02d", prefix, seq_len(n_lag) - 1)
  
  output_terms <- c(
    sprintf("%.16g * (%s)", weights[1], input_expr),
    sprintf("%.16g * %s", weights[-1], q)
  )
  
  output_formula <- as.formula(
    paste(output_name, "~", paste(output_terms, collapse = " + "))
  )
  
  queue_updates <- vector("list", n_lag)
  
  queue_updates[[1]] <- as.formula(
    paste(q[1], "~", input_expr)
  )
  
  if (n_lag > 1) {
    for (i in 2:n_lag) {
      queue_updates[[i]] <- as.formula(
        paste(q[i], "~", q[i - 1])
      )
    }
  }
  
  c(list(output_formula), queue_updates)
}

w_hosp <- gamma_pmf_discrete(mean = 7, sd = 4, max_lag = 40)
w_death <- gamma_pmf_discrete(mean = 18, sd = 7, max_lag = 40)

flows <- c(
  list(
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
    
    cum_double_vac ~ cum_double_vac + double_vac,
    cum_booster ~ cum_booster + booster_shot,
    
    inc ~ incS + incv2 + incv3,
    inc_leq2dose ~ incS + incv2,
    inc_booster ~ incv3,
    
    cum_inc_leq2dose ~ cum_inc_leq2dose + inc_leq2dose,
    cum_inc_booster ~ cum_inc_booster + inc_booster,
    cum_inc_total ~ cum_inc_total + inc,
    
    newR_total ~ asymp_recov1 + asymp_recov2 + asymp_recov3 +
      symp_recov1 + symp_recov2 + symp_recov3,
    
    sero_inc_total ~ inc,
    serop_total ~ R / N,
    beta_thing ~ beta,
    total_inc ~ R1 + R2 + R3,
    
    symp_leq2dose ~ inc_symp1 + inc_symp2,
    symp_booster ~ inc_symp3,
    symp_total ~ symp_leq2dose + symp_booster,
    
    hosp_from_symp_leq2dose ~ p_hosp1 * inc_symp1 + p_hosp2 * inc_symp2,
    hosp_from_symp_booster ~ p_hosp3 * inc_symp3,
    hosp_from_symp_total ~ hosp_from_symp_leq2dose + hosp_from_symp_booster,
    
    death_from_symp_leq2dose ~ p_death1 * inc_symp1 + p_death2 * inc_symp2,
    death_from_symp_booster ~ p_death3 * inc_symp3,
    death_from_symp_total ~ death_from_symp_leq2dose + death_from_symp_booster
  ),
  
  make_delay_queue(
    prefix = "H_leq2dose",
    input_expr = "hosp_from_symp_leq2dose",
    weights = w_hosp,
    output_name = "hosp_leq2dose"
  ),
  
  make_delay_queue(
    prefix = "H_booster",
    input_expr = "hosp_from_symp_booster",
    weights = w_hosp,
    output_name = "hosp_booster"
  ),
  
  make_delay_queue(
    prefix = "D_leq2dose",
    input_expr = "death_from_symp_leq2dose",
    weights = w_death,
    output_name = "death_leq2dose"
  ),
  
  make_delay_queue(
    prefix = "D_booster",
    input_expr = "death_from_symp_booster",
    weights = w_death,
    output_name = "death_booster"
  ),
  
  list(
    hosp_total ~ hosp_leq2dose + hosp_booster,
    death_total ~ death_leq2dose + death_booster
  )
)

saveVars(flows)