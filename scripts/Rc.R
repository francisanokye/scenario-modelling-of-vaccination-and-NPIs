suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(zoo)
  library(ggthemes)
  library(patchwork)
  library(shellpipes)
  library(tidyverse)
  library(macpan2)
  library(grid)
})

theme_set(theme_bw())

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ---- inputs ----
calibrator <- rdsRead()

# ---- model constants ----
kappa1 <- params[["kappa1"]]
kappa2 <- params[["kappa2"]]
kappa3 <- params[["kappa3"]]

gamma_i <- params[["gamma_i"]]
gamma_a <- params[["gamma_a"]]

mu   <- params[["mu"]]
zeta <- params[["zeta"]]

bracket_term <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)

# ---- time-varying vaccination-cohort proportions ----
baseRc <- (mp_trajectory_sd(calibrator,conf.int=TRUE, back_transform=TRUE) 
	|>filter(matrix %in% c("p1", "p2", "p3","beta_thing")
		, time >= time_omicron_start
	) 
	|> select(time, matrix, value, conf.low,conf.high) 
	|> pivot_longer(-c(time,matrix),names_to = "est", values_to = "value")
	|> transmute(NULL
		, time 
		, vartype = paste0(matrix,"_",est)
		, value
	)
	|> pivot_wider(names_from = vartype, values_from = value)
	|> transmute(NULL
		, date = burn_in_start + (time - 1)
		, beta = beta_thing_value
		, beta_lwr = beta_thing_conf.low
		, beta_upr = beta_thing_conf.high
		, p1 = p1_value
		, p2 = p2_value
		, p3 = p3_value
#		, susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3
		, susceptibility = 0.15 * kappa1 + 0.85 * kappa2 + 0 * kappa3
		, psum = p1 + p2 + p3
		, mult_t = bracket_term * susceptibility
		, Rc = beta*mult_t
		, Rc_lwr = beta_lwr * mult_t
		, Rc_upr = beta_upr * mult_t
	) 
	|> filter(date >= omicron_start, date <= analysis_end) 
)

saveVars(baseRc)

