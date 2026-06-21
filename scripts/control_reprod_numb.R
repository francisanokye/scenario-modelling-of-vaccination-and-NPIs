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
calibrator     <- readRDS("calibrate.rds")
seroprevdata   <- readRDS("seroprevdata.rds")

beta_ALS4only  <- rdsRead("beta_als4")
beta_noALS     <- rdsRead("policy_neutral")

upper_plot_time <- 300
offset0         <- 150
anchor_start    <- as.Date("2021-12-15")
anchor_end      <- as.Date("2022-05-22")
cut_date        <- as.Date("2022-03-15")

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
		, time >= offset0
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
		, date = anchor_start + time -offset0  
	 	, beta = beta_thing_value
		, beta_lwr = beta_thing_conf.low
		, beta_upr = beta_thing_conf.high
		, p1 = p1_value
		, p2 = p2_value
		, p3 = p3_value
		, susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3
#  	, susceptibility = 0.15 * kappa1 + 0.85 * kappa2 + 0 * kappa3
		, psum = p1 + p2 + p3
		, mult_t = bracket_term * susceptibility
		, Rc = beta*mult_t
		, Rc_lwr = beta_lwr * mult_t
		, Rc_upr = beta_upr * mult_t
  		) 
	|> filter(date >= anchor_start, date <= anchor_end) 
)

print(dim(baseRc))

newRc <- function(beta_vec,scenario){
	Rcdat <- (baseRc
		|> select(date,beta,Rc,Rc_lwr,Rc_upr)
		|> mutate(beta_ratio = beta_vec/beta
			, new_beta = beta_ratio* beta
			, Rc = beta_ratio * Rc
			, Rc_lwr = beta_ratio * Rc_lwr
			, Rc_upr = beta_ratio * Rc_upr
			, scenario = scenario
		)
	)
}


noALS_Rc <- newRc(beta_noALS[-c(1:(offset0-1))],scenario="noALS")
ALS4only_Rc <- newRc(beta_ALS4only[-c(1:(offset0-1))],scenario="strict")

comboRc <- (baseRc
	|> select(date,Rc,Rc_lwr,Rc_upr)
	|> mutate(scenario = "base")
	|> bind_rows(noALS_Rc)
	|> bind_rows(ALS4only_Rc)
)

print(comboRc)

scenario_levels <- c("base", "noALS", "strict")
scenario_labels <- c(
  "base"   = "Actual NPIs",
  "noALS"  = "No NPIs",
  "strict" = "Strict NPIs"
)
scenario_colors <- c(
  "base"   = "black",
  "noALS"  = "red",
  "strict" = "#00C5CD"
)

gg <- (ggplot(comboRc, aes(date, Rc, color = scenario, fill = scenario))
       + geom_ribbon(aes(ymin = Rc_lwr, ymax = Rc_upr), alpha = 0.2, colour = NA)
       + geom_line(linewidth = 2)
       + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2)
       + geom_vline(
         xintercept = cut_date,
         linetype = "dashed",
         linewidth = 3
       )
       + scale_color_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_fill_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_x_date(
         limits = c(anchor_start, anchor_end),
         date_breaks = "2 weeks",
         date_labels = "%d %b %Y",
         expand = c(0.01, 0)
       )
       + labs(
         title = "Reproduction number for counterfactual scenarios (timevarying prop)",
         y = "Reproduction number",
         x = NULL,
         tag = "C"
       )
       + theme_clean()
       + theme(
         axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
         axis.text.y = element_text(size = 30),
         axis.title.y = element_text(size = 30, color = "black"),
         legend.position = c(0.85, 0.65),
         legend.justification = c(0.5, 0.5),
         legend.direction = "vertical",
         legend.box = "vertical",
         legend.title = element_text(size = 25),
         legend.text = element_text(size = 25),
         legend.background = element_blank(),
         plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
         plot.background = element_blank(),
         plot.tag = element_text(size = 35, face = "bold"),
         plot.tag.position = c(0.01, 0.98)
       )
)

png("../figures/timevarying_props_reprod_numb.png", width = 5000, height = 2000, res = 300, bg = "white",type = "cairo")
print(gg)
dev.off()
